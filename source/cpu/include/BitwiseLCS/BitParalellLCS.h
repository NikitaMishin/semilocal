
#include <bitset>
#include "../LCS/LCS.h"
#include "unordered_set"


/// Helper macros for stringification
#define TO_STRING_HELPER(X)   #X
#define TO_STRING(X)          TO_STRING_HELPER(X)
#define GCC_UNROLL_MAX_SIZE 256 /* TODO move to Cmake */
#define COMPARISON_LOOP_UNROLL 5

// Define loop unrolling depending on the compiler
#if defined(__ICC) || defined(__ICL)
#define UNROLL_LOOP(n)      _Pragma(TO_STRING(unroll (n)))
#elif defined(__clang__)
#define UNROLL_LOOP(n)      _Pragma(TO_STRING(unroll (n)))
#elif defined(__GNUC__) && !defined(__clang__)
#define UNROLL_LOOP(n)      _Pragma(TO_STRING(GCC unroll (GCC_UNROLL_MAX_SIZE)))
#elif defined(_MSC_BUILD)
#pragma message ("Microsoft Visual C++ (MSVC) detected: Loop unrolling not supported!")
#define UNROLL_LOOP(n)
#else
#warning "Unknown compiler: Loop unrolling not supported!"
#define UNROLL_LOOP(n)
#endif

// Define loop unrolling depending on the compiler
#if defined(__ICC) || defined(__ICL)
#define PRIVATE_UNROLL_LOOP(n)      _Pragma(TO_STRING(unroll (n)))
#elif defined(__clang__)
#define PRIVATE_UNROLL_LOOP(n)      _Pragma(TO_STRING(unroll (n)))
#elif defined(__GNUC__) && !defined(__clang__)
#define PRIVATE_UNROLL_LOOP(n)      _Pragma(TO_STRING(GCC unroll (n)))
#elif defined(_MSC_BUILD)
#pragma message ("Microsoft Visual C++ (MSVC) detected: Loop unrolling not supported!")
#define PRIVATE_UNROLL_LOOP(n)
#else
#warning "Unknown compiler: Loop unrolling not supported!"
#define PRIVATE_UNROLL_LOOP(n)
#endif


template<class UnsignedMachineWord>
constexpr UnsignedMachineWord getBraidMaskOnesMSB(int alignment, int bitsPerSymbol) {
    UnsignedMachineWord singleStrand = UnsignedMachineWord(1);
    UnsignedMachineWord braidOnes = UnsignedMachineWord(1);
    for (int i = 0; i < sizeof(UnsignedMachineWord) * 8 - alignment; i += bitsPerSymbol) braidOnes |= (singleStrand << i);
    return braidOnes;
}

namespace bit_parallel::lcs {

    struct AlphabetConverter {
        template<class MachineWordType>
        struct ConversionDetails {
            int initalNumSymbols;
            int bitsPerSymbol;
            int machineWordSize;


            std::optional<std::pair<int, int>> getNthPositionOfSymbol(int i) {
                auto machineWordId = i / getNumberOfEmbeddingsInWord();
                if (machineWordId >= getCompressedSizeInWords()) return std::nullopt;

                auto withinWordId = i % getNumberOfEmbeddingsInWord();
                return std::make_pair(machineWordId, withinWordId);
            }


            inline int getWordAlignmentInBits() const { return machineWordSize - (machineWordSize / bitsPerSymbol) * bitsPerSymbol; }


            inline int getNumberOfEmbeddingsInWord() const { return machineWordSize / bitsPerSymbol; }

            inline MachineWordType getAlignmentMask() const {
                auto alignmentInBits = getWordAlignmentInBits();
                if (alignmentInBits == 0) return MachineWordType(0);
                MachineWordType mask = MachineWordType((~MachineWordType(0))) - (MachineWordType((~MachineWordType(0))) >> alignmentInBits);
                return mask;
            }

            inline MachineWordType getPaddingMaskOfLastWord(bool isLSB) const {
                auto activeSymb = initalNumSymbols % getNumberOfEmbeddingsInWord();
                if (activeSymb == 0) return MachineWordType(0);
                auto paddedSymb = getNumberOfEmbeddingsInWord() - activeSymb;
                auto align = getWordAlignmentInBits();
                if (isLSB) {
                    return (MachineWordType(~MachineWordType(0)) >> align) - ((MachineWordType(1) << activeSymb * bitsPerSymbol) - 1);
                } else {
//                    std::cout<<"HUA:"<<std::bitset<16>(((MachineWordType(1) << (align + paddedSymb * bitsPerSymbol)) - 1))<<'-'<<std::bitset<16>(getAlignmentMask(isLSB))<<std::endl;
                    return ((MachineWordType(1) << (paddedSymb * bitsPerSymbol)) - 1);
                }


            }

            inline int getPaddingOfLastWordInSymbols() const {
                auto activeEmbeddingsInLastWord = initalNumSymbols % getNumberOfEmbeddingsInWord();
                if (activeEmbeddingsInLastWord == 0) return activeEmbeddingsInLastWord;
                return getNumberOfEmbeddingsInWord() - activeEmbeddingsInLastWord;
            }


            inline int getCompressedSizeInWords() const {
                return int(std::ceil(initalNumSymbols * 1.0 / getNumberOfEmbeddingsInWord()));
            }

            inline int getCompressedSizeInBytes() const {
                return getCompressedSizeInWords() * (machineWordSize / sizeof(MachineWordType));
            }

        };

        using Number = int;
        template<class Input>
        using Map = std::unordered_map<Input, Number>;
        template<class Output>
        using ReverseMap = std::unordered_map<Number, Output>;

        /**
        * Maps alphabet symbols to numbers from interval 0 to alphabetSize
        * @tparam Input  any comparable type and hashable
        * @tparam Output number type. Unsigned int/U long,size_t,short
        * @param alphabetSet
        * @return
        */
        template<class Input>
        std::pair<Map<Input>, ReverseMap<Input>>
        encodeAlphabet(const std::unordered_set<Input> &alphabetSet) {
            Map<Input> mapperForward;
            ReverseMap<Input> mapperReverse;
            Number encoder = 0;
            for (Input s: alphabetSet) {
                mapperForward.insert({s, encoder});
                mapperReverse.insert({encoder, s});
                encoder++;
            }
            return {mapperForward, mapperReverse};
        }

        template<class Input, class Output, bool IS_LSB_ENCODING, bool PACK_REVERSE>
        std::pair<Output *, ConversionDetails<Output>> packSequence(const Input *a, int size, Map<Input> mapperForward, int alphabetSize) {
            ConversionDetails<Output> conversionDetails{
                .initalNumSymbols=size,
                .bitsPerSymbol=std::max(1, int(std::ceil(log2(alphabetSize)))),
                .machineWordSize=sizeof(Output) * 8
            };

            auto bytesNeeded = conversionDetails.getCompressedSizeInBytes();
            auto bitsetArray = static_cast<Output *> (aligned_alloc(sizeof(Output), bytesNeeded));

//            auto alignment = conversionDetails.getWordAlignmentInBits();
            auto total = 0;
            auto symbolsInWord = conversionDetails.getNumberOfEmbeddingsInWord();
            auto bitsPerSymbol = conversionDetails.bitsPerSymbol;
            auto n = conversionDetails.getCompressedSizeInWords();

            for (int i = 0; i < n; ++i, total += symbolsInWord) {
                Output word = 0;
                for (int symbol = 0; symbol < symbolsInWord; symbol++) {
                    auto character = (total + symbol < size) ? mapperForward[a[i * symbolsInWord + symbol]] : Output(0);
                    if constexpr(IS_LSB_ENCODING) {
//                        std::cout<<"CHAR_LSB"<<character<<std::endl;
                        auto offset = (bitsPerSymbol * symbol);
                        word |= (total + symbol < size) ? Output(character) << offset : Output(0) << offset;
                    } else {
//                        std::cout<<"CHAR_MBS"<<character<<std::endl;

//                        if (total + symbol < size) {
//                            auto offset = (bitsPerSymbol * (symbolsInWord - 1 - symbol) + residue);
//                            word |= Output(character) << offset;
//                        } else {
//                            auto singleOffset = bitsPerSymbol;
//                            auto offset = (bitsPerSymbol * (symbolsInWord - 1 - symbol) + residue);
//                            word |= (Output(0)  << offset);
//                        }
                        auto offset = (bitsPerSymbol * (symbolsInWord - 1 - symbol) );
                        word |= Output(character) << offset;
                    }
                }
                if constexpr(!PACK_REVERSE) {
                    bitsetArray[i] = word;
                } else {
                    bitsetArray[n - 1 - i] = word;
                }

            }

            return {bitsetArray, conversionDetails};
        }


        /**
        * Encode given sequence with Input type symbols to bit array packed in Output type
         * according to alphabet_mapping and amount of bits per symbol
        * @tparam Input
        * @tparam Output
        * @param a
         * @return (packed array, amount of packed elements, size of initial sequence)
         */
        template<class Input, class Output>
        std::pair<Output *, int> packSequenceLSB(const Input *a, int size, Map<Input> mapperForward, int alphabetSize) {
            auto bitsPerSymbol = std::max(1, int(std::ceil(log2(alphabetSize))));
            auto shift = bitsPerSymbol;
            auto wordSizeInBits = sizeof(Output) * 8;
            auto symbolsInWord = int(wordSizeInBits / shift);

            auto bytesNeeded = int(std::ceil(size * 1.0 / symbolsInWord) * sizeof(Output));

            auto bitsetArray = static_cast<Output *> (aligned_alloc(sizeof(Output), bytesNeeded));
            auto n = bytesNeeded / sizeof(Output);

            auto residue = wordSizeInBits % bitsPerSymbol; //TODO

            //    fill bitset
            auto total = 0;
            for (int i = 0; i < n; ++i, total += symbolsInWord) {
                Output word = 0;
                for (int symbol = 0; symbol < symbolsInWord; symbol++) {
                    auto offset = (shift * symbol + residue);
                    if (total + symbol < size) {
                        auto character = mapperForward[a[i * symbolsInWord + symbol]];
                        word |= Output(character) << offset;
                    } else {
                        word |= Output(0) << offset;
                    }
                }
                bitsetArray[i] = word;
            }

            return {bitsetArray, n};
        }

        /*
        *
        * Encode given sequence with Input type symbols to  reversed bit array packed in Output type
        * according to alphabet_mapping and amount of bits per symbol.
        * I.e abcde -> edcba -> to bits for each symbol
        * @return (packed array, amount of packed elements)
        */
        template<class Input, class Output>
        std::pair<Output *, int> packMSBReverse(const Input *a, int size, Map<Input> &mapperForward, int alphabetSize, bool isDefaultZero = true) {

            auto bitsPerSymbol = std::max(1, int(std::ceil(log2(alphabetSize))));
            auto shift = bitsPerSymbol;
            auto wordSizeInBits = sizeof(Output) * 8;
            auto symbolsInWord = int(wordSizeInBits / shift);

            auto bytesNeeded = int(std::ceil(size * 1.0 / symbolsInWord) * sizeof(Output));

            auto bitsetArray = static_cast<Output *> (aligned_alloc(sizeof(Output), bytesNeeded));
            auto n = bytesNeeded / sizeof(Output);
            auto residue = wordSizeInBits % bitsPerSymbol;


            auto defaultMachineWord = (isDefaultZero) ? Output(0) : ~Output(0);

            // fill first guy
            for (int i = 0; i < n - 1; ++i) {
                Output word = Output(0);
                for (int symbol = 0; symbol < symbolsInWord; symbol++) {
                    auto character = mapperForward[a[(i * symbolsInWord + symbol)]];
                    auto offset = (bitsPerSymbol * (symbolsInWord - symbol - 1) + residue);
                    word |= Output(character) << offset;
                }
                bitsetArray[n - 1 - i] = word;
            }

            //    fill last
            for (int i = n - 1; i < n; ++i) {
                Output word = Output(0);
                for (int symbol = 0; (n - 1) * symbolsInWord + symbol < size; symbol++) {
                    auto character = mapperForward[a[(i * symbolsInWord + symbol)]];
                    word |= Output(character) << (bitsPerSymbol * (symbolsInWord - symbol - 1) + residue);
                }
                bitsetArray[n - 1 - i] = word;
            }
            return {bitsetArray, n};
        }


        /**
        * Decode given packed sequence in bits in numeric type Input to symbol sequence of InputType
        * according to decode and amount of bits per symbol
         */
        template<class Input, class Output>
        void unpackLSB(Input const *a, int n, int totalSymbols, const ReverseMap<Output> &decoder, Output *sequence) {
            auto alphabetSize = decoder.size();
            auto bitsPerSymbol = int(std::ceil(log2(alphabetSize)));
            auto shift = bitsPerSymbol;
            auto wordSizeInBits = sizeof(Input) * 8;
            auto symbolsInWord = int(wordSizeInBits / shift);
            auto mask = (Input(1) << bitsPerSymbol) - 1;

            auto ptrIndex = 0;
            for (int i = 0; i < n - 1; ++i) {
                auto word = a[i];
                for (int symbol = 0; symbol < symbolsInWord; symbol++) {
                    auto key = Number((word >> shift * symbol) & mask);
                    sequence[ptrIndex] = decoder[key];
                    ptrIndex++;
                }
            }
            for (int i = n - 1; i < n; ++i) {
                auto word = a[i];
                for (int symbol = 0; i * symbolsInWord + symbol < totalSymbols; symbol++) {
                    auto key = Number((word >> shift * symbol) & mask);
                    sequence[ptrIndex] = decoder[key];
                    ptrIndex++;
                }
            }
        }

        /**
         * Decode given reversed packed sequence in bits in numeric type Input to symbol vector of type Output
         * according to decode and amount of bits per symbol
         */
        template<class Input, class Output>
        void unpackMSBReverse(Input const *a, int n, int totalSymbols, const ReverseMap<Output> &decoder) {
//            TODO
        }


    };

    template<class UnderlyingMachineWord, class InputType>
    class BitParallel : public ::lcs::LCSStrategy<InputType> {
    public:
    };

    /**
     *
     * @tparam UnderlyingMachineWord  should be unsigned type
     * @tparam InputType
     */
    template<class UnderlyingMachineWord, class InputType>
    class Hyyro : public BitParallel<UnderlyingMachineWord, InputType> {
    public:
        int compute(const InputType *a, int aSize, const InputType *b, int bSize) override {
            std::unordered_set<InputType> alphabet;
            for (int i = 0; i < aSize; ++i) alphabet.insert(a[i]);
            std::unordered_map<InputType, UnderlyingMachineWord *> lookup;
            auto words = preprocess(b, bSize, alphabet, lookup);
            return hyyroMagic(a, aSize, words, lookup);
        }

    private:
        /**
         * Build a lookup table for a given string & alphabet
         */
        int preprocess(const InputType *str, int size, std::unordered_set<InputType> &alphabet, std::unordered_map<int, UnderlyingMachineWord *> &lookup) {

            auto wordSize = (sizeof(UnderlyingMachineWord) * 8);
            auto sizeWords = int(std::ceil(size * 1.0 / wordSize));
            auto mod = size % wordSize;

            if (mod == 0) mod = wordSize;

            for (auto &symbol: alphabet) {
                auto vectorB = new UnderlyingMachineWord[sizeWords];
                auto ptr = 0;

                for (int i = 0; i < sizeWords - 1; i++) {
                    auto word = UnderlyingMachineWord(0);
                    for (int j = 0; j < wordSize; j++) {
                        auto bit = UnderlyingMachineWord(str[ptr] == symbol);
                        word |= (bit << j);
                        ptr++;
                    }
                    vectorB[i] = word;
                }

                auto word = UnderlyingMachineWord(0);
                for (int j = 0; j < mod; ++j) {
                    auto bit = UnderlyingMachineWord(str[ptr] == symbol);
                    word |= (bit << j);
                    ptr++;
                }
                vectorB[sizeWords - 1] = word;
                lookup[symbol] = vectorB;
            }
            return sizeWords;
        }

        int hyyroMagic(const InputType *strA, int sizeA, int sizeBWord, std::unordered_map<InputType, UnderlyingMachineWord *> &lookup) {

            auto vectorV = new UnderlyingMachineWord[sizeBWord];
            for (int i = 0; i < sizeBWord; ++i) vectorV[i] = ~UnderlyingMachineWord(0);


            auto bitFlagSum = new UnderlyingMachineWord[sizeBWord];
            bitFlagSum[0] = UnderlyingMachineWord(0);

            // for each symbol
            for (int i = 0; i < sizeA; ++i) {
                const auto table = lookup[strA[i]];

                for (int j = 0; j < sizeBWord - 1; ++j) {
                    auto oldV = vectorV[j];
                    auto p = table[j] & oldV;
                    auto withOffset = UnderlyingMachineWord(bitFlagSum[j] + oldV + p);

                    vectorV[j] = (oldV ^ p) | withOffset;
                    bitFlagSum[j + 1] = withOffset < oldV;
                }

                for (int j = sizeBWord - 1; j < sizeBWord; ++j) {
                    auto oldV = vectorV[j];
                    auto p = table[j] & oldV;
                    auto withOffset = bitFlagSum[j] + oldV + p;
                    vectorV[j] = (oldV ^ p) | withOffset;
                }
            }

            auto score = 0;
            for (int i1 = 0; i1 < sizeBWord; ++i1) {
                //  Brian Kernighan’s Algorithm
                int counter = 0;
                UnderlyingMachineWord number = vectorV[i1];
                //  LogNumber
                while (number) {
                    number &= (number - 1);
                    counter++;
                }
                score += sizeof(UnderlyingMachineWord) * 8 - counter;
            }
            for (auto&[_, v]: lookup) delete[] v;
            delete[] vectorV;
            delete[] bitFlagSum;
            return score;
        }
    };


    template<class UnsignedMachineWord, class InputType>
    class AbstractBitSemiLocalSticky : public BitParallel<UnsignedMachineWord, InputType> {
    public:


    protected:

        /**
        * Process upper triangles of squares of size sizeof(Input)*8 that lies on specific antidiagonal
        *  First version
        */
        inline void
        loopUpperHalfBinary(int lowerBound, int upperBound, int shift, int lEdge, int tEdge,
                            UnsignedMachineWord activeBits, UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands, UnsignedMachineWord *aReverse,
                            UnsignedMachineWord *b) {

            for (int j = lowerBound; j < upperBound; ++j) {

                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                UnsignedMachineWord lStrandCap = lStrand >> shift;
                auto symbolA = aReverse[lEdge + j];
                auto symbolB = b[tEdge + j];
                UnsignedMachineWord cond = activeBits & (((~(symbolA >> shift)) ^ symbolB) | (((~(lStrandCap)) & tStrand)));//todo AUTO
                UnsignedMachineWord invCond = ~cond;

                tStrands[tEdge + j] = (invCond & tStrand) | (cond & lStrandCap);
                tStrand = tStrand << shift;

                cond = cond << shift;
                invCond = ~cond;
                lStrands[lEdge + j] = (invCond & lStrand) | (cond & tStrand);
            }
        }


        inline void loopLowerHalfBinary(int lowerBound, int upperBound, int shift, int lEdge, int tEdge,
                                        UnsignedMachineWord activeBits, UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands,
                                        UnsignedMachineWord *aReverse,
                                        UnsignedMachineWord *b) {

            for (int j = lowerBound; j < upperBound; ++j) {

                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                UnsignedMachineWord lStrandCap = lStrand << (shift + 1);
                auto symbolA = aReverse[lEdge + j];
                auto symbolB = b[tEdge + j];
                UnsignedMachineWord cond = activeBits & (((~(symbolA << (shift + 1))) ^ symbolB) | (((~(lStrandCap)) & tStrand)));
                UnsignedMachineWord invCond = ~cond;

                tStrands[tEdge + j] = (invCond & tStrand) | (cond & lStrandCap);

                tStrand = tStrand >> (shift + 1);

                cond = cond >> (shift + 1);
                invCond = ~cond;
                lStrands[lEdge + j] = (invCond & lStrand) | (cond & tStrand);
            }
        }

        /**
        * Process upper triangles of squares of size sizeof(Input)*8 that lies on specific antidiagonal
        */
        inline void
        loopUpperHalfBinaryMpi(int lowerBound, int upperBound, int shift, int lEdge, int tEdge,
                               UnsignedMachineWord activeBits, UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands, UnsignedMachineWord *aReverse,
                               UnsignedMachineWord *b) {

#pragma omp  for simd schedule(static) aligned(lStrands, tStrands, aReverse, b:sizeof(UnsignedMachineWord)*8)
            for (int j = lowerBound; j < upperBound; ++j) {

                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                auto lStrandCap = lStrand >> shift;
                UnsignedMachineWord symbolA = aReverse[lEdge + j];
                UnsignedMachineWord symbolB = b[tEdge + j];

                UnsignedMachineWord cond = activeBits & (((~(symbolA >> shift)) ^ symbolB) | (((~(lStrandCap)) & tStrand)));
                UnsignedMachineWord invCond = ~cond;

                tStrands[tEdge + j] = (invCond & tStrand) | (cond & lStrandCap);
                tStrand = tStrand << shift;

                cond = cond << shift;
                invCond = ~cond;
                lStrands[lEdge + j] = (invCond & lStrand) | (cond & tStrand);
            }
        }

        /**
         * Same as loopUpperHalfBinary but without offset requirement
         */
        inline void loopCenterHalfBinaryMpi(int lowerBound, int upperBound, int lEdge, int tEdge,
                                            UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands, UnsignedMachineWord *aReverse,
                                            UnsignedMachineWord *b) {
#pragma omp  for simd schedule(static) aligned(lStrands, tStrands, aReverse, b:sizeof(UnsignedMachineWord)*8)
            for (int j = lowerBound; j < upperBound; ++j) {


                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                UnsignedMachineWord cond = ((~(aReverse[lEdge + j] ^ b[tEdge + j])) | ((~lStrand) & tStrand));
                UnsignedMachineWord revCombingCond = ~cond;

                lStrands[lEdge + j] = (revCombingCond & lStrand) | (cond & tStrand);
                tStrands[tEdge + j] = (revCombingCond & tStrand) | (cond & lStrand);
            }
        }

        inline void loopLowerHalfBinaryMpi(int lowerBound, int upperBound, int shift, int lEdge, int tEdge,
                                           UnsignedMachineWord activeBits, UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands,
                                           UnsignedMachineWord *aReverse,
                                           UnsignedMachineWord *b) {

#pragma omp  for simd schedule(static) aligned(lStrands, tStrands, aReverse, b:sizeof(UnsignedMachineWord)*8)
            for (int j = lowerBound; j < upperBound; ++j) {

                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                UnsignedMachineWord lStrandCap = lStrand << (shift + 1);
                auto symbolA = aReverse[lEdge + j];
                auto symbolB = b[tEdge + j];
                UnsignedMachineWord cond = activeBits & (((~(symbolA << (shift + 1))) ^ symbolB) | (((~(lStrandCap)) & tStrand)));
                UnsignedMachineWord invCond = ~cond;

                tStrands[tEdge + j] = (invCond & tStrand) | (cond & lStrandCap);

                tStrand = tStrand >> (shift + 1);

                cond = cond >> (shift + 1);
                invCond = ~cond;
                lStrands[lEdge + j] = (invCond & lStrand) | (cond & tStrand);
            }
        }


        template<int BITS_PER_STRAND, bool USE_FIRST_FORMULA>
        inline void processAntidiagonale(int lowerBound, int upperBound, int lEdge, int tEdge,
                                         UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands, UnsignedMachineWord const *aReverse,
                                         UnsignedMachineWord const *b) {

            auto constexpr IS_BINARY = BITS_PER_STRAND == 1;
            constexpr  int ALIGN_SIZE = sizeof(UnsignedMachineWord) * 8 - (sizeof(UnsignedMachineWord) * 8 / BITS_PER_STRAND) * BITS_PER_STRAND; // 2
            UnsignedMachineWord constexpr singleStrand = UnsignedMachineWord(1); //
            constexpr int size = sizeof(UnsignedMachineWord) * 8 - BITS_PER_STRAND - ALIGN_SIZE;// 8-2-3 = 5

            constexpr int activeNumBits = (sizeof(UnsignedMachineWord) * 8) - ALIGN_SIZE - BITS_PER_STRAND;//CHANGE om bits
            constexpr UnsignedMachineWord braidOnesMask = getBraidMaskOnesMSB<UnsignedMachineWord>(ALIGN_SIZE, BITS_PER_STRAND);


//#pragma omp   for  simd schedule(static)  aligned(lStrands, tStrands:sizeof(UnsignedMachineWord)*8) aligned(aReverse, b:sizeof(UnsignedMachineWord)*8)
            for (int j = lowerBound; j < upperBound; ++j) {
                std::cout<<"iteration "<<j<<std::endl;
                UnsignedMachineWord tStrandCap;
                UnsignedMachineWord lStrandCap, cond;

                //load
                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                auto symbolA = aReverse[lEdge + j];
                auto symbolB = b[tEdge + j];

                UnsignedMachineWord mask = singleStrand;

                UNROLL_LOOP(sizeof(UnsignedMachineWord) * 8 + 1)
                for (int shift = activeNumBits; shift > 0; shift -= BITS_PER_STRAND) {
                    // 15 operations inside cycle
                    // could be reduced to 14 if we store not a but ~a in memory

                    // shift = 5
                    // shift

                    lStrandCap = lStrand >> shift;
                    tStrandCap = tStrand << shift;

                    if constexpr(!IS_BINARY) {
                        std::cout<<"YU"<<std::endl;
                        //reduction block
                        cond = ~(((symbolA >> shift)) ^ symbolB);
                        UnsignedMachineWord eq = cond;
                        PRIVATE_UNROLL_LOOP(COMPARISON_LOOP_UNROLL)
                        for (int i = 1; i < BITS_PER_STRAND; i++) cond &= (eq >> i);
                        std::cout<<"oldL"<<std::bitset<16>(lStrand)<<std::endl;
                        std::cout<<"oldR"<<std::bitset<16>(tStrand)<<std::endl;

                        tStrand = (lStrandCap | ((braidOnesMask ^ (mask)))) & (tStrand | (cond & (mask)));
                        lStrand = tStrandCap ^ (tStrand << shift) ^ lStrand;
//                        std::cout<<"mask"<<std::bitset<16>(braidOnesMask ^ mask)<<","<<std::bitset<16>(cond & (mask))<<",!"<<std::bitset<16>(lStrandCap)<<std::endl;
                        std::cout<<"L"<<std::bitset<16>(lStrand)<<std::endl;
                        std::cout<<"R"<<std::bitset<16>(tStrand)<<std::endl;
                        std::cout<<"invmadsk"<<std::bitset<16>((braidOnesMask ^ (mask))>>ALIGN_SIZE)<<std::endl;
                    } else {
                        if constexpr(USE_FIRST_FORMULA) {
                            /**
                             * The idea is as follows.
                             * to process some antidiagonal that have been built upon Input cells (that contains a batch of strands) we have to
                             * process each such square (Input \times Input) in the antidiagonal fashion.
                             * While processing each antidiagonal of such square, we need to implement the following logic:
                             * in each iteration swap only those bit-strands that  active in iteration &  satisfy the combing condition.
                             * This logic can be implemented by several formulas, here is presented one of it.
                             * There is 20 operations inside cycle
                             */
                            cond = mask & ((~(((symbolA >> shift)) ^ symbolB)) | (((~(lStrandCap)) & tStrand)));
                            UnsignedMachineWord invCond = ~cond;

                            UnsignedMachineWord tStrandShifted = tStrand << shift;
                            tStrand = (invCond & tStrand) | (cond & lStrandCap);

                            cond <<= shift;
                            invCond = ~cond;
                            lStrand = (invCond & lStrand) | (cond & tStrandShifted);
                        } else {
                            /**
                             * The idea is as follows.
                             * to process some antidiagonal that have been built upon Input cells (that contains a batch of strands) we have to
                             * process each such square (Input \times Input) in the antidiagonal fashion.
                             * While processing each antidiagonal of such square, we need to implement the following logic:
                             * in each iteration swap only those bit-strands that  active in iteration &  satisfy the combing condition.
                             * This logic can be implemented by several formulas, here is presented one of it.
                             * There is 15 operations inside cycle (for fast formula)
                             */
                            cond = ~((symbolA >> shift) ^ symbolB);
                            tStrand = (lStrandCap | (~mask)) & (tStrand | (cond & mask));
                            lStrand = tStrandCap ^ (tStrand << shift) ^ lStrand;
                        }
                    }
                    mask = (mask << BITS_PER_STRAND) | singleStrand;
                }

                if constexpr(!IS_BINARY) {
                    cond = ~((symbolA) ^ symbolB);
                    UnsignedMachineWord eq = cond;

                    std::cout<<"SCConda"<<std::bitset<16>(cond)<<std::endl;
                    std::cout<<"A"<<std::bitset<16>(symbolA)<<std::endl;
                    std::cout<<"B"<<std::bitset<16>(symbolB)<<std::endl;


                    PRIVATE_UNROLL_LOOP(COMPARISON_LOOP_UNROLL)
                    for (int i = 1; i < BITS_PER_STRAND; i++) cond &= (eq >> i);
                    std::cout<<"CConda"<<std::bitset<16>(cond)<<std::endl;

                    std::cout<<"oldL"<<std::bitset<16>(lStrand)<<std::endl;
                    std::cout<<"oldR"<<std::bitset<16>(tStrand)<<std::endl;

                    lStrandCap = lStrand;
//                    tStrandCap = tStrand;
                    tStrandCap = tStrand;


                    tStrand = (lStrandCap) & (tStrand | (cond & (braidOnesMask)));//TODO simplify mask>>ALIGN

                    lStrand = tStrandCap ^ (tStrand) ^ lStrand;
//                    std::cout<<"mask"<<std::bitset<16>(mask)<<","<<std::bitset<16>(cond & (mask>>ALIGN_SIZE))<<",!"<<std::bitset<16>(lStrandCap)<<std::endl;
                    std::cout<<"L"<<std::bitset<16>(lStrand)<<std::endl;
                    std::cout<<"RR"<<std::bitset<16>(tStrand)<<std::endl;
//                    std::cout<<"RRMMMASK"<<std::bitset<16>(mask)<<std::endl;


                    mask = braidOnesMask;

                } else {
                    if constexpr(USE_FIRST_FORMULA) {
                        // center
                        cond = (~(symbolA ^ symbolB));
                        cond = (cond | ((~lStrand) & tStrand));
                        UnsignedMachineWord invCond = ~cond;
                        UnsignedMachineWord tStrandShifted = tStrand;
                        tStrand = (invCond & tStrand) | (cond & lStrand);
                        lStrand = (invCond & lStrand) | (cond & tStrandShifted);
                    } else {
                        // center, no shifts
                        cond = ~((symbolA ^ symbolB));
                        lStrandCap = lStrand;
                        tStrandCap = tStrand;

                        tStrand = (lStrandCap | (~mask)) & (tStrand | (cond & mask));
                        lStrand = tStrandCap ^ (tStrand) ^ lStrand;
                    }
                    mask = ~UnsignedMachineWord(0);
                }

                UNROLL_LOOP(sizeof(UnsignedMachineWord) * 8 + 1)
//                for (int shift = BITS_PER_STRAND - ALIGN_SIZE; shift < activeNumBits + 1 - ALIGN_SIZE; shift += BITS_PER_STRAND) {
                    for (int shift = BITS_PER_STRAND; shift < activeNumBits + 1; shift += BITS_PER_STRAND) {

                        mask <<= BITS_PER_STRAND;//SOMEWHWRE ALIGN TODO

                    lStrandCap = lStrand << shift;
                    tStrandCap = tStrand >> shift;

                    if constexpr(!IS_BINARY) {
                        //reduction block
                        cond = ~((symbolA << shift) ^ (symbolB));
                        UnsignedMachineWord eq = cond;
                        std::cout<<"Cond"<<std::bitset<16>(cond)<<std::endl;

                        PRIVATE_UNROLL_LOOP(COMPARISON_LOOP_UNROLL)
                        for (int i = 1; i < BITS_PER_STRAND; i++) cond &= (eq >> i);
                        std::cout<<"Cond"<<std::bitset<16>(cond)<<std::endl;
                        lStrandCap = lStrand << (shift);
                        tStrandCap = tStrand >> (shift);

                        std::cout<<"oldL"<<std::bitset<16>(lStrand)<<std::endl;
                        std::cout<<"oldR"<<std::bitset<16>(tStrand)<<std::endl;

                        tStrand = (lStrandCap | ((braidOnesMask ^ mask))) & (tStrand | ( ((cond & mask))));
                        lStrand = (tStrandCap ^ (tStrand >> (shift)) ^ lStrand);

                        std::cout<<"M"<<std::bitset<16>(mask)<<std::endl;
                        std::cout<<"T"<<std::bitset<16>(((braidOnesMask) ^ mask))<<","<<std::bitset<16>(braidOnesMask)<<","<<std::bitset<16>(((cond & mask)))<<std::endl;
                        std::cout<<"L"<<std::bitset<16>(lStrand)<<std::endl;
                        std::cout<<"R"<<std::bitset<16>(tStrand)<<std::endl;
                    } else {
                        if constexpr(USE_FIRST_FORMULA) {
                            lStrandCap = lStrand << (shift);
                            cond = ~(((symbolA << ((shift))) ^ symbolB)); // NOT A XOR B = NOT (A XOR B)// ECONOMY

                            cond = mask & (cond | (((~(lStrandCap)) & tStrand)));
                            UnsignedMachineWord invCond = ~cond;

                            UnsignedMachineWord tStrandShifted = tStrand >> shift;
                            tStrand = (invCond & tStrand) | (cond & lStrandCap);
                            cond >>= shift;
                            invCond = ~cond;

                            lStrand = (invCond & lStrand) | (cond & tStrandShifted);
                        } else {
                            lStrandCap = lStrand << shift;
                            tStrandCap = tStrand >> shift;
                            cond = ~(((symbolA << (shift)) ^ symbolB));
                            tStrand = (lStrandCap | (~mask)) & (tStrand | (cond & mask));
                            lStrand = tStrandCap ^ (tStrand >> shift) ^ lStrand;
                        }
                    }
                }

                // store
                lStrands[lEdge + j] = lStrand;
                tStrands[tEdge + j] = tStrand;
            }
        }


        template<class Input>
        inline void processAntidiagonal(int lowerBound, int upperBound, int lEdge, int tEdge,
                                        Input *lStrands, Input *tStrands, Input const *aReverse, Input const *b,
                                        int residue, int bitsPerStrand, Input braidOnes) {
            Input singleStrand = Input(1) << residue;
            int size = sizeof(Input) * 8 - residue - bitsPerStrand;

#pragma omp   for  simd schedule(static)  aligned(lStrands, tStrands:sizeof(Input)*8) aligned(aReverse, b:sizeof(Input)*8)
            for (int j = lowerBound; j < upperBound; ++j) {

                Input lStrandCap, cond, tStrandCap, eq;
                //load phase
                Input lStrand = lStrands[lEdge + j];
                Input tStrand = tStrands[tEdge + j];
                Input symbolA = aReverse[lEdge + j];
                Input symbolB = b[tEdge + j];

                Input mask = singleStrand;

                // manual say 256 just for complete
#pragma GCC unroll  256
                for (int shift = size; shift > 0; shift -= bitsPerStrand) {


                    lStrandCap = lStrand >> shift;
                    tStrandCap = tStrand << shift;

                    //reduction block
                    cond = ~(((symbolA >> shift)) ^ symbolB);
                    eq = cond;
#pragma GCC unroll  10
                    for (int i = 1; i < bitsPerStrand; i++) {
                        cond &= (eq >> i);
                    }

                    tStrand = (lStrandCap | (braidOnes ^ mask)) & (tStrand | (cond & mask));
                    lStrand = tStrandCap ^ (tStrand << shift) ^ lStrand;

                    mask = (mask << bitsPerStrand) | singleStrand;

                }

                cond = ~((symbolA) ^ symbolB);
                eq = cond;

#pragma GCC unroll  10
                for (int i = 1; i < bitsPerStrand; i++) cond &= (eq >> i);

                lStrandCap = lStrand;
                tStrandCap = tStrand;


                tStrand = (lStrandCap | (braidOnes ^ braidOnes)) & (tStrand | (cond & braidOnes));


                lStrand = tStrandCap ^ tStrand ^ lStrand;

                mask = braidOnes << bitsPerStrand;

                //lower half
#pragma GCC unroll 256
                for (int shift = bitsPerStrand; shift <= size; shift += bitsPerStrand) {

                    //reduction block
                    cond = ~((symbolA << shift) ^ symbolB);
                    eq = cond;
#pragma GCC unroll  10
                    for (int i = 1; i < bitsPerStrand; i++) cond &= (eq >> i);

                    lStrandCap = lStrand << shift;
                    tStrandCap = tStrand >> shift;

                    tStrand = (lStrandCap | (braidOnes ^ mask)) & (tStrand | (cond & mask));
                    lStrand = (tStrandCap ^ (tStrand >> shift) ^ lStrand);

                    mask <<= bitsPerStrand;
                }

                lStrands[lEdge + j] = lStrand;
                tStrands[tEdge + j] = tStrand;
            }
        }


    };


    template<class UnsignedMachineWord, class InputType, bool IS_FILL_MACHINE_WORD = true, bool IS_ONLY_BINARY_ALPHABET = true, bool USE_NAIVE_LAYOUT = false, bool USE_FAST_FORMULA = true>
    class BitSemiLocalSticky : public AbstractBitSemiLocalSticky<UnsignedMachineWord, InputType> {

    private:
        int threadsNum;

    public:
        BitSemiLocalSticky(int threads) : threadsNum(threads) {}

        int compute(const InputType *a, int aSize, const InputType *b, int bSize) override {
            if (aSize > bSize) return compute(b, bSize, a, aSize);

            auto converter = AlphabetConverter();
            std::unordered_set<InputType> alphabet;
            for (int i = 0; i < aSize; ++i) alphabet.insert(a[i]);
            for (int i = 0; i < bSize; ++i) alphabet.insert(b[i]);
            alphabet.insert({7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27});
            auto[forwardMapper, inverseMapper] = converter.encodeAlphabet(alphabet);
            std::cout<<"SET"<<std::endl;
            for (auto[k, v]: forwardMapper) {
                std::cout<<k<<"->"<<std::bitset<3>(v)<<std::endl;
            }

            auto[packedRevA, packedADetails] = converter.packSequence<InputType, UnsignedMachineWord, false, true>(a, aSize, forwardMapper, alphabet.size());
            auto[packedB, packedBDetails] = converter.packSequence<InputType, UnsignedMachineWord, true, false>(b, bSize, forwardMapper, alphabet.size());
            int res = 0;

            if constexpr(IS_ONLY_BINARY_ALPHABET) {
                if constexpr(USE_NAIVE_LAYOUT) {
                    res = llcs2SymbolNaiveCombing(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                } else {
                    switch (packedADetails.bitsPerSymbol) {
                        case 1:
                            res = llcs2SymbolSmartCombing<USE_FAST_FORMULA, 1>(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                            break;
                        case 2:
                            std::cout<<"4 symbol alphabet"<<std::endl;
//                            res = llcs_nary_symbol_smart_combing(packedRevA,packedADetails.getCompressedSizeInWords(),packedB,packedBDetails.getCompressedSizeInWords(),packedADetails.initalNumSymbols,2);
                            res = llcs2SymbolSmartCombing<USE_FAST_FORMULA, 2>(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                            break;
                        case 3:
                            std::cout<<"8 symbol alphabet"<<std::endl;
                            res = llcs2SymbolSmartCombing<USE_FAST_FORMULA, 3>(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                            break;
                        case 4:
                            std::cout<<"16 symbol alphabet"<<std::endl;

                            res = llcs2SymbolSmartCombing<USE_FAST_FORMULA, 4>(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                            break;
                        case 5:
                            std::cout<<"32 symbol alphabet"<<std::endl;

                            res = llcs2SymbolSmartCombing<USE_FAST_FORMULA, 5>(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                            break;
                        case 6:
                            res = llcs2SymbolSmartCombing<USE_FAST_FORMULA, 6>(packedRevA, packedADetails, packedB, packedBDetails, threadsNum);
                            break;
                        default:
                            exit(-2);
                    }
                }
            } else {
                throw std::runtime_error("Not implemented yet");
            }

            delete[] packedRevA;
            delete[] packedB;
            return res;
        }

    private:
        /**
         * This is the first non-optimized version of algorithm
         * Several condition should be satisfied:
         * 1) Since both strings are binary the elements should be compressed and stored in bits
         * 2) The first string assumed be less then the second one
         * 3) First string stored in reverse order to allow us to access elements in cell processing routine in sequential manner
         * 4) Input, that store pack of bits, should be unsigned to eliminate problem with signed shift
         * 5) Size of strings should be divisible by sizeof(Input)*8, if no, see details in paper or implementation
         * Algorithm follows the idea of iterative combing  algorithm but strands have only 0 and 1 number.
         * To see  cell processing routine see documentation of methods that used within the algorithm.
         * @tparam Input
         * @param aReverse
         * @param aSize
         * @param b
         * @param bSize
         * @param aTotalSymbols
         * @param threads_num
         * @return
         */
        int llcs2SymbolNaiveCombing(UnsignedMachineWord *aReverse, AlphabetConverter::ConversionDetails<UnsignedMachineWord> aDetails,
                                    UnsignedMachineWord *b, AlphabetConverter::ConversionDetails<UnsignedMachineWord> bDetails, int thdsNum) {
            // also stored in the reverse order
            auto *lStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), aDetails.getCompressedSizeInBytes()));
            auto *tStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), bDetails.getCompressedSizeInBytes()));

            auto m = aDetails.getCompressedSizeInWords(), n = bDetails.getCompressedSizeInWords();
            // total amount of strands that at the end hit right border of grid
            int disBraid = 0;


            auto numDiag = m + n - 1;

            auto totalSameLengthDiag = numDiag - (m) - (m - 1);

            auto[bottom, right] =  getPreciseBottomBorder(aDetails, bDetails);

            auto bBorderStart = *bDetails.getNthPositionOfSymbol(bottom.first);
            auto bBorderEnd = *bDetails.getNthPositionOfSymbol(bottom.second - 1);
            auto rBorder = *aDetails.getNthPositionOfSymbol(right.second);

            UnsignedMachineWord braidOnes = UnsignedMachineWord(1) << aDetails.getWordAlignmentInBits();
            for (int i = 0; i < aDetails.machineWordSize; i += aDetails.bitsPerSymbol) braidOnes |= (braidOnes << i);


            auto bitsPerChar = aDetails.bitsPerSymbol;
            lStrands[0] = (~aDetails.getPaddingMaskOfLastWord(false)) & (~aDetails.getAlignmentMask());
            tStrands[n - 1] = (bDetails.getPaddingMaskOfLastWord(true)) & (~bDetails.getAlignmentMask());

#pragma omp parallel num_threads(thdsNum)  default(none) shared(lStrands, tStrands, aReverse, b, m, n, disBraid, totalSameLengthDiag, aDetails, bDetails, rBorder, bBorderStart, bBorderEnd, braidOnes)
            {
                initStep(lStrands, m, tStrands, n, braidOnes);

                UnsignedMachineWord mask;
                auto upperBound = (sizeof(UnsignedMachineWord) * 8) - 1;

                //process first triangle in 0,0 cube
                mask = UnsignedMachineWord(1);
                UnsignedMachineWord maskR = UnsignedMachineWord(1) << (sizeof(UnsignedMachineWord) * 8 - 1);
                int bitsShift = (sizeof(UnsignedMachineWord) * 8 - 1);

                //PHASE 0:Process first triangle
#pragma omp single
                {
                    for (int insideDiagNum = 0; insideDiagNum <= upperBound; ++insideDiagNum, bitsShift--) {
                        this->loopUpperHalfBinary(0, 1, bitsShift, m - 1, 0, mask, lStrands, tStrands, aReverse, b);
                        mask = (mask << 1) | UnsignedMachineWord(1);
                    }
                }

                //PHASE 1: Process diagonal till fill big left triangle,
                for (int curDiagCubeLen = 1; curDiagCubeLen < m; curDiagCubeLen++) {
                    //to process current
                    bitsShift = (sizeof(UnsignedMachineWord) * 8 - 1);
                    mask = UnsignedMachineWord(1);
                    //to process previous
                    UnsignedMachineWord maskPrev = ~static_cast<UnsignedMachineWord>(0);

                    //process previous size/2 - 1 cubes and current size/2 -1  cubes
                    for (int insideDiagNum = 0; insideDiagNum < upperBound; ++insideDiagNum, bitsShift--) {
                        //update mask of prev move
                        maskPrev <<= 1;

                        this->loopUpperHalfBinaryMpi(0, curDiagCubeLen + 1, bitsShift, m - 1 - curDiagCubeLen, 0, mask, lStrands, tStrands, aReverse,
                                                     b);

                        this->loopLowerHalfBinaryMpi(0, curDiagCubeLen, insideDiagNum, m - 1 - curDiagCubeLen + 1, 0, maskPrev, lStrands, tStrands,
                                                     aReverse, b);


                        //update mask of current move
                        mask = (mask << 1) | UnsignedMachineWord(1);
                    }

                    this->loopCenterHalfBinaryMpi(0, curDiagCubeLen + 1, m - 1 - curDiagCubeLen, 0, lStrands, tStrands, aReverse, b);

                }

                //PHASE 2
                for (int k = 0; k < totalSameLengthDiag; ++k) {

                    //to process current
                    bitsShift = (sizeof(UnsignedMachineWord) * 8 - 1);
                    mask = UnsignedMachineWord(1);

                    //to process previous
                    UnsignedMachineWord maskPrev = ~UnsignedMachineWord(0);

                    for (int insideDiagNum = 0; insideDiagNum < upperBound; ++insideDiagNum, bitsShift--) {
                        //update mask of prev move
                        maskPrev <<= 1;

                        this->loopUpperHalfBinaryMpi(0, m, bitsShift, 0, k + 1, mask, lStrands, tStrands, aReverse, b);

                        this->loopLowerHalfBinaryMpi(0, m, insideDiagNum, 0, k, maskPrev, lStrands, tStrands, aReverse, b);


                        //update mask of current move
                        mask = (mask << 1) | UnsignedMachineWord(1);
                    }

                    this->loopCenterHalfBinaryMpi(0, m, 0, k + 1, lStrands, tStrands, aReverse, b);
                }


                auto startJ = totalSameLengthDiag + 1;
                for (int curDiagCubeLen = m - 1; curDiagCubeLen >= 1; curDiagCubeLen--, startJ++) {

                    //to process current
                    bitsShift = (sizeof(UnsignedMachineWord) * 8 - 1);
                    mask = UnsignedMachineWord(1);

                    //to process previous
                    UnsignedMachineWord maskPrev = ~UnsignedMachineWord(0);

                    //process previous size/2 - 1 cubes and current size/2 -1  cubes
                    for (int insideDiagNum = 0; insideDiagNum < upperBound; ++insideDiagNum, bitsShift--) {
                        //update mask of prev move
                        maskPrev <<= 1;

                        this->loopUpperHalfBinaryMpi(0, curDiagCubeLen, bitsShift, 0, startJ, mask, lStrands, tStrands, aReverse, b);

                        this->loopLowerHalfBinaryMpi(0, curDiagCubeLen + 1, insideDiagNum, 0, startJ - 1, maskPrev, lStrands, tStrands, aReverse, b);
                        //update mask of current move
                        mask = (mask << 1) | UnsignedMachineWord(1);

                    }

                    this->loopCenterHalfBinaryMpi(0, curDiagCubeLen, 0, startJ, lStrands, tStrands, aReverse, b);
                }


                //process last triangle in position  m-1, n-1  cube
                mask = ~UnsignedMachineWord(0);
                maskR = mask;
#pragma omp single
                {

                    for (int insideDiagNum = 0; insideDiagNum < upperBound; ++insideDiagNum) {
                        mask = mask << 1;
                        this->loopLowerHalfBinary(0, 1, insideDiagNum, 0, n - 1, mask, lStrands, tStrands, aReverse, b);
                    }

                }

                int ones = (m < n) ? countOnes(lStrands, rBorder.first + 1, m) : countOnes(tStrands, bBorderStart.first + 1, bBorderEnd.first);

#pragma omp atomic update
                disBraid += ones;

            }

            if ((m < n)) {
                disBraid += countOnesInWord(lStrands[0], bitsPerChar * (rBorder.second) + bDetails.getWordAlignmentInBits());
                disBraid += countOnesInWord(tStrands[n - 1], bBorderEnd.second * bitsPerChar + 1);
                disBraid = aDetails.initalNumSymbols - disBraid + bDetails.getPaddingOfLastWordInSymbols();
            } else {
                if (bBorderStart.first == bBorderEnd.first) {
                    //TODO
                    disBraid += countOnesInWord(tStrands[bBorderStart.first], bBorderStart.second * bitsPerChar, (bBorderEnd.second) * bitsPerChar + 1);
                } else {
                    disBraid += countOnesInWord(tStrands[bBorderStart.first], bBorderStart.second * bitsPerChar);
                    disBraid += countOnesInWord(tStrands[bBorderEnd.first], 0, bBorderEnd.second * bitsPerChar + 1);
                }
                disBraid += countOnesInWord(lStrands[0], 0, bitsPerChar * (rBorder.second) + aDetails.getWordAlignmentInBits()); // offset
            }

            delete[] lStrands;
            delete[] tStrands;

            return disBraid;
        }

        std::pair<std::pair<int, int>, std::pair<int, int>> getPreciseBottomBorder(AlphabetConverter::ConversionDetails<UnsignedMachineWord> &aDetails,
                                                                                   AlphabetConverter::ConversionDetails<UnsignedMachineWord> &bDetails) {
            auto extCharSizeA = aDetails.getNumberOfEmbeddingsInWord() * aDetails.getCompressedSizeInWords();//wrt to added paddings
            auto extCharSizeB = bDetails.getNumberOfEmbeddingsInWord() * bDetails.getCompressedSizeInWords();

            auto lStart = aDetails.getPaddingOfLastWordInSymbols();
            auto lEnd = std::min(lStart + bDetails.initalNumSymbols, extCharSizeB);
            auto rEnd = 0;
            if (lEnd == extCharSizeB) {
                rEnd = std::max(lStart + bDetails.initalNumSymbols, extCharSizeB) - lEnd;
            }
            return std::make_pair(std::make_pair(lStart, lEnd), std::make_pair(0, rEnd));
        }


        template<bool USE_FAST_COND, int BITS_PER_SYMBOL>
        int llcs2SymbolSmartCombing(UnsignedMachineWord *aReverse, AlphabetConverter::ConversionDetails<UnsignedMachineWord> aDetails,
                                    UnsignedMachineWord *b, AlphabetConverter::ConversionDetails<UnsignedMachineWord> bDetails, int thdsNum) {

            // also stored in the reverse order
            auto *lStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), aDetails.getCompressedSizeInBytes()));
            auto *tStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), bDetails.getCompressedSizeInBytes()));

            auto m = aDetails.getCompressedSizeInWords(), n = bDetails.getCompressedSizeInWords();
            // total amount of strands that at the end hit right border of grid
            int disBraid = 0;


            auto numDiag = m + n - 1;

            auto totalSameLengthDiag = numDiag - (m - 1) - (m - 1);

            auto[bottom, right] =  getPreciseBottomBorder(aDetails, bDetails);

            auto bBorderStart = *bDetails.getNthPositionOfSymbol(bottom.first);
            auto bBorderEnd = *bDetails.getNthPositionOfSymbol(bottom.second - 1);
            auto rBorder = *aDetails.getNthPositionOfSymbol(right.second);

            auto constexpr ALIGN_SIZE = sizeof(UnsignedMachineWord) * 8 - (sizeof(UnsignedMachineWord) * 8 / BITS_PER_SYMBOL) * BITS_PER_SYMBOL;
            std::cout<<ALIGN_SIZE<<","<<BITS_PER_SYMBOL<<std::endl;
            auto constexpr braidOnes = getBraidMaskOnesMSB<UnsignedMachineWord>(ALIGN_SIZE, BITS_PER_SYMBOL);
            std::cout<<"HUI"<<std::bitset<16>(braidOnes)<<std::endl;

            auto bitsPerChar = aDetails.bitsPerSymbol;
            lStrands[0] = ((~aDetails.getPaddingMaskOfLastWord(false)) & (~aDetails.getAlignmentMask())) & braidOnes;
            tStrands[n - 1] = (bDetails.getPaddingMaskOfLastWord(true)) & (~bDetails.getAlignmentMask())  & (braidOnes);

            std::cout<<"maskA"<<std::bitset<16>(aDetails.getPaddingMaskOfLastWord(false))<<"AD"<< std::bitset<16>(aDetails.getAlignmentMask())<<std::endl;
            std::cout<<"braid:"<<std::bitset<16>(braidOnes)<<"lStrand:"<<std::bitset<16>(lStrands[0])<<"tStrand:"<<std::bitset<16>(tStrands[n-1])<<","<<std::endl;
            std::cout<<"LAST:"<<std::bitset<16>(bDetails.getPaddingMaskOfLastWord(true))<<"C"<<std::bitset<16>(bDetails.getAlignmentMask())<<","<<std::endl;

#pragma omp parallel num_threads(1)  default(none) shared(lStrands, tStrands, aReverse, b, m, n, disBraid, totalSameLengthDiag, aDetails, bDetails, rBorder, bBorderStart, bBorderEnd, braidOnes,std::cout)
            {
                initStep(lStrands, m, tStrands, n, braidOnes);
                // phase 1: process upper left triangle

                std::cout<<"A:"<<std::endl;
                for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(aReverse[i])<<",";
                std::cout<<"B:"<<std::endl;
                for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(b[i])<<",";
                std::cout<<"Lstrands:"<<std::endl;
                for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(lStrands[i])<<",";
                std::cout<<"tStradns:"<<std::endl;
                for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(tStrands[i])<<",";
                std::cout<<std::endl;


                for (int diagLen = 0; diagLen < m - 1; diagLen++) {
                    this->template processAntidiagonale<BITS_PER_SYMBOL, !USE_FAST_COND>(0, diagLen + 1, m - 1 - diagLen, 0, lStrands, tStrands, aReverse, b);
                    std::cout<<"1stphase:"<<std::endl;
                    std::cout<<"Lstrands:"<<std::endl;
                    for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(lStrands[i])<<",";
                    std::cout<<"tStradns:"<<std::endl;
                    for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(tStrands[i])<<",";
                    std::cout<<std::endl;
                }

                // phase2: process parallelogram
                for (int k = 0; k < totalSameLengthDiag; k++) {
                    this->template processAntidiagonale<BITS_PER_SYMBOL, !USE_FAST_COND>(0, m, 0, k, lStrands, tStrands, aReverse, b);
                    std::cout<<"2stphase:"<<totalSameLengthDiag<<std::endl;
                    std::cout<<"Lstrands:"<<std::endl;
                    for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(lStrands[i])<<",";
                    std::cout<<"tStradns:"<<std::endl;
                    for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(tStrands[i])<<",";
                    std::cout<<std::endl;
                }

                auto startJ = totalSameLengthDiag;

                // phase:3: lower-right triangle
                for (int diagLen = m - 1; diagLen >= 1; diagLen--) {
                    this->template processAntidiagonale<BITS_PER_SYMBOL, !USE_FAST_COND>(0, diagLen, 0, startJ, lStrands, tStrands, aReverse, b);
                    startJ++;

                    std::cout<<"3stphase:"<<totalSameLengthDiag<<std::endl;
                    std::cout<<"Lstrands:"<<std::endl;
                    for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(lStrands[i])<<",";
                    std::cout<<"tStradns:"<<std::endl;
                    for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(tStrands[i])<<",";
                    std::cout<<std::endl;
                }

                std::cout<<"Borders:"<<rBorder.first + 1<<","<< m<< "or "<<bBorderStart.first + 1<<","<< bBorderEnd.first<<","<<std::endl;

                std::cout<<" countOnes(lStrands," <<rBorder.first + 1<<", "<< m<<")"<<std::endl;
                int ones = (m <= n) ? countOnes(lStrands, rBorder.first + 1, m) : countOnes(tStrands, bBorderStart.first + 1, bBorderEnd.first);

#pragma omp atomic update
                disBraid += ones;
            }

            if ((m <= n)) {
                std::cout<<"!Borders:"<<rBorder.first<<":"<<  rBorder.second<<","<<(bBorderEnd.second)<<std::endl;
                std::cout<<"Count:"<<bitsPerChar * (rBorder.second)<<bDetails.getWordAlignmentInBits()<<std::endl;
                disBraid += countOnesInWord(lStrands[0], bitsPerChar * (rBorder.second));
                std::cout<<"countOnesInWord(lStrands[0],"<<bitsPerChar * (rBorder.second) + bDetails.getWordAlignmentInBits()<<")"<<std::endl;
                disBraid += countOnesInWord(tStrands[n - 1], (bBorderEnd.second + 1)*bitsPerChar);
                std::cout<<"countOnesInWord(tStrands[n - 1]"<<(bBorderEnd.second + 1)*bitsPerChar<<")"<<std::endl;
                std::cout<<"Count:"<< (bBorderEnd.second +1)* (bitsPerChar)<<std::endl;
                std::cout<<"PADDINGL"<<bDetails.getPaddingOfLastWordInSymbols()<<std::endl;
                disBraid = aDetails.initalNumSymbols - disBraid + bDetails.getPaddingOfLastWordInSymbols();
            } else {
                if (bBorderStart.first == bBorderEnd.first) {
                    //TODO
                    std::cout<<"==Borders:"<<bBorderStart.first<<":"<<  bBorderStart.second<<","<<(bBorderEnd.second)<<std::endl;
                    disBraid += countOnesInWord(tStrands[bBorderStart.first], bBorderStart.second * bitsPerChar, (bBorderEnd.second+1) * bitsPerChar);
                } else {
                    disBraid += countOnesInWord(tStrands[bBorderStart.first], bBorderStart.second * bitsPerChar);
                    disBraid += countOnesInWord(tStrands[bBorderEnd.first], 0, (bBorderEnd.second+1) * (bitsPerChar));
                }
                disBraid += countOnesInWord(lStrands[0], 0, bitsPerChar * (rBorder.second) ); // offset
            }

            std::cout<<"finale:"<<std::endl;
            std::cout<<"A:"<<std::endl;
            for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(aReverse[i])<<",";
            std::cout<<"B:"<<std::endl;
            for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(b[i])<<",";
            std::cout<<"Lstrands:"<<std::endl;
            for(int i = 0;i<m;i++) std::cout<<std::bitset<16>(lStrands[i])<<",";
            std::cout<<"tStradns:"<<std::endl;
            for(int i = 0;i<n;i++) std::cout<<std::bitset<16>(tStrands[i])<<",";
            std::cout<<std::endl;
            delete[] lStrands;
            delete[] tStrands;

            return disBraid;
        }

        void initStep(UnsignedMachineWord *lStrands, int m, UnsignedMachineWord *tStrands, int n, UnsignedMachineWord bradOne) {
            // Initialization step, strands that hit the left grid border all have number 1; strands that hit top grid border are 0.
#pragma omp  for simd schedule(static) aligned(lStrands:sizeof(UnsignedMachineWord)*8)
            for (int k = 1; k < m; ++k) lStrands[k] = bradOne;
#pragma omp  for simd schedule(static) aligned(tStrands:sizeof(UnsignedMachineWord)*8)
            for (int k = 0; k < n - 1; ++k) tStrands[k] = UnsignedMachineWord(0);
        }

        int countOnes(UnsignedMachineWord *strands, int safe_interval_from, int safe_interval_to) {
            auto ones = 0;
            int localSum = 0;
#pragma omp for  simd schedule(static)  aligned(strands:sizeof(UnsignedMachineWord)*8)
            for (int i1 = safe_interval_from; i1 < safe_interval_to; ++i1) {
                //  Brian Kernighan’s Algorithm
                int counter = 0;
                UnsignedMachineWord number = strands[i1];
                //  LogNumber
                while (number) {
                    number &= (number - 1);
                    counter++;
                }
                localSum += counter;
            }
            return localSum;
        }

        int countOnesInWord(UnsignedMachineWord number, int fromInclusive, int toExclusive = sizeof(UnsignedMachineWord) * 8) {

            if (fromInclusive == toExclusive) {
                return 0;
            }
            if (fromInclusive > toExclusive) {
                std::cout << "fromInclusive > toExclusive" << fromInclusive << ">" << toExclusive << std::endl;
//                exit(-1);
                return 0;
            }
            UnsignedMachineWord maskTakeAllFrom = ~UnsignedMachineWord(0);
            if (fromInclusive != 0) maskTakeAllFrom &= ~((UnsignedMachineWord(1) << (fromInclusive)) - 1);
            UnsignedMachineWord maskTakeAllTo = ~UnsignedMachineWord(0);
            if (toExclusive != sizeof(UnsignedMachineWord) * 8) maskTakeAllTo &= (UnsignedMachineWord(1) << (toExclusive)) - 1;
//            std::cout<<"number"<<std::bitset<16>(number)<<","<<std::endl;
            number = number & maskTakeAllFrom & maskTakeAllTo;
//            std::cout<<std::bitset<16>(maskTakeAllFrom)<<","<<std::bitset<16>(maskTakeAllTo)<<","<<fromInclusive<<","<<toExclusive<<std::endl;
//            std::cout<<std::bitset<16>(number)<<","<<std::endl;
            auto counter = 0;
            while (number) {
                number &= (number - 1);
                counter++;
            }
            return counter;
        }


        int llcs_nary_symbol_smart_combing(UnsignedMachineWord *a_reverse, int a_size, UnsignedMachineWord *b, int b_size,
                                           int a_total_symbols, int bits_per_symbol, int threads_num = 1) {

            std::cout << "EEEE" << ',' << b_size << std::endl;
            auto *l_strands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), sizeof(UnsignedMachineWord) * a_size));
            auto *t_strands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), sizeof(UnsignedMachineWord) * b_size));

            auto m = a_size, n = b_size;

            // total amount of strands that at the end hit right border of grid
            int dis_braid = 0;

            auto num_diag = m + n - 1;

            auto total_same_length_diag = num_diag - (m - 1) - (m - 1);

            int residue = (sizeof(UnsignedMachineWord) * 8) % bits_per_symbol;
            UnsignedMachineWord braid_ones = UnsignedMachineWord(1) << residue;
            for (int i = 0; i < sizeof(UnsignedMachineWord) * 8; i += bits_per_symbol) braid_ones |= (braid_ones << i);

#pragma omp parallel num_threads(threads_num)  default(none) shared(std::cout, residue, bits_per_symbol, braid_ones, l_strands, t_strands, a_reverse, b, m, n, dis_braid, total_same_length_diag)
            {

                // Initialization step, strands that hit the left grid border all have number 1; strands that hit top grid border are 0.
#pragma omp  for simd schedule(static) aligned(l_strands:sizeof(UnsignedMachineWord)*8)
                for (int k = 0; k < m; ++k) l_strands[k] = braid_ones;
#pragma omp  for simd schedule(static) aligned(t_strands:sizeof(UnsignedMachineWord)*8)
                for (int k = 0; k < n; ++k) t_strands[k] = UnsignedMachineWord(0);

                // phase 1: process upper left triangle
                for (int diag_len = 0; diag_len < m - 1; diag_len++) {
                    this->processAntidiagonal(0, diag_len + 1, m - 1 - diag_len, 0, l_strands, t_strands, a_reverse, b,
                                        residue, bits_per_symbol, braid_ones);
                }

                // phase2: process parallelogram
                for (int k = 0; k < total_same_length_diag; k++) {
                    this->processAntidiagonal(0, m, 0, k, l_strands, t_strands, a_reverse, b, residue, bits_per_symbol, braid_ones);
                }

                auto start_j = total_same_length_diag;

                // phase:3: lower-right triangle
                for (int diag_len = m - 1; diag_len >= 1; diag_len--) {
                    this->processAntidiagonal(0, diag_len, 0, start_j, l_strands, t_strands, a_reverse, b,
                                        residue, bits_per_symbol, braid_ones);
                    start_j++;
                }

#pragma omp for  simd schedule(static) reduction(+:dis_braid)  aligned(l_strands:sizeof(UnsignedMachineWord)*8)
                for (int i1 = 0; i1 < m; ++i1) {
//
                    //  Brian Kernighan’s Algorithm
                    int counter = 0;
                    UnsignedMachineWord number = l_strands[i1];
                    //  LogNumber
                    while (number) {
                        number &= (number - 1);
                        counter++;
                    }
                    dis_braid += counter;
                }
            }

            return a_total_symbols - dis_braid;
        }

    };


}