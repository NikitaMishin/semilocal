
#include <bitset>
#include "../LCS/LCS.h"
#include "unordered_set"

namespace bit_parallel::lcs {

    struct AlphabetConverter {

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
            auto bitsPerSymbol = std::max(1,int(std::ceil(log2(alphabetSize))));
            auto shift = bitsPerSymbol;
            auto wordSizeInBits = sizeof(Output) * 8;
            auto symbolsInWord = int(wordSizeInBits / shift);

            auto bytesNeeded = int(std::ceil(size * 1.0 / symbolsInWord) * sizeof(Output));

            auto bitsetArray = static_cast<Output *> (aligned_alloc(sizeof(Output), bytesNeeded));
            auto n = bytesNeeded / sizeof(Output);

            auto residue = wordSizeInBits % bitsPerSymbol;

            //    fill bitset
            for (int i = 0; i < n - 1; ++i) {
                Output word = 0;
                for (int symbol = 0; symbol < symbolsInWord; symbol++) {
                    word |= Output(mapperForward[a[i * symbolsInWord + symbol]]) << (shift * symbol + residue);
                }
                bitsetArray[i] = word;
            }

            //    fill last
            for (int i = n - 1; i < n; ++i) {
                Output word = 0;
                for (int symbol = 0; i * symbolsInWord + symbol < size; symbol++) {
                    word |= Output(mapperForward[a[i * symbolsInWord + symbol]]) << (shift * symbol + residue);
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
        std::pair<Output *, int> packMSBReverse(const Input *a, int size, Map<Input> &mapperForward, int alphabetSize) {

            auto bitsPerSymbol = std::max(1,int(std::ceil(log2(alphabetSize))));
            auto shift = bitsPerSymbol;
            auto wordSizeInBits = sizeof(Output) * 8;
            auto symbolsInWord = int(wordSizeInBits / shift);

            auto bytesNeeded = int(std::ceil(size * 1.0 / symbolsInWord) * sizeof(Output));

            auto bitsetArray = static_cast<Output *> (aligned_alloc(sizeof(Output), bytesNeeded));
            auto n = bytesNeeded / sizeof(Output);
            auto residue = wordSizeInBits % bitsPerSymbol;



            // fill first guy
            for (int i = 0; i < n - 1; ++i) {
                Output word = Output(0);
                for (int symbol = 0; symbol < symbolsInWord; symbol++) {
                    word |= Output(mapperForward[a[(i * symbolsInWord + symbol)]]) << (bitsPerSymbol * (symbolsInWord - symbol - 1) + residue);
                }
                bitsetArray[n - 1 - i] = word;
            }

            //    fill last
            for (int i = n - 1; i < n; ++i) {
                Output word = 0;
                for (int symbol = 0; (n - 1) * symbolsInWord + symbol < size; symbol++) {
                    word |= Output(mapperForward[a[(i * symbolsInWord + symbol)]]) << (bitsPerSymbol * (symbolsInWord - symbol - 1) + residue);
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


        inline void processAntidiagFormula1(int lowerBound, int upperBound, int lEdge, int tEdge,
                                            UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands, UnsignedMachineWord const *aReverse,
                                            UnsignedMachineWord const *b) {

            const int strandsPerWord = sizeof(UnsignedMachineWord) * 8 - 1;

            /**
             * The idea is as follows.
             * to process some antidiagonal that have been built upon Input cells (that contains a batch of strands) we have to
             * process each such square (Input \times Input) in the antidiagonal fashion.
             * While processing each antidiagonal of such square, we need to implement the following logic:
             * in each iteration swap only those bit-strands that  active in iteration &  satisfy the combing condition.
             * This logic can be implemented by several formulas, here is presented one of it.
             * There is 20 operations inside cycle
             */
#pragma omp   for  simd schedule(static)  aligned(lStrands, tStrands:sizeof(UnsignedMachineWord)*8) aligned(aReverse, b:sizeof(UnsignedMachineWord)*8)
            for (int j = lowerBound; j < upperBound; ++j) {

                UnsignedMachineWord lStrandCap, cond, invCond, tStrandShifted;
                //load phase
                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                auto symbolA = aReverse[lEdge + j];
                auto symbolB = b[tEdge + j];

                auto mask = UnsignedMachineWord(1);

                // manual say 256 just for complete
#pragma GCC unroll  256
                for (int shift = strandsPerWord; shift > 0; shift--) {


                    lStrandCap = lStrand >> shift;

                    cond = mask & ((~(((symbolA >> shift)) ^ symbolB)) | (((~(lStrandCap)) & tStrand)));
                    invCond = ~cond;

                    tStrandShifted = tStrand << shift;
                    tStrand = (invCond & tStrand) | (cond & lStrandCap);

                    cond <<= shift;
                    invCond = ~cond;

                    lStrand = (invCond & lStrand) | (cond & tStrandShifted);

                    mask = (mask << 1) | UnsignedMachineWord(1);
                }

                // center
                cond = (~(symbolA ^ symbolB));
                cond = (cond | ((~lStrand) & tStrand));
                invCond = ~cond;
                tStrandShifted = tStrand;
                tStrand = (invCond & tStrand) | (cond & lStrand);
                lStrand = (invCond & lStrand) | (cond & tStrandShifted);

                mask = ~UnsignedMachineWord(0);

                //lower half
#pragma GCC unroll 256
                for (int shift = 1; shift < strandsPerWord + 1; shift++) {
                    mask <<= 1;

                    lStrandCap = lStrand << (shift);
                    cond = ~(((symbolA << ((shift))) ^ symbolB)); // NOT A XOR B = NOT (A XOR B)// ECONOMY

                    cond = mask & (cond | (((~(lStrandCap)) & tStrand)));
                    invCond = ~cond;

                    tStrandShifted = tStrand >> shift;
                    tStrand = (invCond & tStrand) | (cond & lStrandCap);
                    cond >>= shift;
                    invCond = ~cond;

                    lStrand = (invCond & lStrand) | (cond & tStrandShifted);
                }

                lStrands[lEdge + j] = lStrand;
                tStrands[tEdge + j] = tStrand;
            }
        }


        inline void processAntidiagFormula2(int lowerBound, int upperBound, int lEdge, int tEdge,
                                            UnsignedMachineWord *lStrands, UnsignedMachineWord *tStrands, UnsignedMachineWord const *aReverse,
                                            UnsignedMachineWord const *b) {

            const int strandsPerWord = sizeof(UnsignedMachineWord) * 8 - 1;

            /**
             * The idea is as follows.
             * to process some antidiagonal that have been built upon Input cells (that contains a batch of strands) we have to
             * process each such square (Input \times Input) in the antidiagonal fashion.
             * While processing each antidiagonal of such square, we need to implement the following logic:
             * in each iteration swap only those bit-strands that  active in iteration &  satisfy the combing condition.
             * This logic can be implemented by several formulas, here is presented one of it.
             * There is 15 operations inside cycle
             */
#pragma omp   for  simd schedule(static)  aligned(lStrands, tStrands:sizeof(UnsignedMachineWord)*8) aligned(aReverse, b:sizeof(UnsignedMachineWord)*8)
            for (int j = lowerBound; j < upperBound; ++j) {

                UnsignedMachineWord tStrandCap;
                UnsignedMachineWord lStrandCap, cond;

                //load
                auto lStrand = lStrands[lEdge + j];
                auto tStrand = tStrands[tEdge + j];
                auto symbolA = aReverse[lEdge + j];
                auto symbolB = b[tEdge + j];

                auto mask = UnsignedMachineWord(1);

                // manual say 256 just for complete
#pragma GCC unroll  256
                for (int shift = strandsPerWord; shift > 0; shift--) {
                    // 15 operations inside cycle
                    // could be reduced to 14 if we store not a but ~a in memory

                    lStrandCap = lStrand >> shift;
                    tStrandCap = tStrand << shift;

                    cond = ~((symbolA >> shift) ^ symbolB);

                    tStrand = (lStrandCap | (~mask)) & (tStrand | (cond & mask));
                    lStrand = tStrandCap ^ (tStrand << shift) ^ lStrand;

                    mask = (mask << 1) | UnsignedMachineWord(1);
                }

                // center, no shifts
                cond = ~((symbolA ^ symbolB));
                lStrandCap = lStrand;
                tStrandCap = tStrand;

                tStrand = (lStrandCap | (~mask)) & (tStrand | (cond & mask));
                lStrand = tStrandCap ^ (tStrand) ^ lStrand;

                mask = ~UnsignedMachineWord(0);

                //lower half
#pragma GCC unroll 256
                for (int shift = 1; shift < strandsPerWord + 1; shift++) {
                    mask <<= 1;

                    lStrandCap = lStrand << shift;
                    tStrandCap = tStrand >> shift;
                    cond = ~(((symbolA << (shift)) ^ symbolB));
                    tStrand = (lStrandCap | (~mask)) & (tStrand | (cond & mask));
                    lStrand = tStrandCap ^ (tStrand >> shift) ^ lStrand;
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
            auto[forwardMapper, inverseMapper] = converter.encodeAlphabet(alphabet);
            auto[packedRevA, packedSizeA] = converter.packMSBReverse<InputType, UnsignedMachineWord>(a, aSize, forwardMapper, alphabet.size());
            auto[packedB, packedSizeB] = converter.packSequenceLSB<InputType, UnsignedMachineWord>(b, bSize, forwardMapper, alphabet.size());
            int res = 0;
            if constexpr(IS_ONLY_BINARY_ALPHABET) {
                if constexpr(USE_NAIVE_LAYOUT) {
                    res = llcs2SymbolNaiveCombing(packedRevA, packedSizeA, packedB, packedSizeB, aSize, threadsNum);
                } else {
                    res = llcs2SymbolSmartCombing<USE_FAST_FORMULA>(packedRevA, packedSizeA, packedB, packedSizeB, aSize, threadsNum);
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
        int llcs2SymbolNaiveCombing(UnsignedMachineWord *aReverse, int aSize, UnsignedMachineWord *b, int bSize, int aTotalSymbols, int thdsNum) {
            // also stored in the reverse order
            auto *lStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), sizeof(UnsignedMachineWord) * aSize));
            auto *tStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), sizeof(UnsignedMachineWord) * bSize));

            auto m = aSize, n = bSize;

            // total amount of strands that at the end hit right border of grid
            int disBraid = 0;

            auto numDiag = m + n - 1;

            auto totalSameLengthDiag = numDiag - (m) - (m - 1);

#pragma omp parallel num_threads(thdsNum)  default(none) shared(lStrands, tStrands, aReverse, b, m, n, disBraid, totalSameLengthDiag)
            {
                initStep(lStrands, m, tStrands, n);


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

                auto countOnes = this->countOnes(lStrands, m);
#pragma omp atomic update
                disBraid += countOnes;
            }
            delete[] lStrands;
            delete[] tStrands;

            return aTotalSymbols - disBraid;
        }


        template<bool USE_FAST_COND,bool CONSIDER_BOUNDS_CASE = false>
        int llcs2SymbolSmartCombing(UnsignedMachineWord *aReverse, int aSize, UnsignedMachineWord *b, int bSize, int aTotalSymbols, int thdsNum) {

            // also stored in the reverse order
            auto *lStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), sizeof(UnsignedMachineWord) * aSize));
            auto *tStrands = static_cast<UnsignedMachineWord *> (aligned_alloc(sizeof(UnsignedMachineWord), sizeof(UnsignedMachineWord) * bSize));

            auto m = aSize, n = bSize;

            // total amount of strands that at the end hit right border of grid
            int disBraid = 0;

            auto numDiag = m + n - 1;

            auto totalSameLengthDiag = numDiag - (m - 1) - (m - 1);
#pragma omp parallel num_threads(thdsNum)  default(none) shared(lStrands, tStrands, aReverse, b, m, n, disBraid, totalSameLengthDiag, std::cout)
            {
                initStep(lStrands, m, tStrands, n);
                // phase 1: process upper left triangle
                for (int diagLen = 0; diagLen < m - 1; diagLen++) {
                    if constexpr(!USE_FAST_COND) {
                        this->processAntidiagFormula1(0, diagLen + 1, m - 1 - diagLen, 0, lStrands, tStrands, aReverse, b);
                    } else {
                        this->processAntidiagFormula2(0, diagLen + 1, m - 1 - diagLen, 0, lStrands, tStrands, aReverse, b);
                    }
                }

                // phase2: process parallelogram
                for (int k = 0; k < totalSameLengthDiag; k++) {
                    if constexpr(!USE_FAST_COND) {
                        this->processAntidiagFormula1(0, m, 0, k, lStrands, tStrands, aReverse, b);
                    } else {
                        this->processAntidiagFormula2(0, m, 0, k, lStrands, tStrands, aReverse, b);
                    }
                }

                auto startJ = totalSameLengthDiag;

                // phase:3: lower-right triangle
                for (int diagLen = m - 1; diagLen >= 1; diagLen--) {
                    if constexpr(!USE_FAST_COND) {
                        this->processAntidiagFormula1(0, diagLen, 0, startJ, lStrands, tStrands, aReverse, b);
                    } else {
                        this->processAntidiagFormula2(0, diagLen, 0, startJ, lStrands, tStrands, aReverse, b);
                    }
                    startJ++;
                }

                auto countOnes = this->countOnes(lStrands, m);
#pragma omp atomic update
                disBraid += countOnes;
            }

            delete[] lStrands;
            delete[] tStrands;
            return aTotalSymbols - disBraid;
        }

        void initStep(UnsignedMachineWord *lStrands, int m, UnsignedMachineWord *tStrands, int n) {
            // Initialization step, strands that hit the left grid border all have number 1; strands that hit top grid border are 0.
#pragma omp  for simd schedule(static) aligned(lStrands:sizeof(UnsignedMachineWord)*8)
            for (int k = 0; k < m; ++k) lStrands[k] = ~UnsignedMachineWord(0);
#pragma omp  for simd schedule(static) aligned(tStrands:sizeof(UnsignedMachineWord)*8)
            for (int k = 0; k < n; ++k) tStrands[k] = UnsignedMachineWord(0);

        }

        int countOnes(UnsignedMachineWord *lStrands, int m) {
            auto ones = 0;
            int localSum = 0;
#pragma omp for  simd schedule(static)  aligned(lStrands:sizeof(UnsignedMachineWord)*8)
            for (int i1 = 0; i1 < m; ++i1) {
                //  Brian Kernighan’s Algorithm
                int counter = 0;
                UnsignedMachineWord number = lStrands[i1];
                //  LogNumber
                while (number) {
                    number &= (number - 1);
                    counter++;
                }
                localSum += counter;
            }
            return localSum;
        }


        int llcs_nary_symbol_smart_combing(UnsignedMachineWord *a_reverse, int a_size, UnsignedMachineWord *b, int b_size,
                                           int a_total_symbols, int bits_per_symbol, int threads_num = 1) {

            std::cout << a_size << ',' << b_size << std::endl;
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
                    processAntidiagonal(0, diag_len + 1, m - 1 - diag_len, 0, l_strands, t_strands, a_reverse, b,
                                        residue, bits_per_symbol, braid_ones);
                }

                // phase2: process parallelogram
                for (int k = 0; k < total_same_length_diag; k++) {
                    processAntidiagonal(0, m, 0, k, l_strands, t_strands, a_reverse, b, residue, bits_per_symbol, braid_ones);
                }

                auto start_j = total_same_length_diag;

                // phase:3: lower-right triangle
                for (int diag_len = m - 1; diag_len >= 1; diag_len--) {
                    processAntidiagonal(0, diag_len, 0, start_j, l_strands, t_strands, a_reverse, b,
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