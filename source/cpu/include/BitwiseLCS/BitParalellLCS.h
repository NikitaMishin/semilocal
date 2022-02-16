
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
            auto bitsPerSymbol = int(std::ceil(log2(alphabetSize)));
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
        std::pair<Output *, int> packMSBReverse(const Input *a, int size, const Map<Input> &mapperForward, int alphabetSize) {

            auto bitsPerSymbol = int(std::ceil(log2(alphabetSize)));
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
     * @tparam UnderlyingMachineWord  should be unsigned int
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


}