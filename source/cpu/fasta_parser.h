////
//// Created by nikita on 30.07.2020.
////
//
#ifndef CPU_FASTA_PARSER_H
#define CPU_FASTA_PARSER_H


#include <cmath>
#include <algorithm>
#include <unordered_set>


#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>

std::pair<int,std::pair<std::string,std::string>> parse_input_file(const std::string&  filename){
    std::ifstream file(filename);
    std::string name;
    std::string sequence;
    std::getline(file, name);
    std::getline(file, sequence);
    return std::make_pair(0,std::make_pair(name,sequence));

};



std::vector<int> transform_to_int_vector(std::string & sequence){
    std::unordered_map<char,int> acgt_to_int( {{'a',0},{'c',1},{'g',2},{'t',3}});
    std::unordered_map<char,int> int_to_acgt( {{0,'a'},{1,'c'},{2,'g'},{3,'t'}});

    std::for_each(sequence.begin(), sequence.end(), [](char & c){
        c = ::tolower(c);
    });

    auto vector = new std::vector<int>();
    for (char & i : sequence) {
        vector->push_back(acgt_to_int.count(i)==0 ? 0 : acgt_to_int[i]);
    }

    return *vector;
};



#endif //CPU_FASTA_PARSER_H
