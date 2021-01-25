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
#include "semi_local.h"

std::pair<int,std::pair<std::string,std::string>> parse_input_file(const std::string&  filename){
    std::ifstream file(filename);
    std::string name;
    std::string sequence;
    std::string length;
    std::getline(file, name);
    std::getline(file,length);
    std::getline(file, sequence);
    return std::make_pair(stoi(length), std::make_pair(name,sequence));
};

const long long R = 4294967279;
const long long M = 4294967291;

//static const int length = 1024*8;
long long hash(Permutation &arr, int size) {
    long long hash = 0;
    for (int i =0; i<size;i++){
        hash = (R * hash + arr.get_row_by_col(i)) % M;
    }
    return hash;
}

int * split(std::string str, const std::string& token, int arr_length){
    int *result =  new int[arr_length];
    int i=0;
    while(str.size()){
        int index = str.find(token);
        if(index!=std::string::npos){
            result[i] = stoi(str.substr(0,index));
            i++;
            str = str.substr(index+token.size());
        }else{
            str = "";
        }
    }
    return result;
}



#endif //CPU_FASTA_PARSER_H
