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

std::pair<int,std::pair<std::string,std::string>> parse_input_file(const std::string&  filename){
    std::ifstream file(filename);
    std::string name;
    std::string sequence;
    std::getline(file, name);
    std::getline(file, sequence);
    return std::make_pair(0,std::make_pair(name,sequence));

//    std::ifstream input(filename);
//    if(!input.good()){
//        std::cerr << "Error opening "<<filename<<std::endl;
//        return std::make_pair(-1,std::make_pair("",""));
//    }
//
//    std::string line, name, content;
//
//    while( std::getline( input, line ).good() ){
//        if( line.empty() || line[0] == '>' ){ // Identifier marker
//            if( !name.empty() ){ // Print out what we read from the last entry
//                std::cout << name << " : " << content << std::endl;
//                name.clear();
//            }
//            if( !line.empty() ){
//                name = line.substr(1);
//            }
//            content.clear();
//        } else if( !name.empty() ){
//            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
//                name.clear();
//                content.clear();
//            } else {
//                content += line;
//            }
//        }
//    }
//    return std::make_pair(0,std::make_pair(name,content));
};

std::unordered_map<char,int> acgt_to_int( {{'a',0},{'c',1},{'g',2},{'t',3}});
std::unordered_map<char,int> int_to_acgt( {{0,'a'},{1,'c'},{2,'g'},{3,'t'}});


std::vector<int> transform_to_int_vector(std::string & sequence){

    std::for_each(sequence.begin(), sequence.end(), [](char & c){
        c = ::tolower(c);
    });

    auto vector = new std::vector<int>();
    for (int i = 0; i < sequence.size(); ++i) {
        vector->push_back(acgt_to_int.count(sequence[i]) ? 0 : acgt_to_int[sequence[i]]);
    }

    return *vector;
};



#endif //CPU_FASTA_PARSER_H
