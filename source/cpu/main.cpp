//
// Created by nikita on 27.07.2020.
//

#include "library.h"
#include <string>
#include <iostream>

int main(){
    std::string str_a = "Hellokjsdfnjsdnf world!";
    std::string str_b = "wrsdhfbsdfjsdfkjsdjkfkjsdbfksdfod!";
    std::vector<char>a(str_a.begin(),str_a.end());
    std::vector<char>b(str_b.begin(),str_b.end());

    std::cout<<naive_prefix_lcs<char>(a,b)<<std::endl;
    std::cout << prefix_lcs_seqential<char>(a, b) << std::endl;

    return 0;
}

