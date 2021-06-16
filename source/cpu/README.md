# CPU based implementations for solving semi local lcs problem (proof of concept)

## Implementations

* cpu_impl directory contains run programs that runs several algorithms:

- Input have similar interfaces:
  
    - in case of braid multiplications programs accept several arguments: 
      * int depth how many layers should be paralellized
      * int n  size of matrices that should be multiplied
      * int seed  random seed
    - in case of semi local lcs and lcs  algorithms
      * int thds number of threads 
      * string filepath_a path to file with sequence A
      * string filepath_b path to file with sequence B
         -  both string in the is a integer sequence with delimiter ',' between each symbol see for example directory Dataset with sequences both real and synthetic
- Output:
  * precalc_elapsed_time time for precomputation in ms
  * elapsed_time time of algorithm in ms without precomputation time 
  * hash:
    - in case of lcs -- it simply lcs score
    - in case of semi-local lcs it is hash of permutation matrix
    - in case of braid multiplication - hash of permutation matrix
  * seed used seed for generation  

* prefix_lcs directory contains algorithms for compute lcs score (both simple and bitwise ones)
* semi_local.h contains algorithms for solving semi local lcs
* unit_monge_mult directory contains implementation of steady_ant algorithm
* scirpt.bash allows to run several py files that performs running of tests and form csv with with details of runs
* .py files is a part of testing system for evaluating algorithms
* requirements:

    - G++10.2.0 
    - fopenmp
    - c++11 or c++17 standard

If yoy have any questions or problem feel free to ask via opening new issue or directly emailing to mishinnikitam@protonmail.com since
this repo lack of documentation and under development.

Soon this branch will be refactored and there will be c++ class wrappers around algorithm (also will be gpu implementations) see gpu branch for example of changes
