#!/bin/bash

## braids
#python3 braid_mul_testing_system.py 5
#
#
#
## run for figures 4
#python3 fig4.py Datasets/clean_real_dataset/ fig4_real.csv 5 True
#python3 fig4.py Datasets/synthethic_dataset/binary/ fig4_synthethic_2.csv 5 False
#python3 fig4.py Datasets/synthethic_dataset/5/ fig4_synthethic_5.csv 5 False
#python3 fig4.py Datasets/synthethic_dataset/26/ fig4_synthethic_26.csv 5 False
#
## run for figures 5
#
#python3 fig5.py Datasets/clean_real_dataset/ fig5_real.csv 5 True
#python3 fig5.py Datasets/synthethic_dataset/binary/ fig5_synthethic_2.csv 5 False
#python3 fig5.py Datasets/synthethic_dataset/5/ fig5_synthethic_5.csv 5 False
#python3 fig5.py Datasets/synthethic_dataset/26/ fig5_synthethic_26.csv 5 False
#
#
## run for figures 6 and 7
#python3 fig67.py 16 Datasets/clean_real_dataset/ fig67_real.csv  5 True
#python3 fig67.py 16 Datasets/synthethic_dataset/binary/ fig67_synthethic_2.csv 5 False
#python3 fig67.py 16 Datasets/synthethic_dataset/5/ fig67_synthethic_5.csv 5  False
#python3 fig67.py 16 Datasets/synthethic_dataset/26/ fig67_synthethic_26.csv 5 False

# run for figures 8910 synthehic binary  bit parallel old vs new, formula 1 and 2 and iterative hybrid
python3 fig8.py 16 Datasets/synthethic_dataset/bit_binary/ fig8910.csv 6

python3 fig9.py 16 Datasets/synthethic_dataset/bit_binary/ fig9.csv 2
# figures 9:  new formula 1 new formula 2 fot bit just sequential
# figures 10: comparison of parallel semi-local best vs bit parallel for optimal





# precalc stuff compute later
# multiplication compute later
