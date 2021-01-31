#!/bin/bash

# to draw charts for optimizations
python3 single_threaded_test.py 1 Datasets/synthethic_dataset/binary/ res_binary_single_threaded_test.csv
python3 single_threaded_test.py 1 Datasets/synthethic_dataset/5/ res_5_single_threaded_test.csv
python3 single_threaded_test.py 1 Datasets/synthethic_dataset/26/ res_26_single_threaded_test.csv

# test recursive combing
python3 rec_testing_combing.py  1 Datasets/synthethic_dataset/26/ res_26_rec.csv

# all for braid mult
python3 braid_mul_testing_system.py 7

# test parallel versions
python3 combing_testing_system.py 16 Datasets/synthethic_dataset/binary/ res_par_binary.csv

python3 combing_testing_system.py 16 Datasets/synthethic_dataset/5/ res_par_5.csv

python3 combing_testing_system.py 16 Datasets/synthethic_dataset/26/ res_par_26.csv

# on real data need when find it