#!/bin/bash



# run for figures 4
python3 fig4.py Datasets/clean_real_dataset/ fig4_real.csv 1 True
python3 fig4.py Datasets/synthethic_dataset/binary/ fig4_synthethic_2.csv 1 False
python3 fig4.py Datasets/synthethic_dataset/5/ fig4_synthethic_5.csv 1 False
python3 fig4.py Datasets/synthethic_dataset/26/ fig4_synthethic_26.csv 1 False

# run for figures 5

python3 fig5.py Datasets/clean_real_dataset/ fig5_real.csv 1 True
python3 fig5.py Datasets/synthethic_dataset/binary/ fig5_synthethic_2.csv 1 False
python3 fig5.py Datasets/synthethic_dataset/5/ fig5_synthethic_5.csv 1 False
python3 fig5.py Datasets/synthethic_dataset/26/ fig5_synthethic_26.csv 1 False


# run for figures 6 and 7
python3 fig67.py 4 Datasets/clean_real_dataset/ fig67_real.csv  1 True
python3 fig67.py 4 Datasets/synthethic_dataset/binary/ fig67_synthethic_2.csv 1 False
python3 fig67.py 4 Datasets/synthethic_dataset/5/ fig67_synthethic_5.csv 1  False
python3 fig67.py 4 Datasets/synthethic_dataset/26/ fig67_synthethic_26.csv 1 False

# run for figures 8 synthehic binary  bit parallel old vs new
python3 fig8bit.py 4 Datasets/synthethic_dataset/bit_binary/ fig8_real.csv 1 False

# figures 9:  new formula 1 new formula 2 fot bit just sequential
python3 fig9bit.py 1 Datasets/synthethic_dataset/bit_binary/ fig9_real.csv 1 False

# figures 10: comparison of parallel semi-local best vs bit parallel for optimal TOOD
#python3 fig10bit.py




# precalc stuff compute later
# multiplication compute later
