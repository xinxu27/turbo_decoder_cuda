cat 12_4Blocks_Bian_Max_10000Fs_Iter25_06_27.txt | grep 'Ber' | cut -d '=' -f2 | awk 'NR%25 == 10{print $0}'
