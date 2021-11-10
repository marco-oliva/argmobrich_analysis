#!/bin/sh
readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly mobilome=${argmobrich}/analysis/mobilome

readonly err3152366_1="ERR3152366_1"
readonly err3152366_1_num="4"

readonly err3152367_1="ERR3152367_1"
readonly err3152367_1_num="30"

readonly loman_input=(${err3152366_1} ${err3152366_1_num} ${err3152367_1} ${err3152367_1_num})
readonly indices=(0 2)

for i in ${indices[@]}; do
    j=$(($i+1))
    python loman_gen_mobilome.py "${loman_input[$i]}" "${loman_input[$j]}"
done
