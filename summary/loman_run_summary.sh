#!/bin/sh

readonly user_dir=/home/noyes046/jsettle
readonly loman=${user_dir}/loman_mock

readonly err3152364_file=${loman}/ERR3152364_1
readonly err3152364_parts="1"
readonly err3152364=(${err3152364_file} ${err3152364_parts})

readonly err3152365_file=${loman}/ERR3152365_1
readonly err3152365_parts="29"
readonly err3152365=(${err3152365_file} ${err3152365_parts})

readonly err3152366_file=${loman}/ERR3152366_1
readonly err3152366_parts="4"
readonly err3152366=(${err3152366_file} ${err3152366_parts})

readonly err3152367_file=${loman}/ERR3152367_1
readonly err3152367_parts="30"
readonly err3152367=(${err3152367_file} ${err3152367_parts})

readonly loman_samples=(${err3152364[@]} ${err3152365[@]} ${err3152366[@]} ${err3152367[@]})
readonly indices=(0 2 4 6)

for i in ${indices[@]}; do
    j=$(($i+1))
    python loman_gen_summary.py "${loman_samples[$i]}" "${loman_samples[$j]}"
done
