#!/bin/sh

readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly loman_tmp=${argmobrich}/analysis/datasheets/loman_tmp
readonly loman_fqs=${argmobrich}/loman

readonly err3152364_fq=${loman_fqs}/ERR3152364_1
readonly err3152364_sam=${loman_tmp}/ERR3152364_1
readonly err3152364_parts="1"
readonly err3152364=(${err3152364_fq} ${err3152364_parts} ${err3152364_sam})

readonly err3152365_file=${loman_fqs}/ERR3152365_1
readonly err3152365_sam=${loman_tmp}/ERR3152365_1
readonly err3152365_parts="29"
readonly err3152365=(${err3152365_file} ${err3152365_parts} ${err3152365_sam})

readonly err3152366_fq=${loman_fqs}/ERR3152366_1
readonly err3152366_sam=${loman_tmp}/ERR3152366_1
readonly err3152366_parts="4"
readonly err3152366=(${err3152366_fq} ${err3152366_parts} ${err3152366_sam})

readonly err3152367_file=${loman_fqs}/ERR3152367_1
readonly err3152367_sam=${loman_tmp}/ERR3152367_1
readonly err3152367_parts="30"
readonly err3152367=(${err3152367_file} ${err3152367_parts} ${err3152367_sam})

#readonly loman_samples=(${err3152364[@]} ${err3152365[@]} ${err3152366[@]} ${err3152367[@]})
#readonly indices=(0 2 4 6)
readonly loman_samples=(${err3152365[@]} ${err3152367[@]})
readonly indices=(0 3)

for i in ${indices[@]}; do
    j=$(($i+1))
    k=$(($i+2))
    python loman_gen_on_target.py "${loman_samples[$i]}" "${loman_samples[$j]}" \
                                  "${loman_samples[$k]}"
done
