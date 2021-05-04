#!/bin/sh
readonly user_dir=/home/noyes046/jsettle
readonly argmobrich=${user_dir}/argmobrich
readonly mobilome=${argmobrich}/analysis/mobilome

readonly a01_aclame=${mobilome}/A01_aclame.sam
readonly a01_iceberg=${mobilome}/A01_iceberg.sam
readonly a01_pf=${mobilome}/A01_plasmidfinder.sam

readonly b01_aclame=${mobilome}/B01_aclame.sam
readonly b01_iceberg=${mobilome}/B01_iceberg.sam
readonly b01_pf=${mobilome}/B01_plasmidfinder.sam

readonly c01_aclame=${mobilome}/C01_aclame.sam
readonly c01_iceberg=${mobilome}/C01_iceberg.sam
readonly c01_pf=${mobilome}/C01_plasmidfinder.sam

readonly d01_aclame=${mobilome}/MOCK_D01_aclame.sam
readonly d01_iceberg=${mobilome}/MOCK_D01_iceberg.sam
readonly d01_pf=${mobilome}/MOCK_D01_plasmidfinder.sam

readonly e01_aclame=${mobilome}/MOCK_E01_aclame.sam
readonly e01_iceberg=${mobilome}/MOCK_E01_iceberg.sam
readonly e01_pf=${mobilome}/MOCK_E01_plasmidfinder.sam

readonly f01_aclame=${mobilome}/MOCK_F01_aclame.sam
readonly f01_iceberg=${mobilome}/MOCK_F01_iceberg.sam
readonly f01_pf=${mobilome}/MOCK_F01_plasmidfinder.sam

readonly pilot_mobilome_sams=(${a01_aclame} ${a01_iceberg} ${a01_pf} ${b01_aclame} ${b01_iceberg} ${b01_pf} ${c01_aclame} ${c01_iceberg} ${c01_pf} ${d01_aclame} ${d01_iceberg} ${d01_pf} ${e01_aclame} ${e01_iceberg} ${e01_pf} ${f01_aclame} ${f01_iceberg} ${f01_pf} )
readonly indices=(0 3 6 9 12 15)

for i in ${indices[@]}; do
    j=$(($i+1))
    k=$(($i+2))
    python gen_mobilome.py "${pilot_mobilome_sams[$i]}" "${pilot_mobilome_sams[$j]}" "${pilot_mobilome_sams[$k]}"
done
