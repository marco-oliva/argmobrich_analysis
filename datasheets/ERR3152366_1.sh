datasheet_dir=/home/noyes046/jsettle/argmobrich/analysis/datasheets
align_script=${datasheet_dir}/run_loman_align.sh
datasheet_script=${datasheet_dir}/gen_datasheet.py
tmp_out=${datasheet_dir}/loman_tmp
data_out=${datasheet_dir}/loman_data

#Align and generate datasheets for single digit indices (because of file naming)
indices=(01 02 03 04)
for i in "${indices[@]}"; do
    ${align_script} "/home/noyes046/jsettle/loman_mock/ERR3152366_1.part-${i}.fastq" "ERR3152366_1.part-${i}"
    python ${datasheet_script} "${tmp_out}/ERR3152366_1.part-${i}_megares_mapped.sam" "${tmp_out}/ERR3152366_1.part-${i}_unmapped" "${data_out}"
done
