datasheet_dir=/home/noyes046/jsettle/argmobrich/analysis/datasheets
align_script=${datasheet_dir}/run_loman_align.sh
datasheet_script=${datasheet_dir}/gen_datasheet.py
tmp_out=${datasheet_dir}/loman_tmp
data_out=${datasheet_dir}/loman_data

#Align and generate datasheets
#${align_script} "/home/noyes046/jsettle/loman_mock/ERR3152364_1.fastq" "ERR3152364_1"
python ${datasheet_script} "${tmp_out}/ERR3152364_1_megares_mapped.sam" "${tmp_out}/ERR3152364_1_unmapped" "${data_out}"
