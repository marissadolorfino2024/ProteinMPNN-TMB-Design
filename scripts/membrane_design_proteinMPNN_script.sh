#!/bin/bash

source /data/brussel/vo/000/bvo00014/venvs/mlfold/bin/activate

folder_with_pdbs="fulldes_bbs_pdb"

output_dir="parsed_chains_temp0.3"
if [ ! -d $output_dir ]
then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"
chains_to_design="A"
fixed_positions="1 2 3 4 5 6 14 15 16 17 18 19 20 32 33 34 40 41 42 43 44 45 46 47 48 49 50 62 63 64 94 95 96 102 103 104 105 106 107 107 109 110 111 112 122 123"


python /data/brussel/vo/000/bvo00014/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

python /data/brussel/vo/000/bvo00014/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

python /data/brussel/vo/000/bvo00014/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

python /data/brussel/vo/000/bvo00014/ProteinMPNN/protein_mpnn_run.py \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 10 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
