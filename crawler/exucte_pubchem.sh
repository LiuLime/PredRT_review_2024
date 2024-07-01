# nohup python3 pubchem.py > "output.log" 2>&1 &
# nohup python3 pubchem.py --load_ckp > "output.log" 2>&1 &
nohup python3 pubchem.py --load_ckp --data_path /Users/liuyuting/Git_Liu/RT_review_2023/01_raw_data/predret_raw.csv --query_column inchi --query_namespace inchi --output_name pubchem_query02_predret.json --load_ckp_file /Users/liuyuting/Git_Liu/RT_review_2023/02_collate_dataset/collate_data/pubchem_query02_predret.json> "output_predret_pubchem_1218.log" 2>&1 &
# nohup python3 pubchem.py --load_ckp --data_path /Users/liuyuting/Git_Liu/RT_review_2023/02_collate_dataset/collate_data/03_MoNA_202311_LC_craw.csv --query_column C@inchi --query_namespace inchi --output_name pubchem_query03_MoNA.json --load_ckp_file /Users/liuyuting/Git_Liu/RT_review_2023/02_collate_dataset/collate_data/pubchem_query02_predret.json> "output_MoNA_pubchem_1215.log" 2>&1 &
