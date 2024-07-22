# PredRT_review_2024
source code for analyzing current states of retention time datasets till 2024

## 1. Datasets source:
- 01 [SMRT](https://doi.org/10.6084/m9.figshare.8038913.)
  Accessed 11 Nov, 2023
  Ref: Domingo-Almenara X, Guijas C, Billings E, Montenegro-Burke JR, Uritboonthai W, Aisporna AE, Chen E, Benton HP, Siuzdak G: The METLIN small molecule dataset for machine learning-based retention time prediction. Nature communications 2019, 10(1):5811.

- 02 Massbank_japan
  "MassBank-data-2023.11"
  Ref: Horai H, Arita M, Kanaya S, Nihei Y, Ikeda T, Suwa K, Ojima Y, Tanaka K, Tanaka S, Aoshima K: MassBank: a public repository for sharing mass spectral data for life sciences. Journal of mass spectrometry 2010, 45(7):703-714.

- 03 [MoNA-MassBank of North America](https://mona.fiehnlab.ucdavis.edu/)
  "MoNA-export-Experimental_Spectra.json", Accessed 11 Nov, 2023.
  
- 04 Predret
  Accessed 11 Nov, 2023
  Ref: Stanstrup J, Neumann S, Vrhovsek U: PredRet: prediction of retention time by direct mapping between multiple chromatographic systems. Analytical chemistry 2015, 87(18):9421-9428.

## 2. Abnormal data detection
### purpose:
`./AbnormalRecord/AbnormalDataDetector.py` is used for detecting inconsistent records including invalid molecular and inconsistent object, the method is by comparing recorded inchi with database identifier/name transformed inchi.    
For predret dataset:
- if they matched, they will not be marked as abnormal
- if they are not match firstly, but could match after removing stereo information, they will be marked as `inconsistent_stereoconformation: orgin_inchi&database_inchi` in output dataframe in `Abnormal@inconsistent_type` column;
- if they are not match even after removing stereo informatino, they will be marked as `inconsistent_object: orgin_inchi&database_inchi` in `Abnormal@inconsistent_type` column;
- The orgin inchi and database identifier searched inchi will be filled in `Abnormal@inconsistent_type_origin_inchi` and `Abnormal@inconsistent_type_database_inchi` columns relatively;
- if inchi id could not be transformed into `molecular` object by rdkit tool, they will be marked as `invalid inchi to molecular` in `Abnormal@invalid_molecular` column.    
For SMRT dataset:
- if inchi id could not be transformed into `molecular` object by rdkit tool, they will be marked as `invalid inchi to molecular` in `Abnormal@invalid_molecular` column;
- if there are multiple duplicated isomer inchi presented, the record will be mared as 'True' in `Abnormal@indistinguishable_isomer` column;
- if it is 'True' in `Abnormal@indistinguishable_isomer` column, inchi blasted result in pubchem will be filled in `Abnormal@indistinguishable_isomer_PubchemInChI` and `Abnormal@indistinguishable_isomer_IUAPCName` columns.
### Usage:
1.please put dataframe for detection in `collate_dataset` folder    
Input Example:     
file name: `SMRT.csv`,`predret.csv`    
file format example:    

| ID@db_id | C@name | C@pubchem | C@chebi | C@kegg | C@cas | C@inchi                                         | C@inchikey                  | CH@rt | CH@rt_second | CH@rt_no_unit | MS@instrument | CHECK@have_origin_inchi | CHECK@have_origin_db_id |
|----------|--------|-----------|---------|--------|-------|-------------------------------------------------|-----------------------------|-------|--------------|---------------|---------------|-------------------------|-------------------------|
| SMRT_0   |        | 5139      |         |        |       | InChI=1S/C3H8N2S/c1-2-6-3(4)5/h2H2,1H3,(H3,4,5) | VFIZBHJTOHUOEK-UHFFFAOYSA-N | 93.5  | 93.5         |               | LC            | yes                     | yes                     |
|          |        |           |         |        |       |                                                 |                             |       |              |               |               |                         |                         |

running:    
first time:
```shell
python3 AbnormalDataDetector.py  
```
running this if process was broke down before finish:    
```shell
python3 AbnormalDataDetector.py  --load_ckp --load_file_name SMRT_check_temp.csv --start_point XX  # XX-> number for start rows
```

output:   `SMRT_check_temp.csv`,`predret_check_temp.csv` in `./AbnormalRecord/03_result/` path

## 3. Clean data analysis
`_Filters.py`       
purpose: Filter data records by inchikey in collate dataframe, classify instrument type (GC or LC),
and count unique compounds by inchikey)    
- screen out records with inchikey (records with valid molecular)
- classify records by instrument type GC and LC
- count unique compound records according to instrument type
- print summary information

`Sorter.py`:    
Input: same dataframe format as mentioned above. Records with valid `C@inchikey` will be selected to processs here.  
Output: Table 1,2,3,4   
Count RT records in SMRT, MassBank, MoNA and Predret datasets    
- Screen out records with inchikey (exclude records with invalid identifier) [_Filters.py]    
- classify RT records into RT records with unit(min or sec) or no unity fill in two columns "CH@rt_second", "CH@rt_no_unit" 
    Notice: min was transformed into sec by multiplying 60
    [_Filters.py]
- classify records into LC (liquid-chromatography) and GC(gas-chromatography) instrument type, summarize records individually.
  [_Filters.py]

output files:
--------------    
- save_path: './04_result/'
- 01_SMRT_LC.csv: final clean data with valid inchi and inchikey
- 02_MBK_LC.csv, 03_MoNA_LC.csv, 04_predret_LC:
  final clean data with valid inchi and inchikey, clear unit type "CH@rt_second" & "CH@rt_no_unit", and LC type only
- 02_MBK_GC.csv, 03_MoNA_GC.csv, 04_predret_GC.csv:
  final clean data with valid inchi and inchikey, clear unit type "CH@rt_second" & "CH@rt_no_unit", and GC type only
- 02_MBK_LC_summary.csv, 03_MoNA_LC_summary.csv, 04_predret_LC_summary.csv:
  Summary of Massbank data source, number of compounds, number of unique compounds using LC instrument type
- 02_MBK_GC_summary.csv, 03_MoNA_GC_summary.csv, 04_predret_GC_summary.csv:
  Summary of Massbank data source, number of compounds, number of unique compounds using GC instrument type

Notice
----------    
If origin records provide indistinguishable inchi for isormers, then their inchi and inchikey will be the same in cleaned data,
please keep this point in mind.

Usage:    
```shell
python3 Sorter.py  
```

`ChemDraw.py`:      
Purpose:    
Fig.2A: Ratio of compound superclass in SMRT, MassBank, MoNA, PredRet datasets (LC records only).   
Fig.2B: Upset figure show intersections across datasets (LC records only).    

input:    
Please conduct  `Sorter.py` first, this script will utilize the output dataframe saved by `Sorter.py`.    

output:    
- save_path: './04_result/'
- Figure2A_pie.csv: statistical analysis result
- Figure2A_pie.png: ratio of compound superclass across datasets
- Figure2B_upset.png

Usage:    
```shell
python3 ChemDraw.py  
```
