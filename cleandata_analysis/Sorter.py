"""Table 1,2,3,4
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
please keep this point in mind

@ Liu
"""

import sys

sys.path.append(".")
import os
import pandas as pd
import numpy as np
from _Filters import Filter, print_collate_info


def read_file(dataset, save_path=None):
    with open(os.path.join(save_path, dataset), "r") as file:
        data = pd.read_csv(file, header=0, sep=",")
        data = data.replace(["N/A", "NA", "na", "n/a", "null", np.nan], None)
    return data


def process_dataset(dataset, save_path=None, filename=None, inst=None):
    """process instrument screening and rt unit check"""
    # chose records with inchikey
    F = Filter(dataset)
    dataset = F.keep_valid_molecular(dataset)
    # filter instrument
    dataset_inst = F.filter_df_by_instrument(dataset, inst)
    # unify rt unit
    dataset_rt_inst = F.filter_df_by_rt(dataset_inst)
    # save abstract of dataset
    dataset_rt_inst.to_csv(os.path.join(save_path, f"{filename}_{inst}.csv"), sep=",", index=False)
    print_collate_info(dataset_rt_inst).to_csv(os.path.join(save_path, f"{filename}_{inst}_summary.csv"), sep=",",
                                               index=True)


data_path = "../02_collate_dataset/collate_data/"
save_path = "04_result/"
dataset = ["01_SMRT_202309_LC.csv", "02_Massbank_202311.csv", "03_MoNA_202311.csv", "04_predret_202311.csv"]
smrt, mbk, mona, predret = [read_file(dataset[i], data_path) for i in range(len(dataset))]

"""SMRT dataset"""
print("SMRT dataset process------------------------")
# records with inchikey
F = Filter(smrt)
smrt = F.keep_valid_molecular(smrt)
# print duplicate records
duplicates = smrt[smrt.duplicated('C@inchikey', keep=False)]
print("duplicate records:", duplicates["C@inchikey"], duplicates.shape[0])
print_smrt_df = print_collate_info(smrt)
smrt.to_csv(os.path.join(save_path, "01_SMRT_LC.csv"), sep=",", index=False)

"""massbank dataset"""
print("Massbank dataset process--------------------")
process_dataset(mbk, save_path, filename="02_MBK", inst="lc")
process_dataset(mbk, save_path, filename="02_MBK", inst="gc")

"""MoNA dataset"""
print("MoNA dataset process------------------------")
process_dataset(mona, save_path, filename="03_MoNA", inst="lc")
process_dataset(mona, save_path, filename="03_MoNA", inst="gc")

"""predret dataset"""
print("PredRet dataset process---------------------")
process_dataset(predret, save_path, filename="04_predret", inst="lc")
process_dataset(predret, save_path, filename="04_predret", inst="gc")
