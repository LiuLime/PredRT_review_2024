"""Figure 2
Detect abnormal data

output:


@Liu 2023/12/21
"""
import os
import sys

sys.path.append("..")
sys.path.append("../02_collate_dataset/")
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from utils import _nameConverter, MySQL

import time
from _BaseAbnormalDetector import abnormalType
import argparse


def get_parser():
    parser = argparse.ArgumentParser(description="detect abnormal data")

    parser.add_argument(
        "--load_ckp",
        default=False,
        help="restart from checkpoint",
        action="store_true",
    )
    parser.add_argument(
        "--start_point",
        default=0,
        help="restart from checkpoint",
        type=int
    )
    parser.add_argument(
        "--load_file_name",
        default="",
        help="the file name for loading ckp, option [smrt],[mbk],[mona],[predret]",
        type=str
    )
    return parser.parse_args()


def read_file(dataset, save_path=None):
    with open(os.path.join(save_path, dataset), "r") as file:
        data = pd.read_csv(file, header=0, sep=",")
        data = data.replace(["N/A", "NA", "na", "n/a", "null", np.nan], None)
    return data


class process:
    def __init__(self, df: pd.DataFrame, name: str, save_path, args):
        self.name = name
        self.df = df
        self.df_dict = df.to_dict(orient='records')
        self.save_path = save_path
        self.args = args
        self.n = self.args.start_point
        self.converter = _nameConverter.FetchInchiFromDatabase()

    def save_checkpoint(self):
        dataframe = pd.DataFrame(self.df_dict)
        dataframe.to_csv(os.path.join(self.save_path, f"{self.name}_check_temp.csv"), index=False, sep=",")

    def df_query_to_mysql(self, namespace):

        match namespace:
            case "name":
                query_list = self.df["C@name"].dropna().unique().tolist()
            case "pubchem":
                query_list = self.df["C@pubchem"].dropna().unique().tolist()
            case _:
                raise KeyError("not valid namespace")

        print(f"{len(query_list)} queries waiting...")
        for idx, query in enumerate(query_list):
            self.converter.convert_identifier_to_inchi(query=namespace, content=query)
            if idx % 50 == 0:
                print(f"processed idx {idx}")

    def df_process(self):
        count = self.n
        for idx, entry in enumerate(self.df_dict):
            if idx < self.n:
                continue

            A = abnormalType(entry)
            # check invalid molecular
            entry["Abnormal@invalid_molecular"] = A.invalid_molecular()

            # check inconsistent object
            # errortype, origin_inchi, smiles_inchi, database_inchi, origin_iupac, database_iupac = A.inconsistent_object()
            errortype, origin_inchi, smiles_inchi, database_inchi = A.inconsistent_object()
            entry["Abnormal@inconsistent"] = errortype
            entry["Abnormal@inconsistent_origin_inchi"] = origin_inchi
            entry["Abnormal@inconsistent_smiles_inchi"] = smiles_inchi
            entry["Abnormal@inconsistent_db_inchi"] = database_inchi
            # entry["Abnormal@inconsistent_origin_inchi_iupac"] = origin_iupac
            # entry["Abnormal@inconsistent_db_inchi_iupac"] = database_iupac

            count += 1
            if count % 100 == 0:
                self.save_checkpoint()
                print(f"save checkpoint at {count}, restart point at {count}")

        final_df = pd.DataFrame(self.df_dict)
        final_df.to_csv(os.path.join(save_path, f"{self.name}_check.csv"), index=False, sep=",")


#
# class smrtAbnormal(abnormalType):
#     def __init__(self, entry: pd.Series):
#         super().__init__(entry)
#         self.entry = entry
#
#     def revise_invalid_molecular(self):
#         return pcp.Compound.from_cid(self.entry["C@pubchem"]).inchi[0]
#
#     def invalid_molecular(self):
#         # have inchi no inchikey-> inchi to mol invalid
#         message = "invalid inchi to molecular"
#         if self.inchi and not self.inchikey:
#             mol = Chem.MolFromInchi(self.inchi)
#             if not mol:
#                 revise_inchi = self.revise_invalid_molecular()
#                 return message, revise_inchi
#         return None, None
#
#     @staticmethod
#     def _fetch_pubchem(id):
#         pubchem_inchi, pubchem_iupac = None, None
#         try:
#             time.sleep(0.5)
#             c = pcp.get_properties(properties=["InChI", "IUPACName"], identifier=id, namespace="cid")[0]
#             if c is None:
#                 print(f"id get null property")
#             else:
#                 pubchem_inchi = c.get("InChI")
#                 pubchem_iupac = c.get("IUPACName")
#         except Exception as e:
#             print(f"isomer pubchem blast fail for {e}")
#         return pubchem_inchi, pubchem_iupac
#
#     @staticmethod
#     def indistinguishable_isomer_represent(df):
#         # same inchi but different pubchem id and corresponding IUPAC name, pubchem blast inchi
#         duplicates = df[df.duplicated('C@inchi', keep=False)]
#         for idx, record in duplicates.iterrows():
#             pubchem_inchi, pubchem_iupac = smrtAbnormal._fetch_pubchem(record["C@pubchem"])
#             duplicates.loc[idx, "Abnormal@isomer_PubchemInChI"] = pubchem_inchi
#             duplicates.loc[idx, "Abnormal@isomer_PubchemIUAPCName"] = pubchem_iupac
#
#         duplicate_temp = duplicates.loc[:, ["C@pubchem",
#                                             "Abnormal@isomer_PubchemInChI",
#                                             "Abnormal@isomer_PubchemIUAPCName"]]
#         df["Abnormal@indistinguishable_isomer"] = df.index.isin(duplicates.index)
#         df = df.merge(duplicate_temp, on="C@pubchem", how="left")
#         return df


# read files
dataset_path = "../02_collate_dataset/collate_data/"
save_path = "./03_result/"
dataset = ["01_SMRT_202309_LC.csv", "02_Massbank_202311.csv", "03_MoNA_202311.csv", "04_predret_202311.csv"]
set_name = ["SMRT", "MassBank", "MoNA", "Predret"]
smrt, mbk, mona, predret = [read_file(dataset[i], dataset_path) for i in range(len(dataset))]
args = get_parser()
if args.load_ckp:
    dataset_path = "./03_result/"
    load_file = args.load_file_name
    if load_file == "smrt":
        smrt = read_file("SMRT_check_temp.csv", dataset_path)
    if load_file == "mbk":
        mbk = read_file("MBK_check_temp.csv", dataset_path)
    if load_file == "mona":
        mona = read_file("MoNA_check_temp.csv", dataset_path)
    if load_file == "predret":
        predret = read_file("predret_check_temp.csv", dataset_path)

# smrt
# print("SMRT starts process----------------------------")

# smrt_invalid, revise_inchi = [], []
#
# for idx, entry in smrt.iterrows():
#     A = smrtAbnormal(entry)
#     # check invalid molecular
#     invalid_ans, inchi_ans = A.invalid_molecular()
#     smrt_invalid.append(invalid_ans)
#     revise_inchi.append(inchi_ans)
#
# smrt["Abnormal@invalid_molecular"] = smrt_invalid
# smrt["Abnormal@invalid_molecular_PubchemInChI"] = revise_inchi
# A = smrtAbnormal(pd.Series())
# smrt = A.indistinguishable_isomer_represent(smrt)
# smrt.to_csv(os.path.join(save_path, "01_SMRT_check2.csv"), index=False, sep=",")


# Predret
print("Predret starts process-----------------------------")
p = process(predret, "predret", save_path, args)
# p.df_query_to_mysql(namespace="name")
# p.df_query_to_mysql(namespace="pubchem")
p.df_process()
print("Predret finish process-----------------------------")
