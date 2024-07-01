"""
Detect abnormal data records
@ Liu
"""

import os
import sys

sys.path.append("..")
sys.path.append("../02_collate_dataset/")
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from utils import _nameConverter
import pubchempy as pcp
import time
from utils import sendURL
from urllib.parse import quote


class abnormalType:
    def __init__(self, entry: pd.Series | dict) -> None:
        """
        detect abnormal data records including invalid molecular and inconsistent object
        :param entry: pd.Series, iterrows
        """

        self.entry = entry
        self.inchi = None
        if entry["C@inchi"]:
            if entry["C@inchi"].startswith("InChI="):
                self.inchi = entry["C@inchi"]

        self.inchikey = entry["C@inchikey"]
        self.smiles = entry["C@smiles"]
        self.cas = entry["C@cas"]
        self.pubchem = entry["C@pubchem"]
        self.kegg = entry["C@kegg"]
        self.chebi = entry["C@chebi"]
        self.name = entry["C@name"]
        self.have_origin_inchi = entry["CHECK@have_origin_inchi"]
        self.have_db_id = entry["CHECK@have_origin_db_id"]
        self.have_name = "yes" if self.name else "no"
        self.converter = _nameConverter.Converter()
        self.FD = _nameConverter.FetchInchiFromDatabase()

    # check SMRT, MBK, MoNA, Predret
    def invalid_molecular(self):
        # have inchi no inchikey-> inchi to mol invalid
        # have smiles no inchi no inchikey -> smiles to mol invalid
        if self.inchi and not self.inchikey:
            mol = Chem.MolFromInchi(self.inchi)
            if not mol:
                return "invalid inchi to molecular"
        elif self.smiles and not self.inchi and not self.inchikey:
            mol = Chem.MolFromSmiles(self.smiles)
            if not mol:
                return "invalid smiles to molecular"
        return

    # check Predret, SMRT
    @staticmethod
    def _std_inchi(inchi):
        """remove stereo conformation
        :param: inchi
        :return: standard inchi without conformation info
        """
        if inchi:
            mol = Chem.MolFromInchi(inchi)
            if mol:
                try:
                    smiles = Chem.MolToSmiles(mol, isomericSmiles=False)
                    std_inchi = Chem.MolToInchi(Chem.MolFromSmiles(smiles))
                    return std_inchi
                except Exception:
                    print(f"fail to STD {inchi}")
        return None

    def inconsistent_object(self):
        """check valid inchi with database id blast inchi (inchi and database id both available cases) """

        if not self.inchi:  # check valid inchi records
            # return None, None, None, None, None, None
            return None, None, None, None

        if self.have_origin_inchi == "yes" or self.smiles:  # check existance of chemical identifier(smiles or inchi)
            if self.have_db_id == "yes" or self.have_name == "yes":  # check existance of database identifier(database id or chemical name)
                pass
            else:
                # return None, None, None, None, None, None
                return None, None, None, None
        else:
            # return None, None, None, None, None, None
            return None, None, None, None

        # records with both chemical and databse identifier will be processed

        origin_inchi, smiles_inchi, database_inchi, errortype = None, None, None, None
        # origin_inchi_iupac, origin_database_inchi_iupac = None, None

        if self.inchi and self.have_origin_inchi == "yes":
            origin_inchi = self.inchi
        if self.smiles and not origin_inchi:
            smiles_inchi = self.converter.smilesToinchi(self.smiles)
        if self.have_db_id == "yes" or self.have_name == "yes":
            query_list = {"cas": self.cas, "pubchem": self.pubchem, "kegg": self.kegg, "chebi": self.chebi,
                          "name": self.name}
            for idx, identifier in query_list.items():
                if identifier:
                    database_inchi = self.FD.convert_identifier_to_inchi(idx, identifier)
                    break

        inchi_set = {origin_inchi, smiles_inchi, database_inchi}
        try:
            inchi_set.remove(None)
        except KeyError:
            pass

        if len(inchi_set) <= 1:
            pass
        elif len(inchi_set) == 2:
            if origin_inchi and database_inchi:
                # origin_inchi_iupac = self.FD.inchiToIUPAC(origin_inchi)
                # origin_database_inchi_iupac = self.FD.inchiToIUPAC(database_inchi)
                if self._std_inchi(origin_inchi) != self._std_inchi(database_inchi):
                    errortype = "inconsistent_object: orgin_inchi&database_inchi"
                else:
                    errortype = "inconsistent_stereoconformation: orgin_inchi&database_inchi"
            elif smiles_inchi and database_inchi:
                # origin_inchi_iupac = self.FD.inchiToIUPAC(smiles_inchi)
                # origin_database_inchi_iupac = self.FD.inchiToIUPAC(database_inchi)
                if self._std_inchi(smiles_inchi) != self._std_inchi(database_inchi):
                    errortype = "inconsistent_object: smiles_inchi&database_inchi"
                else:
                    errortype = "inconsistent_stereoconformation: smiles_inchi&database_inchi"

        # return errortype, origin_inchi, smiles_inchi, database_inchi, origin_inchi_iupac, origin_database_inchi_iupac
        return errortype, origin_inchi, smiles_inchi, database_inchi