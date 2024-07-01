"""
Fetch information
-----------------
    * fetch compound classification from Classyfire Batch API
      http://cfb.fiehnlab.ucdavis.edu/

    * convert pubchem, cas, chemspider, chebi, kegg, IUPACname to inchi from pubchempy library and following website:
      https://cactus.nci.nih.gov/chemical/structure
      https://www.kegg.jp/kegg/rest/keggapi.html

Return
-----------------
    class FetchClassyfire:
    * superclass and class of compound in [str|list|pd.Series] format
    * save new classyfire json file including new crawler information if inchikey not exist in avaliable ../crawler/classyfire.json

    class FetchDatabase:
    * fetch inchi id from pubchem, cas, chemspider, chebi, kegg, name (IUPAC Name)

class
------------------
    * class FetchClassyfire
    * class FetchDatabase


@ Liu
"""
import sys
import time

sys.path.append("..")
import os
import subprocess
import pandas as pd
from typing import Any, Tuple, Union, List
import json

from utils import sendURL
from crawler import classyfire


def read_saved_classyfire():
    with open("../config.json", "r") as c:
        config = json.load(c)
    with open(config["classyfire_file"], "r") as c:
        classyfire_info = json.load(c)
    return config, classyfire_info


def save_classyfire(data, name):
    with open(name, "w") as j:
        json.dump(data, j)


class FetchClassyfire:
    """
    Example usage
    ------------------

        Examples
        --------
        >>> inchikey_str = "VFIZBHJTOHUOEK-UHFFFAOYSA-N"
        >>> inchikey_list = ["VFIZBHJTOHUOEK-UHFFFAOYSA-N", "BYGQDRYSEFUKGO-UHFFFAOYSA-N"]
        >>> inchikey_series = pd.Series(["VFIZBHJTOHUOEK-UHFFFAOYSA-N", "BYGQDRYSEFUKGO-UHFFFAOYSA-N"])

        For a single string inchikey

        >>> classifier = FetchClassyfire()
        >>> result = classifier.fetchClassyfireBatch(inchikey_str)

        For a list of inchikeys

        >>> classifier_batch1 = FetchClassyfire()
        >>> result_batch_list = classifier_batch1.fetchClassyfireBatch(inchikey_list)

        For a pd.Series of inchikeys

        >>> classifier_batch2 = FetchClassyfire()
        >>> result_batch_series = classifier_batch2.fetchClassyfireBatch(inchikey_series)
    """

    def __init__(self):

        self.config, self.classyfire = read_saved_classyfire()
        self.entry = None

    def inchikey_avaliable(self, Oneinchikey: str) -> Union[dict, None]:
        """
        check inchikey avaliability in classyfire.json, not exist inchikey will fetch info by get_url in classyfire.py
        """
        if Oneinchikey in self.classyfire:
            self.entry = self.classyfire.get(Oneinchikey)
        else:

            self.entry = classyfire.get_url(Oneinchikey)
            if self.entry:
                self.classyfire[Oneinchikey] = self.entry
                save_classyfire(self.classyfire, self.config["classyfire_file"])
        return self.entry

    def _fetch_superclass(self) -> Union[str, None]:
        return self.entry.get("superclass")

    def _fetch_class(self) -> Union[str, None]:
        return self.entry.get("class")

    def fetchClassyfire(self, inchikey) -> Tuple[Union[str, None], Union[str, None]]:
        match inchikey:
            case str():
                self.entry = self.inchikey_avaliable(inchikey)
                if self.entry:
                    return self._fetch_superclass(), self._fetch_class()
                else:
                    return None, None
            case None:
                return None, None
            case _:
                raise TypeError("fetchClassyfire is for single string of inchikey.")

    def fetchClassyfireBatch(self, inchikey) -> Tuple[Union[str, pd.Series, None], Union[str, pd.Series, None]]:
        match inchikey:
            case list():
                inchikey = pd.Series(inchikey)
                return inchikey.apply(self.fetchClassyfire).tolist()
            case pd.Series():
                return inchikey.apply(self.fetchClassyfire)
            case _:
                raise TypeError(
                    "fetch_classify_batch is for list or pd.Series of inchikeys, string see 'fetchClassyfireBatch'.")


