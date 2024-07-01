"""
Convert smiles or inchi into inchikey

Notice
--------

**The stereochemistry is KEPT!**

In open-source database, some records point out the stereo structure in InChI identifier however most of them did
not, one possible reason may cause by uncertain annotation or using a mixture reference standards for isomers.
Considering the importance of stereo information for isomer isolation, we kept the isomeric representation in
original records. As example showed below, with or without isomeric information will generate different InChIKey.
**So, we want to notice that the statistical analysis of unique compounds by inchikey will count R,S or mixture isomers as
different unique compounds.**

Examples
--------
* (1S)-6-Chloro-1-phenyl-2,3,4,5-tetrahydro-1H-3-benzazepine-7,8-diol
InChI=1S/C16H16ClNO2/c17-15-11-6-7-18-9-13(10-4-2-1-3-5-10)12(11)8-14(19)16(15)20/h1-5,8,13,18-20H,6-7,9H2/t13-/m0/s1
c1ccc(cc1)[C@@H]2CNCCc3c2cc(c(c3Cl)O)O
GHWJEDJMOVUXEC-ZDUSSCGKSA-N

* (1R)-6-Chloro-1-phenyl-2,3,4,5-tetrahydro-1H-3-benzazepine-7,8-diol
c1ccc(cc1)[C@H]2CNCCc3c2cc(c(c3Cl)O)O
InChI=1S/C16H16ClNO2/c17-15-11-6-7-18-9-13(10-4-2-1-3-5-10)12(11)8-14(19)16(15)20/h1-5,8,13,18-20H,6-7,9H2/t13-/m1/s1
GHWJEDJMOVUXEC-CYBMUJFWSA-N

* 6-Chloro-1-phenyl-2,3,4,5-tetrahydro-1H-3-benzazepine-7,8-diol
InChI=1S/C16H16ClNO2/c17-15-11-6-7-18-9-13(10-4-2-1-3-5-10)12(11)8-14(19)16(15)20/h1-5,8,13,18-20H,6-7,9H2
c1ccc(cc1)C2CNCCc3c2cc(c(c3Cl)O)O
GHWJEDJMOVUXEC-UHFFFAOYSA-N

@Liu
"""

from rdkit import Chem
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.inchi import rdinchi
import pandas as pd
from typing import Union, List, Tuple
import sys
import warnings

sys.path.append("..")

from utils import sendURL, MySQL

from urllib.parse import quote

import json
from bs4 import BeautifulSoup as BS

sys.path.append("../02_collate_dataset")


class GenerateIdentifier:
    def __init__(self, **kwargs) -> None:
        """
        generating inchi or inchikey dealing with single cases

        :param kwargs:
        """

        self.FD = FetchInchiFromDatabase()
        self.converter = Converter()

        self.idx = kwargs.get("idx")
        self.inchi = kwargs.get("inchi")
        self.smiles = kwargs.get("smiles")
        self.cas = kwargs.get("cas")
        self.pubchem = kwargs.get("pubchem")
        self.kegg = kwargs.get("kegg")
        self.chebi = kwargs.get("chebi")
        self.name = kwargs.get("name")

    def generate_inchi(self) -> str | None:
        """fetch inchi identifier from database via database id or smiles via RDkit"""
        inchi = None
        if self.inchi:
            inchi = self.converter.unify_inchi(self.inchi)
        elif self.smiles:
            inchi = self.converter.smilesToinchi(self.smiles)
        else:  # query database
            query_list = {"cas": self.cas, "pubchem": self.pubchem, "kegg": self.kegg, "chebi": self.chebi,
                          "name": self.name}
            for idx, identifier in query_list.items():
                if identifier:
                    inchi = self.FD.convert_identifier_to_inchi(idx, identifier)
                    break

        if inchi is None:
            message = f"{self.idx} has no valid identifier"
            warnings.warn(message, UserWarning)
        return inchi

    def generate_inchikey_from_inchi(self, ans_inchi) -> str | None:
        """Single convertion query for generating inchikey from inchi

        :returns: InChIkey string by convert InChI
        """
        return self.converter.inchiToinchikey(ans_inchi)


class Converter:
    """
    Convert smiles to inchi, inchi to inchikey by RDKit
    invalid inchi or smiles will raise Userwarnings
    """

    def __init__(self):
        pass

    @staticmethod
    def unify_inchi(inchi):
        if inchi is not None:
            inchi = str(inchi)
            if inchi.startswith("InChI="):
                pass
            elif inchi.startswith("('InChI="):  # deal with abnormal format
                try:
                    inchi_tuple = eval(inchi)  # convert to tuple
                    inchi = inchi_tuple[0]
                except Exception:
                    pass
            else:
                pass
        return inchi

    def check_validity(self, represent, identifier: str) -> bool:
        # print(f"check {identifier} validity")
        try:
            match represent:
                case "inchi":
                    mol = Chem.MolFromInchi(identifier)
                case "smiles":
                    mol = Chem.MolFromSmiles(identifier)
                case _:
                    raise ValueError("please input inchi or smiles")
            if mol:
                return True
            else:
                warnings.warn(f"{identifier} to invalid molecule, return NoneType", UserWarning)
                return False
        except Exception as e:
            print(e)
            return False

    def inchiToinchikey(self, inchi: Union[str, tuple]) -> Union[str, None]:
        inchikey = None
        inchi = self.unify_inchi(inchi)
        match inchi:
            case str():
                if self.check_validity("inchi", inchi):
                    inchikey = rdinchi.InchiToInchiKey(inchi)
            case None:
                pass
            case _:
                print(type(inchi), inchi)
                raise TypeError("inchiToinchikey is for single string, list or pd.Series see inchiToinchikeyBatch.")
        return inchikey

    def inchiToinchikeyBatch(self, inchi_batch: Union[List[str], pd.Series]) \
            -> Union[List[str], pd.Series]:
        match inchi_batch:
            case list():
                return [self.inchiToinchikey(i) for i in inchi_batch]
            case pd.Series():
                return inchi_batch.apply(self.inchiToinchikey)
            case _:
                raise TypeError("inchiToinchikeyBatch is for list or pd.Series, single string see inchiToinchikey.")

    def smilesToinchi(self, smiles: str) -> str | None:
        match smiles:
            case str():
                if self.check_validity("smiles", smiles):
                    mol = Chem.MolFromSmiles(smiles)
                    inchi = rdinchi.MolToInchi(mol)
                    return inchi
            case None:
                pass
            case _:
                print(type(smiles), smiles)
                raise TypeError("smilesToinchi is for single string, list or pd.Series see smilesToinchikeyBatch.")
        return None

    def smilesToinchiBatch(self, smiles_batch: Union[List[str], pd.Series]) \
            -> Union[List[str], pd.Series]:
        match smiles_batch:
            case list():
                return [self.smilesToinchi(s) for s in smiles_batch]
            case pd.Series():
                return smiles_batch.apply(self.smilesToinchi)
            case _:
                raise TypeError("smilesToinchikeyBatch is for list or pd.Series, single string see smilesToinchikey.")


class FetchInchiFromDatabase:
    """
    Retrieve InChI identifier through database id.

    Parameters
    ----------
    query: chose from "pubchem", "chebi","cas", "kegg", "chebi", "iupacname"
    content: database identifier, see 'Content Example'

    Return
    ----------
    InChI string

    Content Example:
    ----------
        PUBCHEM:
            cid:16720 or sid:9860 or 16720 (only input int will recognize as cid identifier)
        CHEBI:
            CHEBI:53003 or 53003
        CAS:
            41394-05-2
        KEGG:
            C10930
        IUPACNAME:
            4-amino-3-methyl-6-phenyl-1,2,4-triazin-5-one

    Examples
    --------
    >>> F = FetchInchiFromDatabase()
    >>> F.open_sql()
    >>> return_inchi = F.convert_identifier_to_inchi("cas", "41394-05-2")
    >>> F.close_sql()
    """

    def __init__(self):
        self.query = None
        self.content = None

        self.kegg_headers = {
            'accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7',
            'accept-language': 'zh-CN,zh;q=0.9,en;q=0.8,ja;q=0.7',
            'cache-control': 'max-age=0',
            'sec-ch-ua': '"Not_A Brand";v="8", "Chromium";v="120", "Google Chrome";v="120"',
            'sec-ch-ua-mobile': '?0',
            'sec-ch-ua-platform': '"macOS"',
            'sec-fetch-dest': 'document',
            'sec-fetch-mode': 'navigate',
            'sec-fetch-site': 'none',
            'sec-fetch-user': '?1',
            'upgrade-insecure-requests': '1',
            'user-agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36',
        }

        self.config = self.open_config()
        self.MySQL = MySQL.sql()

    def open_config(self):
        with open("../config.json", "r") as c:
            config = json.load(c)
        return config

    def convert_identifier_to_inchi(self, query: str, content: str | None) -> str | None:
        """
        Match query of pubchem, chebi, cas, kegg and name (IUPAC Name recommend)

        Parameter:
            * query: str, "pubchem","chebi","cas","kegg","name"
            * content: str, identifier
            * properties: list,the fetched content from pubchem into MySQL, default=["InChI", "IUPACName", "IsomericSMILES", "MolecularFormula"]

        Returns:
            * InChI: str
        """
        self.query = query
        self.content = content
        if not self.content:
            return

        match self.query:
            case "pubchem":
                return self.pubchemToinchi()
            case "chebi":
                return self.chebiToinchi()
            case "cas":
                return self.casToinchi()
            case "kegg":
                return self.keggToinchi()
            case "name":
                return self.nameToinchi()

    def pubchemToinchi(self, namespace=None) -> str | None:
        """
         Queries PubChem based on the given namespace. If the record is not in the local MySQL database,
         it searches PubChem, stores the result in MySQL, and returns the InChI string.

         :param namespace: One of 'cid', 'sid', or 'name' to specify the query type.
         :return: The InChI string or None if not found.
         """

        # deal with CID:XXXX or SID:XXXX content
        if namespace is None:
            namespace, self.content = self.std_content("pubchem")
        # deal with float type pubchem cid
        self.floatToint()
        cid, sid, inchi, iupac, isomericsmiles, formula, synomys = None, None, None, None, None, None, None

        match namespace:
            case "cid":
                table, column = "cid", "cid"
                cid = self.content
            case "sid":
                table, column = "cid", "sid"
                sid = self.content
            case "name":
                table, column = "synomys", "synomys"
                synomys = self.content
            case _:
                raise KeyError("invalid namespace. Use 'cid','sid' or 'name'")

        # if record not exist in mysql db, search in pubchem and store in mysql
        with self.MySQL as sql:
            if not sql.check_exist(table=table, column=column, query=self.content):
                cid, inchi, iupac, isomericsmiles, formula = query_pubchem(namespace=namespace, content=self.content)
                sql.store(cid=cid,
                          InChI=inchi,
                          IUPACName=iupac,
                          IsomericSMILES=isomericsmiles,
                          MolecularFormula=formula,
                          sid=sid,
                          synomys=synomys)
            else:
                try:
                    if namespace == "name":
                        cid = sql.search("synomys", "synomys", "cid", synomys)[0]
                        inchi = sql.search("cid", "cid", "InChI", cid)[0]
                    if namespace == "sid":
                        inchi = sql.search("cid", "sid", "InChI", sid)[0]
                    if namespace == "cid":
                        inchi = sql.search("cid", "cid", "InChI", cid)[0]
                except Exception:
                    print(f"{self.content} InChI is null")
        return inchi

    def casToinchi(self):
        return self.nameToinchi()

    def keggToinchi(self):
        url = f"https://rest.kegg.jp/conv/pubchem/cpd:{self.content}"
        try:
            keggTopubchemSID = sendURL.send_url(url, headers=self.kegg_headers).split("pubchem:")[1]
            print(f"convert KEGG ID {self.content} to PUBCHEM SID {keggTopubchemSID} success")
            self.content = keggTopubchemSID
        except AttributeError:
            print(f"convert KEGG ID {self.content} to PUBCHEM SID failed")
            return None
        return self.pubchemToinchi(namespace="sid")

    def chebiToinchi(self):
        self.content = self.std_content("chebi")
        self.floatToint()
        url = f"https://rest.kegg.jp/conv/compound/chebi:{self.content}"
        try:
            chebiTokegg = sendURL.send_url(url, headers=self.kegg_headers).split("cpd:")[1]
            print(f"convert CHEBI ID {self.content} to KEGG ID {chebiTokegg} success")
            self.content = chebiTokegg
        except AttributeError:
            print(f"convert CHEBI ID {self.content} to KEGG ID failed")
            return None
        return self.keggToinchi()

    def nameToinchi(self) -> str | None:
        """
        Blast name in pubchem
        :return: inchi
        """
        return self.pubchemToinchi(namespace="name")

    def floatToint(self):
        """
        Unify pubchem id and chebi id into INT type
        :return: INT type id
        """
        match self.content:
            case str():  # str [name or "1353.0" or "1353"]
                try:
                    self.content = float(self.content)  # "1353.0" or "1353" case
                    self.content = int(self.content)  # convert to "1353" int
                except ValueError:  # name case
                    pass
            case float():  # float ["1353.0"]
                self.content = int(self.content)  # convert to "1353" int
            case int():  # int ["1353"]
                pass
            case _:  # unknown type
                print("unknown id type", type(self.content))

    def std_content(self, db_type):
        match db_type:
            case "pubchem":
                id = str(self.content).lower()
                if "cid" in id:
                    if ":" in id:
                        return "cid", id.split("cid:")[1]
                    else:
                        return "cid", id.split("cid")[1]
                elif "sid" in id:
                    if ":" in id:
                        return "sid", id.split("sid:")[1]
                    else:
                        return "sid", id.split("sid")[1]
                else:
                    return "cid", id

            case "chebi":
                id = str(self.content).lower()
                if ":" in id:
                    return id.split(":")[1]
                elif "chebi" in id:
                    return id.split("chebi")[1]
                else:
                    return id
            case _:
                return self.content


def query_pubchem(namespace, content, properties=None):
    """

    :param namespace: name, cid,sid
    :param content: query content
    :param properties: default=["InChI", "IUPACName", "IsomericSMILES", "MolecularFormula"]
    :return: Tuple(cid, inchi, iupac, isomericsmiles, formula)
    """
    if properties is None:
        properties = ["InChI", "IUPACName", "IsomericSMILES", "MolecularFormula"]
    query_list = ",".join(i for i in properties)
    url = None
    cid, inchi, iupac, isomericsmiles, formula = None, None, None, None, None

    match namespace:
        case "name":
            content = quote(content)
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{content}/property/{query_list}/XML"
        case "cid":
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{content}/property/{query_list}/XML"
        case "sid":
            url_sid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{content}/cids/XML"
            ans_sid = sendURL.send_url(url_sid, rate_limit=0.3)
            if ans_sid:
                bs_sid = BS(ans_sid, "lxml-xml")
                cid = bs_sid.find("CID").get_text()
                url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{query_list}/XML"
        case _:
            raise KeyError("not valid query namespace")
    if url:
        ans = sendURL.send_url(url, rate_limit=0.3)
        if ans:
            try:
                bs = BS(ans, "lxml-xml")
                cid = bs.find("CID").get_text()
                inchi = bs.find("InChI").get_text()
                iupac = bs.find("IUPACName").get_text()
                isomericsmiles = bs.find("IsomericSMILES").get_text()
                formula = bs.find("MolecularFormula").get_text()
            except Exception as e:
                print(e)
    return cid, inchi, iupac, isomericsmiles, formula

# %%
