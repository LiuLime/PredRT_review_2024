"""
Save pubchem search result and classyfire result as MySQL file

@ Liu
"""

import pymysql


class sql:
    """
    Example
    -------
    >>> sql = sql()
    >>> sql.init()
    >>> cid = sql.search("synomys", "synomys","cid","ribitol")[0]
    >>> name = sql.search("cid","cid","InChI",cid)[0]
    >>> sql.close()
    """

    def __init__(self, db="pubchem"):
        self.db = db

    def __enter__(self):
        self.connection = pymysql.Connection(host="localhost",
                                             user="root",
                                             db=self.db,
                                             charset='utf8mb4')
        self.cur = self.connection.cursor()
        self.cur.execute(f"USE {self.db}")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.cur.close()
        self.connection.close()

    def store(self, **kwargs):
        """

        :param kwargs: cid,sid,InChI,IUPACName,synomys,IsomericSMILES,MolecularFormula
        :return:
        """
        pubchem_cid = kwargs.get("cid")
        pubchem_sid = kwargs.get("sid")
        inchi = kwargs.get("InChI")
        iupac = kwargs.get("IUPACName")
        name = kwargs.get("synomys")
        isomericsmiles = kwargs.get("IsomericSMILES")
        formula = kwargs.get("MolecularFormula")
        sql1 = "INSERT IGNORE INTO cid (cid, sid, InChI, IUPACName, IsomericSMILES, MolecularFormula) VALUES (%s, %s,%s,%s,%s,%s)"
        self.cur.execute(sql1, (pubchem_cid, pubchem_sid, inchi, iupac, isomericsmiles, formula))
        if name:
            name = name.lower()
            sql2 = "INSERT IGNORE INTO synomys (synomys,cid) VALUES (%s,%s)"
            self.cur.execute(sql2, (name, pubchem_cid))
        self.connection.commit()

    def check_exist(self, table, column, query):
        """

        :param namespace: cid, sid, InChI,IUPACName,synomys,IsomericSMILES,MolecularFormula
        :param kwargs:
        :return:
        """

        match column:
            case "synomys":
                query = query.lower()
            case _:
                pass
        sql = f"SELECT * FROM {table} WHERE {column} = %s"
        self.cur.execute(sql, (query))
        if self.cur.rowcount == 0:
            return False
        else:
            return True

    def search(self, table, define_column, query_column, query):
        sql = f"SELECT {query_column} FROM {table} WHERE {define_column} = %s"
        self.cur.execute(sql, (query))
        return self.cur.fetchone()
