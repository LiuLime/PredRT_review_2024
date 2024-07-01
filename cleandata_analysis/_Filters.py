"""
Filter data records by having inchikey in collate dataframe, classify instrument type (GC or LC),
and count unique compounds by inchikey)
- screen out records with inchikey (records with valid molecular)
- classify records by instrument type GC and LC
- count unique compound records according to instrument type
- print summary information

@ Liu
"""

from typing import List
from pandas import DataFrame
import pandas as pd
import numpy as np


class Filter:
    def __init__(self, dataframe: pd.DataFrame):
        self.dataframe = dataframe

    @staticmethod
    def unify_na(df: pd.DataFrame) -> DataFrame:
        """unify loss value format
        :param: dataframe
        :return: dataframe"""
        df = df.replace(["N/A", "NA", "na", "n/a", "null", np.nan], None)
        return df

    @staticmethod
    def check_gc(cell):
        """
        Check instrument type

        :param cell: unit cell in dataframe
        :return: GC (gas chromatography) return 0, LC (liquid chromatography) return 1 and null information return 2
        """
        try:
            cell = cell.lower()
            if "gc" in cell:
                return 0  # GC
            else:
                return 1  # LC
        except (TypeError, AttributeError):
            return 2  # null

    """collate Retention Time records into having time unit class and not having tiome unit class"""

    @staticmethod
    def check_rt(entry) -> int | None:
        """
        check value is number and return unit if exist otherwise return NoneType

        :param entry: value for check
        :return: rt with unit return 0; rt without unit return 1, NoneType and illegal character return 2
        """
        try:
            if entry is np.nan:
                return 2

            entry = str(entry)
            float(entry.split()[0])
            u = entry.split()[1]
            return 0
        except IndexError:
            return 1  # no unit
        except (AttributeError, ValueError):  # NoneType->AttributeError
            return 2

    @staticmethod
    def filter_rt_unit(entry) -> int | None:
        """
        return retention time which unit is min or sec, min will transform to sec unit

        :param entry: value for check
        :return: value in second, int or NoneType
        """
        entry = str(entry)
        unit_bool = Filter.check_rt(entry)
        if unit_bool == 0:
            unit = entry.split()[1].lower()
            value = float(entry.split()[0])
            if "m" in unit:
                return value * 60 if value >= 0 else None  # rt_unit is minutes
            elif "s" in unit:
                return value if value >= 0 else None  # rt_unit is seconds
        elif unit_bool == 1 or unit_bool == 2:
            return None  # rt have no unit or no value

    @staticmethod
    def filter_rt_no_unit(entry) -> int | None:
        """
        return retention time which do not have unit

        :param entry: value for check
        :return: int or NoneType
        """
        entry = str(entry)
        unit_bool = Filter.check_rt(entry)
        if unit_bool == 0 or unit_bool == 2:
            return None
        elif unit_bool == 1:  # rt have no unit
            value = float(entry.split()[0])
            return value if value >= 0 else None

    def keep_valid_molecular(self, df: pd.DataFrame) -> pd.DataFrame:
        """keep record with valid inchikey"""
        df = df[df["C@inchikey"].notna()]
        return df

    def filter_df_by_rt(self, df: pd.DataFrame) -> DataFrame:
        """
        classify dataframe CH@rt column into CH@rt_unit and CH@rt_no_unit columns
        :param: Dataframe
        :return: Dataframe
        """
        dataframe_rt = Filter.unify_na(df)
        # if not dataframe_rt["CH@rt_second"].isna().all():
        dataframe_rt.loc[:, "CH@rt_second"] = dataframe_rt["CH@rt"].apply(Filter.filter_rt_unit)
        # if not dataframe_rt["CH@rt_no_unit"].isna().all():
        dataframe_rt.loc[:, "CH@rt_no_unit"] = dataframe_rt["CH@rt"].apply(
                Filter.filter_rt_no_unit)
        dataframe_rt = dataframe_rt.dropna(how="all", subset=["CH@rt", "CH@rt_second", "CH@rt_no_unit"])
        print(
            f"{dataframe_rt.shape[0]} records after rt screening, {len(dataframe_rt['C@inchikey'].unique())} unique inchikey")
        return dataframe_rt

    def filter_df_by_instrument(self, df: pd.DataFrame, inst: str) -> pd.DataFrame:
        """
        classify dataframe according to instrument type lc(liquid-chromatography) or gc(gas-chromatography)
        :param: df: Dataframe
        :param: inst: "gc" or "lc" str
        :return: Dataframe filtered by instrument lc or gc
        """
        dataframe_inst = Filter.unify_na(df)
        if inst == "gc":
            flag = 0
        elif inst == "lc":
            flag = 1
        else:
            raise ValueError("instrument type lc or gc is not given")
        check_inst = dataframe_inst["MS@instrument"].apply(Filter.check_gc)
        record_to_keep = check_inst[check_inst == flag].index
        dataframe_inst = dataframe_inst[dataframe_inst.index.isin(record_to_keep)]
        print(f"{dataframe_inst.shape[0]} records after {inst} instrument screening")
        return dataframe_inst

    def filter_df_by_remove_duplicates(self, df: pd.DataFrame,
                                       subset: list = None) -> pd.DataFrame:
        """
        Remove duplicate records from the dataframe based on a given subset of columns such as [data_source, compound identifier, rt].

        :param df: DataFrame to process.
        :param subset: List of columns to consider for identifying duplicates.
        :return: DataFrame with duplicates removed.
        """
        if subset is None:
            subset = ["CH@rt", "ID@db_source", "C@inchikey"]

        df = df.drop_duplicates(subset=subset)
        print(
            f"{df.shape[0]} records after removing duplicates by {subset}, have {len(df['C@inchikey'].unique())} unique inchikey")
        return df


def print_collate_info(data: pd.DataFrame) -> pd.DataFrame:
    # number of institutes
    MBK_institute = data["ID@db_source"].unique()
    # number of compound records for each institutes
    num_records = data["ID@db_source"].value_counts()
    # number of unique compound records for each institutes
    new_MBK_filtered_unique_record = data.drop_duplicates(
        subset=["ID@db_source", "C@inchikey"])
    num_unique_records = new_MBK_filtered_unique_record["ID@db_source"].value_counts()
    print_info = pd.DataFrame({"num_compounds": num_records, "num_unique_compounds": num_unique_records})
    print(f"small datasets contain {len(MBK_institute)} systems,{data.shape[0]} analytes")
    return print_info
