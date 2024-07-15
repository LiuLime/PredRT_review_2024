"""
Fetching compound class through 「ClassyFire Batch by Fiehn Lab(http://cfb.fiehnlab.ucdavis.edu/)」
Return: json file
Usage: run ./exucte_classyfire.sh in terminal
Parameter:
--load_ckp
--start_point
--data_path
--query_column
--output_name
--load_ckp_file

@Liu 2023/12/13
"""

import sys

sys.path.append("../collate_dataset")
from utils import sendURL
import pandas as pd
import json
import argparse
import numpy as np


def get_parser():
    with open("../config.json", "r") as c:
        config = json.load(c)

    parser = argparse.ArgumentParser(description="download molecular classfication from classyfire")

    parser.add_argument(
        "--load_ckp",
        default=False,
        help="restart from checkpoint",
        action="store_true",
    )
    parser.add_argument(
        "--start_point",
        default=config["classyfire@start_point"],
        help="restart from checkpoint",
        type=int
    )
    parser.add_argument(
        "--data_path",
        default=config["classyfire@path"],
        help="file contains the list of query compound,default=None",
    )
    parser.add_argument(
        "--query_column",
        default=config["classyfire@query_column"],
        help="query content in defined column,default=C@inchikey",
        type=str
    )
    parser.add_argument(
        "--query_namespace",
        default=config["classyfire@query_namespace"],
        help="query classyfire via inchikey from http://cfb.fiehnlab.ucdavis.edu/",
        type=str
    )
    parser.add_argument(
        "--output_name",
        default=config["classyfire@output_name"],
        help="the name of output query json file,default=classyfire_query.json",
        type=str
    )
    parser.add_argument(
        "--load_ckp_file",
        default=config["classyfire@load_ckp_file"],
        help="the name of load query json file,default=classyfire_query.json",
        type=str
    )

    return parser.parse_args()


def read_data_file(args):
    try:
        with open(args.data_path, "r") as S:
            dataset = pd.read_csv(S, header=0, sep=",", )

    except pd.errors.ParserError:
        with open(args.data_path, "r") as S:
            dataset = pd.read_csv(S, header=0, sep=";", )
    return dataset


def get_url(query) -> dict | None |str:
    headers = {
        'authority': 'cfb.fiehnlab.ucdavis.edu',
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
    url = f"http://cfb.fiehnlab.ucdavis.edu/entities/{query}.json"
    ans = sendURL.send_url(url, headers=headers, verify=False)
    if ans:
        data = json.loads(ans)
        kingdom = data.get("kingdom")
        superclass = data.get("superclass")
        _class = data.get("class")
        subclass = data.get("subclass")
        entry = {"inchikey": data.get("_id"),
                 "kingdom": kingdom.get("name", None) if kingdom else None,
                 "superclass": superclass.get("name", None) if superclass else None,
                 "class": _class.get("name", None) if _class else None,
                 "subclass": subclass.get("name", None) if subclass else None,
                 "state": "found"}
        print(f"{query} url success")
    else:
        entry = {"state": "not found"}
        print(f"{query} url failed")
    return entry

def _save_json(classyfire_info, name):
    with open(name, "w") as j:
        json.dump(classyfire_info, j)


def _load_json(name):
    with open(name, "r") as f:
        return json.load(f)


def save_checkpoint_if_needed(classyfire_info, count, args):
    if (count + 1) % 50 == 0:
        print(
            f"save check point at count {count + 1}, {len(classyfire_info)} compounds, reload start_point should be {count + 2}",
            flush=True)
        _save_json(classyfire_info, args.output_name)


def load_checkpoint(args):
    classyfire_info = _load_json(args.load_ckp_file)
    print(f"load {args.load_ckp_file} checkpoint, dataset have {len(classyfire_info)} records")
    return classyfire_info


def main(args):
    start_point = args.start_point
    classyfire_info = load_checkpoint(args) if args.load_ckp else {}
    dataset = read_data_file(args)
    dataset = dataset.replace(["N/A", "NA", "na", "n/a", np.nan], None)
    query_list = dataset[args.query_column].dropna().unique().tolist()


    for i, query_item in enumerate(query_list[start_point:], start=start_point):
        if query_item in classyfire_info:
            continue
        query_result = get_url(query_item)
        if query_result:  # if query is None because of other errors, will not save in json file
            classyfire_info[query_item] = query_result
            save_checkpoint_if_needed(classyfire_info, i, args)

    print(f"dataset save as {args.output_name},have {len(classyfire_info)} records finally")
    _save_json(classyfire_info, args.output_name)


if __name__ == "__main__":
    args = get_parser()
    print("process start")
    main(args)
    print("process finish")
