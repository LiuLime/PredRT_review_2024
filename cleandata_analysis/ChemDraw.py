"""Figure 1
Fig.1A: Ratio of compound superclass in SMRT, MassBank, MoNA, PredRet datasets (LC records only)
Fig.1B: Upset figure show intersections across datasets (LC records only)

output:
    - save_path: './04_result/'
    - Figure1A_pie.csv: statistical analysis result
    - Figure1A_pie.png: ratio of compound superclass across datasets
    - Figure1B_upset.png

@Liu 2023/12/21
"""
import os
import sys

sys.path.append("..")
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import from_contents, plot
import numpy as np
from _Filters import Filter


def read_file(dataset, save_path=None):
    with open(os.path.join(save_path, dataset), "r") as file:
        data = pd.read_csv(file, header=0, sep=",")
        data = data.replace(["N/A", "NA", "na", "n/a", "null", np.nan], None)
    return data


# read files
dataset_path = "04_result/"
save_path = "04_result/"
dataset = ["01_SMRT_LC.csv", "02_MBK_LC.csv", "03_MoNA_LC.csv", "04_predret_LC.csv"]
set_name = ["SMRT", "MassBank", "MoNA", "Predret"]
dataset = [read_file(dataset[i], dataset_path) for i in range(len(dataset))]

# remove records with blank inchikey
F = Filter(pd.DataFrame())
dataset_valid_inchikey = [F.keep_valid_molecular(i) for i in dataset]
# remove duplications according to inchikey
smrt, mbk, mona, predret = [i.drop_duplicates(subset=["C@inchikey"]) for i in dataset_valid_inchikey]

"""Figure 1A. chemical compound superclass pie graph"""

name = [smrt, mbk, mona, predret]
# read superclass information, fill NA by "Unknown"
smrt_superclass, mbk_superclass, mona_superclass, predret_superclass = \
    [pd.DataFrame(n["C@compound_superclass"].fillna("Unknown").value_counts()).reset_index() for n in name]

# statistical analyze
data = smrt_superclass \
    .merge(mbk_superclass, on="C@compound_superclass", how="outer", suffixes=("_smrt", "_mbk")) \
    .merge(mona_superclass, on="C@compound_superclass", how="outer", ) \
    .merge(predret_superclass, on="C@compound_superclass", how="outer", ) \
    .fillna(0) \
    .set_index("C@compound_superclass") \
    .set_axis(set_name, axis=1)
data["Total"] = data.sum(axis=1)
for dataset in set_name:
    data[dataset + '_pct'] = data[dataset] / data['Total']
data.sort_values(by="Total", ascending=False, inplace=True)
data.to_csv(os.path.join(save_path, "Figure1A_pie.csv"), sep=",")

# perpendicular pie figure 横向布局
fig, axes = plt.subplots(nrows=4, ncols=data.shape[0], figsize=(30, 20),
                         gridspec_kw={'height_ratios': [1, 1, 1, 1]}, )

for j, compound in enumerate(data.index):
    for i, dataset in enumerate(set_name):
        ax = axes[i, j]
        ax.pie([data.loc[compound, dataset + '_pct'], 1 - data.loc[compound, dataset + '_pct']],
               labels=None, startangle=90, counterclock=False,
               colors=['#66b3ff', '#fafafa'],
               wedgeprops={"edgecolor": "black", "linewidth": 1})
        ax.set(aspect="equal")
        ax.set_title(f"{int(data.loc[compound, dataset])}", size=30, pad=2, x=0.5, fontfamily="Arial")
        axes[i, j].set_visible(True)
        if i == 3:
            ax.set_xlabel(compound, size=30, labelpad=5, rotation=90, fontfamily="Arial")

plt.tight_layout()
plt.savefig(os.path.join(save_path, "Figure1A_pie.png"))
plt.show()

"""Figure 1B. Upset figure"""
smrt_inchikey = {"SMRT": smrt["C@inchikey"].unique().tolist()}
mbk_inchikey = {"MassBank": mbk["C@inchikey"].unique().tolist()}
mona_inchikey = {"MoNA": mona["C@inchikey"].unique().tolist()}
predret_inchikey = {"PredRet": predret["C@inchikey"].unique().tolist()}
new_dict = {**smrt_inchikey, **mbk_inchikey, **mona_inchikey, **predret_inchikey}

intersect = from_contents(new_dict)
plt.rcParams['font.family'] = "Arial"
upset = plot(intersect, subset_size="count", show_counts="{:d}")
plt.savefig(os.path.join(save_path, "Figure1B_upset.png"), dpi=300)
plt.show()
