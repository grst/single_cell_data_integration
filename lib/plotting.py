import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from collections import OrderedDict

colors = dict()
colors["dataset"] = {
    "azizi_peer_2018": '#1f77b4',
    "azizi_peer_2018_10x": '#ff7f0e',
    "guo_zhang_2018": "#2ca02c",
    "lambrechts_2018_6149_v2": '#d62728',
    "lambrechts_2018_6653" : '#9467bd',
    "savas_loi_2018": '#8c564b',
    "zheng_bileas_2017": '#e377c2',
    "zheng_zhang_2017": '#7f7f7f'
}

colors["cell_type"] = {
 'B cell': '#1f77b4',
 'Mast cells': '#ff7f0e',
 'NK cell': '#2ca02c',
 'T cell CD4+': '#d62728',
 'T cell CD8+': '#9467bd',
 'T cell reg.': '#8c564b',
 'myeloblast-derived': '#e377c2',
 'pDC': '#17becf',
 'stem': '#bcbd22',
 'unknown': '#7f7f7f'
}

colors["clusters"] = OrderedDict([
 ('C0 - GZMK', '#e41a1c'),
 ('C1 - ZNF683', '#377eb8'),
 ('C2 - exhaustion', '#4daf4a'),
 ('C3 - ZNF683/chemokine', '#984ea3'),
 ('C4 - IL7R', '#ff7f00'),
 ('C5 - mitotic', '#ffff33'),
 ('C6 - heat shock', '#a65628'),
 ('C7 - IFIT', '#f781bf'),
 ('C8 - Immunoglobulin', '#999999')
])

colors["batch_patient"] = {'azizi_peer_2018_10x_BC09': '#FFFF00',
 'azizi_peer_2018_10x_BC10': '#1CE6FF',
 'azizi_peer_2018_10x_BC11': '#FF34FF',
 'azizi_peer_2018_BC1': '#FF4A46',
 'azizi_peer_2018_BC2': '#008941',
 'azizi_peer_2018_BC3': '#006FA6',
 'azizi_peer_2018_BC4': '#A30059',
 'azizi_peer_2018_BC5': '#FFDBE5',
 'azizi_peer_2018_BC6': '#7A4900',
 'azizi_peer_2018_BC7': '#0000A6',
 'azizi_peer_2018_BC8': '#63FFAC',
 'guo_zhang_2018_P0616A': '#B79762',
 'guo_zhang_2018_P0616P': '#004D43',
 'guo_zhang_2018_P0617': '#8FB0FF',
 'guo_zhang_2018_P0619': '#997D87',
 'guo_zhang_2018_P0706': '#5A0007',
 'guo_zhang_2018_P0729': '#809693',
 'guo_zhang_2018_P0913': '#FEFFE6',
 'guo_zhang_2018_P1010': '#1B4400',
 'guo_zhang_2018_P1011': '#4FC601',
 'guo_zhang_2018_P1118': '#3B5DFF',
 'guo_zhang_2018_P1120': '#4A3B53',
 'guo_zhang_2018_P1202': '#FF2F80',
 'guo_zhang_2018_P1208': '#61615A',
 'guo_zhang_2018_P1219': '#BA0900',
 'lambrechts_2018_6149_v2_3': '#6B7900',
 'lambrechts_2018_6149_v2_4': '#00C2A0',
 'lambrechts_2018_6149_v2_5': '#FFAA92',
 'lambrechts_2018_6653_6': '#FF90C9',
 'lambrechts_2018_6653_7': '#B903AA',
 'lambrechts_2018_6653_8': '#D16100',
 'savas_loi_2018_0': '#DDEFFF',
 'savas_loi_2018_1': '#000035',
 'zheng_bileas_2017_1': '#7B4F4B',
 'zheng_zhang_2017_P0205': '#A1C299',
 'zheng_zhang_2017_P0322': '#300018',
 'zheng_zhang_2017_P0407': '#0AA6D8',
 'zheng_zhang_2017_P0508': '#013349',
 'zheng_zhang_2017_P1116': '#00846F',
 'zheng_zhang_2017_P1202t': '#372101'}


def make_legend_elements(color_key):
    res = {
        "handles" : [],
        "labels": []
    }
    for label, color in colors[color_key].items():
        res["handles"].append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color, label=label, markersize=15))
        res["labels"].append(label)
    return res


def plot_umap(adata, ax, title="umap", size=0.01, rep="X_umap", color="dataset", **kwargs):
    if color in adata.raw.var_names:
        color_vec = adata.raw.X[:, adata.raw.var_names == color].toarray().flatten()
        color_vec = (color_vec - np.min(color_vec))/np.ptp(color_vec)
    else:
        assert np.all([x in colors[color] for x in np.unique(adata.obs[color].values)]), \
            "the color map contains a color for each distinct value"
        color_vec = np.array([colors[color][x] for x in adata.obs[color].values])
    ax.scatter(adata.obsm[rep][:, 0], adata.obsm[rep][:, 1], s=size,
            c=color_vec, marker='.', **kwargs)
    ax.set_xlabel("UMAP1")
    ax.set_xticks([])
    ax.set_ylabel("UMAP2")
    ax.set_yticks([])
    ax.set_title(title)