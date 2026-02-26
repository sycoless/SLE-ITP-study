import scanpy as sc
# import anndata
from anndata import AnnData
# from scipy import io
# from scipy.sparse import coo_matrix, csr_matrix
import scvelo as scv
import cellrank as cr
import matplotlib.pyplot as plt
import cellrank as cr
import numpy as np
import pandas as pd
import logging, os, sys
import warnings
import copy
from cellrank.estimators import GPCCA
from cellrank.kernels import PseudotimeKernel
from cellrank.kernels import ConnectivityKernel
import scanpy.external as sce
from cellrank.estimators import GPCCA

def compute_prob_HSC(adatas,Num):
    adata = sc.read("hspc_"+adatas+"_fa.pseudotime.h5ad")
    pk = cr.kernels.PseudotimeKernel.from_adata(adata, key="T_fwd")
    ck = ConnectivityKernel(adata)
    g = GPCCA(pk)
    g.fit(n_states=Num, cluster_key="cell_type_update")
    g.plot_macrostates(which="all",save="./hspc_"+adatas+"_fa.pseudotime.macrostates.pdf")
    return adata,g,pk
    
################
# configure 
################
sc.settings.autoshow = False
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, dpi_save=300, format='png', frameon=False, transparent=True, fontsize=10)
plt.rcParams["image.aspect"] = "equal"
plt.rcParams["figure.figsize"] = ([3,3])
colorrs = ["#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#8CD17D","#B6992D","#499894","#E15759","#FF9D9A","#79706E",
           "#D37295","#FABFD2","#B07AA1","#D4A6C8","#9D7660","#E58606", "#5D69B1", "#24796C",
           '#DAA51B', '#000000', '#99C945', '#ED645A']
colorrSS = ["#4E79A7","#A0CBE8","#8CD17D","#499894","#F28E2B","#FFBE7D","#B6992D","#E15759","#FF9D9A","#79706E",
           "#D37295","#FABFD2", '#B07AA1',"#B07AA1","#D4A6C8","#9D7660","#E58606", "#5D69B1", "#24796C","#499894",
           '#DAA51B', '#000000', '#99C945', '#ED645A']

warnings.simplefilter(action='ignore', category=FutureWarning)


os.chdir("./RNA/HSPC/cellrank/HSC_GMP/")

cts = ['HSC', 'MPP', 'MPP_myeloid']
term_stats = ['ProB', 'MKP', 'EryP', 'GMP', 'EBMP']
lineage_columns = ['ProB', 'MKP', 'EryP', 'GMP', 'EBMP']

# List of configurations: (condition, n_clusters, disease_filter, palette, save_prefix)
configs = [
    ("ITP", 7, "ITP", "#9E6D7F", "ITP"),
    ("HC", 5, "HC", "#B2CFE5", "HC"),
    ("SLE-NP", 5, "SLE-NP", "#D97351", "NP"),
    ("SLE-ITP.MKP_high", 5, None, "#8A2A46", "SLE-ITP.MKP_high"),
    ("SLE-ITP.MKP_low", 7, None, "#e5b5b5", "SLE-ITP.MKP_low")
]

def process_condition(cond, n, disease_filter, palette, save_prefix):
    # Compute probabilities
    adata, g, pk = compute_prob_HSC(cond, n)
    g.set_terminal_states(states=term_stats)
    g.compute_fate_probabilities()
    pk.write_to_adata()
    adata.write(f"subset_hspc_{save_prefix}_fa.pseudotime-pk.h5ad")
    
    # Filter adata
    if disease_filter:
        subset = adata[(adata.obs['cell_type_update'].isin(cts)) & (adata.obs['Disease'] == disease_filter)]
    else:
        subset = adata[(adata.obs['cell_type_update'].isin(cts))]
    
    # Plot circular projection
    cr.pl.circular_projection(
        subset, keys=["Disease"], legend_loc="right", lineages=term_stats,
        lineage_order="default", palette=[palette],
        title="", label_distance=1.2, legend_fontsize=8, figsize=(10, 10),
        dpi=400, size=15, save=f'subset.Circular_{save_prefix}_trajectories-group{"22" if "MKP" in save_prefix else ""}.pdf'
    )
    
    # Extract and save transition probabilities
    trans_probs = pd.DataFrame(adata.obsm['lineages_fwd'], columns=lineage_columns)
    trans_probs['Max'] = trans_probs.idxmax(axis=1)
    trans_probs['Cell ID'] = adata.obs.index
    trans_probs.index = adata.obs.index
    trans_probs['Sample'] = adata.obs['Sample'].astype('str')
    trans_probs['Disease'] = adata.obs['Disease'].astype('str')
    trans_probs['cell_type'] = adata.obs['cell_type'].astype('str')
    trans_probs['cell_type_update'] = adata.obs['cell_type_update'].astype('str')
    trans_probs.to_csv(f'subset_{save_prefix}_to_terminal_states.csv')

# Run for all configurations
for config in configs:
    process_condition(*config)


# ====================================
# Proportion
# ====================================
gs = ['HC', 'ITP', 'SLE-ITP.MKP_low', 'SLE-ITP.MKP_high', 'SLE-NP']
dfs = []
for i in gs:
    df = pd.read_csv("subset_"+i+"_to_terminal_states.csv")
    df['source'] = i
    dfs.append(df)
final_df = pd.concat(dfs, ignore_index=True)
print(final_df)

import seaborn as sns
import matplotlib.pyplot as plt
final_df = final_df[final_df['cell_type_update']!="ProB"].copy()
res = final_df.groupby('source').apply(lambda x: (x['Max'] == 'MKP').sum() / len(x))
res = res.reset_index(name='MKP_ratio')

custom_order = ['HC',  'ITP', 'SLE-NP', 'SLE-ITP.MKP_high','SLE-ITP.MKP_low']
res['source_cat'] = pd.Categorical(res['source'], categories=custom_order, ordered=True)
print(res)
res  = res.sort_values(by='source_cat')

ax = res.plot.bar(y='MKP_ratio', figsize=(3,5), grid=False, width=0.8,
                        edgecolor='black', color=["#B2CFE5","#9E6D7F","#D97351","#8A2A46","#e5b5b5"], legend=False
                         )

ax.bar_label(ax.containers[0], fmt='%.2f', fontsize=10, fontweight='bold')

ax.set_ylim(0,0.05)
ax.set_xticklabels(['HC','ITP','SLE-NP','SLE-ITP.MKP high','SLE-ITP.MKP low'], rotation=90, fontsize=12)
ax.set_xlabel('')
ax.set_ylabel('Prop. of HSCs differentiating to MK', fontsize=12)

sns.despine()
plt.tight_layout()
 plt.savefig('figures/subset_Proportion_MK_differentiation.pdf')