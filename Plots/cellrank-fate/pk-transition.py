#!/usr/bin/env python
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
import gc
import sys
import argparse
import pegasus as pg

################
# configure
################
sc.settings.autoshow = False
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=150, dpi_save=300, format='png', frameon=False, transparent=True, fontsize=10)
# sc.set_figure_params(dpi=200, dpi_save=200)
plt.rcParams["image.aspect"] = "equal"
plt.rcParams["figure.figsize"] = ([3,3])
colorrs = ["#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#8CD17D","#B6992D","#499894","#E15759","#FF9D9A","#79706E",
           "#D37295","#FABFD2","#B07AA1","#D4A6C8","#9D7660","#E58606", "#5D69B1", "#24796C",
           '#DAA51B', '#000000', '#99C945', '#ED645A']
colorrSS = ["#4E79A7","#A0CBE8","#8CD17D","#499894","#F28E2B","#FFBE7D","#B6992D","#E15759","#FF9D9A","#79706E",
           "#D37295","#FABFD2", '#B07AA1',"#B07AA1","#D4A6C8","#9D7660","#E58606", "#5D69B1", "#24796C","#499894",
           '#DAA51B', '#000000', '#99C945', '#ED645A']

warnings.simplefilter(action='ignore', category=FutureWarning)

scv.settings.verbosity = 3
scv.settings.set_figure_params(dpi=150, dpi_save=300, format='png', frameon=False, transparent=True, fontsize=8)
sc.settings.set_figure_params(dpi=150, dpi_save=300, format='png', frameon=False, transparent=True, fontsize=8)
cr.settings.verbosity = 2
scv.settings.set_figure_params("scvelo")
os.chdir("./RNA/HSPC/cellrank/HSC_GMP")
parser = argparse.ArgumentParser(description='compute pseudotime and transition matrix')
parser.add_argument('--input','-i',type=str,required=True, help='input anndata file')
parser.add_argument('--groups','-g',type=bool,required=True, help='whether grouping or not')
# parser.add_argument('--count','-c',type=str,required=True, help='anndata path (with count data)')
# parser.add_argument('--anno','-a',type=str,required=True, help='anndata path (with ccell type annotation)')
# parser.add_argument('--ct_loc','-cl',type=str,required=True, help='cell type column)')
# parser.add_argument('--celltype','-ct',type=str,required=True, help='cell type')
# parser.add_argument('--ids','-ids',type=str,required=False, help='ids')
# parser.add_argument('--output','-o',type=str,required=True, help='output loom file path')
args = parser.parse_args()

# --------------------------------------------
# All samples; All celltype paga + fa + dpt 
# --------------------------------------------
adata = sc.read(args.input)
print("successfully load data......")
gs = ['HC', 'ITP', 'SLE-ITP.MKP_low', 'SLE-ITP.MKP_high', 'SLE-NP']

# if args.groups == "True":
bdata = adata.copy()
for i in gs:
    adata = bdata[bdata.obs['group']==i].copy()
    sc.tl.paga(adata, groups='cell_type_update')
    sc.pl.paga(adata, threshold=0.1, show=False,fontsize=5)
    print("start {group} draw graph......".format(group=i))
    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata, color=['cell_type_update'])
    adata.layers["spliced"] = adata.X
    adata.layers["unspliced"] = adata.X
    scv.pp.filter_and_normalize(adata,min_shared_counts=20, n_top_genes=2000,subset_highly_variable=False) #,
    scv.pp.moments(adata)
    print("start {group} diffmap......".format(group=i))
    sc.tl.diffmap(adata, n_comps=100)
    adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type_update']  == 'HSC')[0]
    sc.tl.dpt(adata, n_dcs=100)
    pk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
    print("start {group} compute_transition_matrix......".format(group=i))
    pk.compute_transition_matrix()
    pk.write_to_adata()
    adata.write("hspc_"+i+"_fa.pseudotime.h5ad")
else:
    # for looping each group if 
    sc.tl.paga(adata, groups='cell_type_update')
    sc.pl.paga(adata, threshold=0.1, show=False,fontsize=5)
    print("start draw graph......")
    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata, color=['cell_type_update'])
    adata.layers["spliced"] = adata.X
    adata.layers["unspliced"] = adata.X
    scv.pp.filter_and_normalize(adata,min_shared_counts=20, n_top_genes=2000,subset_highly_variable=False) #,
    scv.pp.moments(adata)
    print("start diffmap......")
    sc.tl.diffmap(adata, n_comps=100)
    adata.uns['iroot'] = np.flatnonzero(adata.obs['cell_type_update']  == 'HSC')[0]
    sc.tl.dpt(adata, n_dcs=100)
    # sc.tl.diffmap(adata, n_comps=100)
    # sc.tl.dpt(adata, n_dcs=100)
    pk = PseudotimeKernel(adata, time_key="dpt_pseudotime")
    print("start compute_transition_matrix......")
    pk.compute_transition_matrix()
    pk.write_to_adata()
    adata.write("hspc_all_fa.pseudotime.h5ad")
        
del adata
del pk
gc.collect()

