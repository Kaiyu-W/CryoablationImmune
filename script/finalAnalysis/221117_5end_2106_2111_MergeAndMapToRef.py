#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 16:00:38 2022

@author: kaiyu
"""

import os
import sys
from os.path import exists
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import anndata as ad

sys.path.append("/home/kaiyu/Desktop/Github/My_scRNA_pipeline")
from utils import *

sc.settings.set_figure_params(dpi=200, facecolor='white')

os.chdir("/home/kaiyu/Desktop/Cryo")
T_2106 = sc.read_h5ad("T_2106.h5ad")
T_2111 = sc.read_h5ad("T_2111.h5ad")
T_3end = sc.read_h5ad("T_3end.h5ad")

T_2111.obs_names = [ str(x).replace('wk', 'wk2111') for x in T_2111.obs_names ]
T_2106.obs_names = [ str(x).replace('wk', 'wk2106') for x in T_2106.obs_names ]


sc.pl.umap(T_3end, color='IFN_Final', legend_loc='on data', size=10, legend_fontsize=8)
sc.pl.tsne(T_2106, color = 'FinalCluster_2', legend_loc='on data', legend_fontsize=8, size=10)

ax = sc.pl.umap(T_2111, color='MainCluster_unique_detail_fromSub', legend_loc='on data', size=10, 
           legend_fontsize=6, frameon=False, show=False)
ax.title.set_size(10)
plt.show()

ax = sc.pl.umap(T_3end, color='IFN_Final', legend_loc='on data', size=10, 
           legend_fontsize=6, frameon=False, show=False)
ax.title.set_size(10)
plt.show()

# Integrate
    # # 1.read data
    # adata_ref = T_2111.copy()
    # adata = T_2106.copy()
    # common_vars = [ x for x in adata.var.index if x in adata_ref.var.index ]
    # adata_ref = adata_ref[:, common_vars]
    # adata = adata[:, common_vars]
    
    # # 2.integrate
    # sc.tl.ingest(adata, adata_ref, obs='MainCluster_unique_detail_fromSub')
    # adata.uns['MainCluster_unique_detail_fromSub_colors'] = adata_ref.uns['MainCluster_unique_detail_fromSub_colors']  # fix colors
    # sc.pl.umap(adata, color=['MainCluster_unique_detail_fromSub', 'FinalCluster_2'], legend_loc='on data', size=10, legend_fontsize=8)
    
    # # 3.combine
    # adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
    # adata_concat.obs.MainCluster_unique_detail_fromSub = adata_concat.obs.MainCluster_unique_detail_fromSub.astype('category')
    # adata_concat.obs.MainCluster_unique_detail_fromSub.cat.reorder_categories(adata_ref.obs.MainCluster_unique_detail_fromSub.cat.categories, inplace=True)  # fix category ordering
    # adata_concat.uns['MainCluster_unique_detail_fromSub_colors'] = adata_ref.uns['MainCluster_unique_detail_fromSub_colors']  # fix category colors
    # sc.pl.umap(adata_concat, color=['batch', 'MainCluster_unique_detail_fromSub'])


# 1.read data
adata_ref = T_3end.raw.to_adata()
adata_ref = adata_ref[
    [ x not in ['NK', 'NK_Proliferative'] for x in adata_ref.obs.IFN_Final]
]
adata1 = T_2111.raw.to_adata()
adata2 = T_2106.raw.to_adata()
adata = ad.concat([adata1,adata2])

qc(adata,0,0)
qc(adata_ref,0,0)
adata_ref = normalize_variable_scale(adata_ref)
adata = normalize_variable_scale(adata)
sc.tl.pca(adata, svd_solver='arpack', random_state=123)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50, random_state=123)
sc.tl.pca(adata_ref, svd_solver='arpack', random_state=123)
del adata_ref.uns['neighbors']['params']
sc.pp.neighbors(adata_ref, n_neighbors=20, n_pcs=50, random_state=123)
adata_ref.obsm['X_umap_0'] = adata_ref.obsm['X_umap']
sc.tl.umap(adata_ref, random_state=123)
sc.pl.umap(adata_ref, color='IFN_Final', legend_loc='on data', size=10, legend_fontsize=8)
adata_ref.obsm['X_umap'] = adata_ref.obsm['X_umap_0']

common_vars = [ x for x in adata.var.index if x in adata_ref.var.index ]
adata_ref = adata_ref[:, common_vars]
adata = adata[:, common_vars]

# 2.integrate
sc.tl.ingest(adata, adata_ref, obs='IFN_Final')
adata.uns['IFN_Final_colors'] = adata_ref.uns['IFN_Final_colors']  # fix colors
sc.pl.umap(adata, color=['IFN_Final'], legend_loc='on data', size=10, legend_fontsize=8)

pd.DataFrame(adata.obsm["X_umap"]).to_csv("T5_2106_2111_X_umap_3end.csv")
adata.obs.to_csv("T5_2106_2111_obs.csv")

# 3.combine
adata_concat = adata_ref.concatenate(adata, batch_categories=['ref', 'new'])
adata_concat.obs.IFN_Final = adata_concat.obs.IFN_Final.astype('category')
adata_concat.obs.IFN_Final.cat.reorder_categories(adata_ref.obs.IFN_Final.cat.categories, inplace=True)  # fix category ordering
adata_concat.uns['IFN_Final_colors'] = adata_ref.uns['IFN_Final_colors']  # fix category colors
sc.pl.umap(adata_concat, color=['batch', 'IFN_Final'])


T_markers_dict = {
    # NK
    'NK' : ["Klrb1c", "Ncr1", "Klrk1", "Nkg7"],
    # T
    'T' : ['Cd3e'],
    'CD4' : ["Cd4"],
    'CD8' : ['Cd8a', 'Cd8b1'],
    'Naive' : ['Sell'],
    # CD8
    'Activated' : ["Cd69", "Il2ra"],
    'Cytotoxic' : ['Gzma', 'Gzmb'],
    'Effector' : ["Cd44"],
    'Memory' : ['Il7r'],
    'Exhausted' : ["Pdcd1", "Ctla4", "Lag3"],
    # CD4
    'Th1' : ['Cxcr3'],
    # 'Th2' : ['Gata3'],
    'Th17' : ['Il17a'],
    'Tfh' : ['Il21r', 'Sh2d1a'],
    'Treg' : ['Foxp3']
}
T_markers = [y for x in T_markers_dict.values() for y in x]
ct = list(T_markers_dict.keys())

sc.pl.dotplot(T_2111, T_markers_dict, groupby='IFN_Final', mean_only_expressed=True, log=True)
# sc.pl.matrixplot(T_2111, T_markers_dict, groupby='IFN_Final', log=True)
# sc.pl.stacked_violin(T_2111, T_markers_dict, groupby='IFN_Final', rotation=45)
# sc.pl.heatmap(T_2111, T_markers_dict, groupby='IFN_Final', log=True, swap_axes=True)

pd.DataFrame(T_2111.obsm["X_umap_3end"]).to_csv("T_2111_X_umap_3end.csv")
T_2111.obs.to_csv("T_2111_obs.csv")

T_2111_obs = pd.read_csv("T_2111_meta_data_pathwayscore.csv",index_col=0)
T_2111.obs = T_2111_obs

pathway_names = [
    "GOBP_RESPONSE_TO_TYPE_I_INTERFERON","GOBP_RESPONSE_TO_INTERFERON_ALPHA",   
    "GOBP_RESPONSE_TO_INTERFERON_BETA","GOBP_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA",
    "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_BETA"
]
sc.pl.embedding(T_2111, 'X_umap_3end', color = pathway_names + ['IFN_Final'], 
                ncols = 2, size=10, wspace=0.4, groups="CD8_Tem_I-IFN", 
                frameon=False, na_in_legend=False)
