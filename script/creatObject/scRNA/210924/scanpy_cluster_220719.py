# -*- coding: utf-8 -*-
"""
Created on Tue Jul 19 14:21:40 2022

@author: PC
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

os.chdir("E://Cryo-TCR/server/210924")
import utils

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=200, facecolor='white', figsize=(8,8))

TNK = sc.read_h5ad("T_sub_real.h5ad")
TNK = utils.get_raw_adata(TNK)

utils.qc(TNK,0,0)
sc.pl.violin(TNK,
             ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], 
             groupby='orig.ident',
             jitter=0, 
             multi_panel=True)

TNK = utils.normalize_variable_scale(TNK, 2500, True, True)

utils.reduce(TNK)

TNK

sc.pl.umap(TNK, color = ['IFN_Final', 'IFN_Final_Cluster'], legend_loc='on data', legend_fontsize=8)
