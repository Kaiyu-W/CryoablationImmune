#!/usr/bin/env python
# coding: utf-8

# # Configure library

# In[1]:


import os
import sys
from os.path import exists
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


# In[2]:


sys.path.append("/home/kaiyu/Desktop/Github/My_scRNA_pipeline/")
from utils import *
sys.path.append("/home/kaiyu/conga/conga-master")
import conga


# In[3]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=200, facecolor='white', figsize=(6,6))


# # Load data

# In[4]:


os.chdir("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/")

T5 = sc.read_h5ad("T_5end.h5ad").raw.to_adata()
# T5.var.index = T5.var['_index']
T5.var.rename(columns = {'_index':'index'}, inplace = True)
T5.var = T5.var.set_index('index')

T_real_cb = pd.read_csv("T5_cellbarcode.csv", index_col=0)
T_real_umap = pd.read_csv("T5_umap.csv", index_col=0)
T_real_meta = pd.read_csv("T5_meta.csv", index_col=0)

T5


# In[5]:


T = T5[T_real_cb.x, :]
T.obsm['X_umap'] = T_real_umap.to_numpy()
T.obs = T_real_meta
T


# In[6]:


colors_dict = {
    "CD4_Naive" : "#709ae1",
    "CD4_T_I-IFN" : "#197ec0",
    "Th1" : "#d2af81",
    "Th17" : "#bd559f",
    "Treg" : "#f7c530",
    "Tfh" : "#46732e",
    "CD8_Naive" : "#d5e4a2",
    "CD8_Tem_I-IFN" : "#ff410d",
    "CD8_Tem" : "#95cc5e",
    "CD8_Tex" : "#d0dfe6",
    "CD8_Tex_Proliferative" : "#f79d1e",
    "T_others" : "#748aa6"
}


# In[7]:


from pandas.api.types import CategoricalDtype
cat = CategoricalDtype(categories=colors_dict.keys())
T.obs['IFN_Final'] = T.obs['IFN_Final'].astype(cat)
T.uns['IFN_Final_colors'] = list(colors_dict.values())


# In[8]:


sc.pl.umap(T, color = ['IFN_Final'], legend_loc='on data', legend_fontsize=8, size = 10)


# In[9]:


T.obs['orig.ident'].value_counts()


# # Load TCR

# In[10]:


# conga data input
os.chdir("/home/kaiyu/Desktop/Cryo/Paper_Result_Plot/TCR_files")

gex_datafile = 'T_5end.h5ad'
gex_datatype = 'h5ad' # other possibilities right now: ['10x_mtx', 'h5ad'] (h5ad from scanpy)
organism = 'mouse'

tcr_AllTCR = 'TCR_T5_all.csv'

clones_tcr_AllTCR =  'TCR_clones_All.tsv'
outfile_prefix = 'TCR_clones_All'

# save GEX
if not exists(gex_datafile):
    T.write_h5ad(gex_datafile)

# save TCR
TCR_files = [
    'CryoTCR_1wk.txt',
    'CryoTCR_2wk.txt',
    'NonCryoTCR_1wk.txt',
    'NonCryoTCR_2wk.txt',
    'CryoTCR_1wk_2.txt',
    'NonCryoTCR_1wk_2.txt',
]
Prefix = [
    'Cryo1wk2111',
    'Cryo2wk2111',
    'NonCryo1wk2111',
    'NonCryo2wk2111',
    'Cryo1wk2106',
    'NonCryo1wk2106',
]
if not exists(tcr_AllTCR):
    data_all = pd.DataFrame()
    for file,prefix in zip(TCR_files,Prefix):
        data = pd.read_csv(file)
#         data.barcode = [ f"{prefix}_{x}" for x in data.barcode ]
#         data.contig_id = [ f"{prefix}_{x}" for x in data.contig_id ]
#         data.raw_clonotype_id = [ f"{x}_{prefix}" for x in data.raw_clonotype_id ]
#         data.raw_consensus_id = [ str(y).replace("consensus", f"{prefix}_consensus") for y in data.raw_consensus_id ]
        print(file,prefix,data.shape)
        data_all = pd.concat([data_all, data], axis=0)
    
    data_all.to_csv(tcr_AllTCR, index=0)
    print(data_all.head())
    
assert exists(tcr_AllTCR)
assert exists(gex_datafile)


# In[11]:


# this creates the TCRdist 'clones file'
if not exists(clones_tcr_AllTCR):
    conga.tcrdist.make_10x_clones_file.make_10x_clones_file( tcr_AllTCR, organism, clones_tcr_AllTCR )

# this command will create another file with the kernel PCs for subsequent reading by conga
if not exists(clones_tcr_AllTCR.split(".tsv")[0] + "_AB.dist_50_kpcs"):
    conga.preprocess.make_tcrdist_kernel_pcs_file_from_clones_file( clones_tcr_AllTCR, organism )


# In[12]:


adata = conga.preprocess.read_dataset( gex_datafile, gex_datatype, clones_tcr_AllTCR )
# os.remove(gex_datafile)

# store the organism info in adata
adata.uns['organism'] = organism

adata


# In[13]:


sc.pl.umap(adata, color = ['IFN_Final', 'orig.ident'], legend_loc='on data', legend_fontsize=8, size = 10)


# In[14]:


adata_raw = adata.copy()


# In[15]:


# adata.obs['Case_Control'] = adata.obs['orig.ident']
# rename_category_not_unique(adata.obs, 'Case_Control', [0,1,2,3])
# adata.obs['Case_Control'] = adata.obs['Case_Control'].cat.as_ordered()
# adata.uns['batch_keys'] = ['Case_Control']
# sc.pl.umap(adata, color = ['MainCluster_byGrouped_Merge','Case_Control'], legend_loc='on data', legend_fontsize=8, size = 10)


# In[15]:


adata_raw.obs['orig.ident'].value_counts()


# # Conga Process

# In[16]:


# get raw gex

# idx = [ x in adata.obs.index for x in T.obs.index ]
# X = adata.uns['raw_adata'][idx, :].X
# var = adata.raw.to_adata().var
# uns = adata.uns
# obs = adata.obs
# obsm = adata.obsm
# obsp = adata.obsp
# adata = anndata.AnnData(X=X, obs=obs,obsm=obsm,obsp=obsp,uns=uns, var=var)
# adata

# adata = adata_raw.raw.to_adata()
adata = adata_raw.copy()

adata.obsm['X_umap_gex'] = adata.obsm['X_umap']
adata.obsm['X_tsne_gex'] = adata.obsm['X_tsne']
adata.obsm['X_pca_gex'] = adata.obsm['X_pca']


# In[17]:


adata.obs['Case_Control'] = adata.obs['orig.ident'].astype("category")
adata.obs['Case_Control'].cat.categories = [0,1,2,3]
adata.obs['Case_Control'] = adata.obs['Case_Control'].cat.as_ordered()
# adata.uns['batch_keys'] = ['Case_Control']
# import re
# adata.obs['Group'] = [ 0 if re.match("^Non.*", x) else 1 for x in adata.obs['orig.ident'] ]
# adata.obs['Time'] = [ 0 if re.match(".*1wk$", x) else 1 for x in adata.obs['orig.ident'] ]
# adata.uns['batch_keys'] = ['Group','Time']


# In[18]:


qc(adata,0,0)

# # do not filtering since that did in R::Seurat

adata = conga.preprocess.filter_and_scale( 
    adata, 
    min_genes_per_cell=0,
    max_genes_per_cell=500000,
    max_percent_mito=1,
#     hvg_batch_key = 'case_control'
)
# adata

# adata.uns['raw_matrix_is_logged'] = True
adata = conga.preprocess.reduce_to_single_cell_per_clone(adata, use_existing_pca_obsm_tag = 'X_pca_gex')
adata


# In[19]:


from pathlib import Path
import random
from scipy.sparse import issparse

def calc_tcrdist_nbrs_cpp(
    adata,
    num_nbrs: int,
    tmpfile_prefix: Optional[str] = None, # old name outfile_prefix, confusing
    n_components_umap: int = 2,
    verbose: bool = True
):
    '''
    stolen from conga.preproces.calc_tcrdist_nbrs_umap_clusters_cpp
    remove the cluster function, only keep that of nbrs computation

    stores in:
    adata.uns['tcr_neighbors']
    adata.obsp['tcr_distances']
    adata.obsp['tcr_connectivities']
    '''
    assert 'conga.util' in sys.modules.keys(), 'No package conga! '+         'Please install it by pipeline in https://github.com/phbradley/conga'
    for obs_tcr in 'va cdr3a vb cdr3b'.split():
        assert obs_tcr in adata.obs.keys()

    if tmpfile_prefix is None:
        tmpfile_prefix = f'tmp_tcrdists{random.randrange(1,10000)}'

    # save tcr info to file for tcrdist_cpp executable
    tcrs_filename = tmpfile_prefix+'_tcrs.tsv'
    adata.obs['va cdr3a vb cdr3b'.split()].to_csv(
        tcrs_filename, sep='\t', index=False)

    if os.name == 'posix':
        exe = Path.joinpath(Path(conga.util.path_to_tcrdist_cpp_bin),
                            'find_neighbors')
    else:
        exe = Path.joinpath(Path(conga.util.path_to_tcrdist_cpp_bin),
                            'find_neighbors.exe')

    if not exists(exe):
        print('need to compile c++ exe:', exe)
        exit(1)

    db_filename = Path.joinpath(
        Path(conga.util.path_to_tcrdist_cpp_db),
        'tcrdist_info_{}.txt'.format(adata.uns['organism']))

    if not exists(db_filename):
        print('need to create database file:', db_filename)
        exit(1)

    outprefix = tmpfile_prefix+'_calc_tcrdist'

    while True:
        # umap runs into trouble if there are too many connected components
        #  in the neighbor graph
        # So we increase num_nbrs inside here if components>2*n_components_umap
        #  based on looking at umap/spectral.py and trying to avoid the call
        #  to component_layout
        #
        # compute the tcrdist neighbors
        cmd = '{} -f {} -n {} -d {} -o {}'.format(
            exe, tcrs_filename, num_nbrs, db_filename, outprefix)

        conga.util.run_command(cmd, verbose=verbose)

        knn_indices_filename = outprefix+'_knn_indices.txt'
        knn_distances_filename = outprefix+'_knn_distances.txt'

        if (not exists(knn_indices_filename) or
            not exists(knn_distances_filename)):
            print('find_neighbors failed:', exists(knn_indices_filename),
                  exists(knn_distances_filename))
            exit(1)

        knn_indices = np.loadtxt(knn_indices_filename, dtype=int)
        knn_distances = np.loadtxt(knn_distances_filename, dtype=float)

        #distances=sc.neighbors.get_sparse_matrix_from_indices_distances_numpy(
        #     knn_indices, knn_distances, adata.shape[0], num_nbrs)

        try: # HACK: the naming of this function changes across scanpy versions
            distances,connectivities= sc.neighbors.compute_connectivities_umap(
                knn_indices, knn_distances, adata.shape[0], num_nbrs)
        except:
            if verbose:
                print('try new name for compute_connectivities_umap')
            distances,connectivities= sc.neighbors._compute_connectivities_umap(
                knn_indices, knn_distances, adata.shape[0], num_nbrs)

        if issparse(connectivities): # I think this is always true
            from scipy.sparse.csgraph import connected_components
            connected_components = connected_components(connectivities)
            number_connected_components = connected_components[0]
            if verbose:
                print('number_connected_components:', number_connected_components)
            if number_connected_components > 2*n_components_umap:
                if verbose:
                    print('tcrdist umap, too many connected components in the',
                          'neighbor graph:', number_connected_components,
                          '--> increasing num_nbrs from', num_nbrs,'to',
                          2*num_nbrs)
                num_nbrs *= 2
                continue
        break # break out of the infinite loop

    ################
    # stash the stuff in adata, stolen from scanpy/neighbors/__init__.py
    #
    adata.uns['tcr_neighbors'] = {}
    adata.uns['tcr_neighbors']['params']={'n_neighbors': num_nbrs, 'method': 'umap'}
    adata.uns['tcr_neighbors']['params']['metric'] = 'tcrdist'# fake metric
    adata.uns['tcr_neighbors']['connectivities_key'] = 'tcr_connectivities'
    adata.uns['tcr_neighbors']['distances_key'] = 'tcr_distances'
    adata.obsp['tcr_distances'] = distances
    adata.obsp['tcr_connectivities'] = connectivities
    print("tcr_network stored in adata.obsp['tcr_connectivities']")

    # cleanup the tmpfiles
    for filename in [knn_indices_filename,
                     knn_distances_filename,
                     tcrs_filename]:
        os.remove(filename)


# In[20]:


# adata = conga.preprocess.cluster_and_tsne_and_umap( 
#     adata,
#     clustering_resolution = 0.5,
#     recompute_pca_gex = False,
#     skip_tsne = True,
# #     skip_tsne = False,
#     clustering_method = 'leiden',
#     n_neighbors = 20,
#     n_gex_pcs = 50
# )

calc_tcrdist_nbrs_cpp(adata, 20)
sc.tl.umap(adata, neighbors_key='tcr_neighbors', n_components=2)
adata.obsm["X_umap_tcr"] = adata.obsm["X_umap"]

sc.tl.leiden(
    adata, resolution=0.4, random_state=0, 
    key_added='clusters_tcr', adjacency=adata.obsp['tcr_connectivities']
)

adata.obs['clusters_tcr'] = adata.obs['clusters_tcr'].astype(int)
conga.preprocess.setup_tcr_cluster_names(adata) #stores in adata.uns
adata.obs['clusters_tcr'] = adata.obs['clusters_tcr'].astype('category')


# In[21]:


sc.pl.embedding(adata, basis='X_umap_tcr', color = ['clusters_tcr','orig.ident'], legend_loc='on data', size=15, legend_fontsize=10)


# In[22]:


# adata.obs['clusters_gex0'] = adata.obs.clusters_gex.astype('category')
adata.obs['clusters_gex'] = adata.obs.IFN_Final.astype(str)

clusters_gex_name = adata.obs['IFN_Final'].cat.categories
adata.uns['clusters_gex_names'] = clusters_gex_name #stores in adata.uns

df_tmp = adata.obs['clusters_gex'].copy()
for i in range(len(clusters_gex_name)):
    cluster_tmp = clusters_gex_name[i]
    df_tmp[df_tmp == cluster_tmp] = i
adata.obs['clusters_gex'] = df_tmp.astype('category')


# In[23]:


# these are the nbrhood sizes, as a fraction of the entire dataset:
nbr_fracs = [0.01, 0.1]

# we use this nbrhood size for computing the nndists
nbr_frac_for_nndists = 0.01

all_nbrs, nndists_gex, nndists_tcr = conga.preprocess.calc_nbrs(
    adata, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)

# stash these in obs array, they are used in a few places...
adata.obs['nndists_gex'] = nndists_gex
adata.obs['nndists_tcr'] = nndists_tcr


# In[24]:


adata.obs['TCR_clusters'] = adata.obs.clusters_tcr.astype('category')
adata.rename_categories('TCR_clusters', adata.uns['clusters_tcr_names'])

sc.pl.embedding(
    adata, basis='X_umap_tcr', color = 'TCR_clusters', 
    legend_loc='on data', size=15, legend_fontsize=10
)


# In[25]:


adata.obs['TCR_clusters'].cat.categories


# In[26]:


adata.uns['TCR_clusters_colors']


# In[27]:


results = conga.correlations.run_graph_vs_graph(
    adata, all_nbrs, outfile_prefix=outfile_prefix)

results.head()


# In[28]:


#put the conga hits on top
adata.obs['log_conga_scores'] = np.maximum(-1*np.log10(adata.obs.conga_scores),0.0)

sc.pl.embedding(
    adata, basis='X_umap_gex', color = 'log_conga_scores', 
    legend_loc='on data', size=20, legend_fontsize=10, vmax=5
)


# In[29]:


# adata.uns['log1p']['base'] = 2
adata.obsm['X_gex_2d'] = adata.obsm['X_umap_gex'] 
adata.obsm['X_tcr_2d'] = adata.obsm['X_umap_tcr']


# In[30]:


# # # plot by UMAP

# # default
# nbrs_gex, nbrs_tcr = all_nbrs[0.1]

# min_cluster_size = 2

# conga.plotting.make_graph_vs_graph_logos(
#     adata,
#     outfile_prefix + '_umap',
#     min_cluster_size,
#     nbrs_gex,
#     nbrs_tcr,
#     include_full_tcr_cluster_names_in_logo_lines = True
# #     clusters_gex_names = adata.uns['clusters_gex_names'],
# #     clusters_tcr_names = adata.uns['clusters_tcr_names']
# )

# # adata.uns['clusters_tcr_names']

# # # or equivalently:
# # #
# # # conga.plotting.make_graph_vs_graph_logos(
# # #     adata,
# # #     outfile_prefix,
# # #     min_cluster_size,
# # #     nbrs_gex,
# # #     nbrs_tcr,
# # #     gex_header_genes = ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
# # #     gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
# # #     logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL','GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1','CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
# # #     gene_logo_width = 6,
# # # )

# # # the batch barplot: (from author)
# # # Right now the batch colors are arranged in increasing order from bottom to top, 
# # # 0 = blue at the bottom, then 1 above that, etc. 
# # # The 'tab10' colormap is used if there are fewer than 11 batches, 
# # # otherwise the 'tab20' colormap is used: https://matplotlib.org/stable/gallery/color/colormap_reference.html
# # # Each batch_key (e.g. 'patient' or 'diagnostic_group' ) gets a bar that is divided left/right with the left, 
# # # thicker part showing the batch composition of the corresponding cluster and the right, 
# # # thinner part showing the batch composition of the full dataset (so you can see enrichment/depletion).


# In[46]:


# set ImageMagick path

os.environ['PATH'] = os.environ['PATH'] + ":" + "/home/kaiyu/anaconda3/envs/cassiopeia/bin/"


# ## CD8 

# In[47]:


if not exists('TCR_clones_CD8'):
    os.makedirs('TCR_clones_CD8')
outfile_prefix_CD8 = 'TCR_clones_CD8/TCR_clones_CD8'

CD8_idx = adata.obs.IFN_Final.isin([
    'CD8_Naive','CD8_Tem_I-IFN', 'CD8_Tem', 'CD8_Tex', 'CD8_Tex_Proliferative'
])

# # plot by UMAP

# default
min_cluster_size = 2

CD8 = adata[CD8_idx].copy()
all_nbrs_CD8, nndists_gex_CD8, nndists_tcr_CD8 = conga.preprocess.calc_nbrs(
    CD8, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)
nbrs_gex_CD8, nbrs_tcr_CD8 = all_nbrs_CD8[0.1]


conga.plotting.make_graph_vs_graph_logos(
    CD8,
    outfile_prefix_CD8 + '_umap',
    min_cluster_size,
    nbrs_gex_CD8,
    nbrs_tcr_CD8,
    include_full_tcr_cluster_names_in_logo_lines = True,
    nocleanup = True
#     clusters_gex_names = adata.uns['clusters_gex_names'],
#     clusters_tcr_names = adata.uns['clusters_tcr_names']
)

# adata.uns['clusters_tcr_names']

# # or equivalently:
# #
# # conga.plotting.make_graph_vs_graph_logos(
# #     adata,
# #     outfile_prefix,
# #     min_cluster_size,
# #     nbrs_gex,
# #     nbrs_tcr,
# #     gex_header_genes = ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
# #     gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
# #     logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL','GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1','CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
# #     gene_logo_width = 6,
# # )

# # the batch barplot: (from author)
# # Right now the batch colors are arranged in increasing order from bottom to top, 
# # 0 = blue at the bottom, then 1 above that, etc. 
# # The 'tab10' colormap is used if there are fewer than 11 batches, 
# # otherwise the 'tab20' colormap is used: https://matplotlib.org/stable/gallery/color/colormap_reference.html
# # Each batch_key (e.g. 'patient' or 'diagnostic_group' ) gets a bar that is divided left/right with the left, 
# # thicker part showing the batch composition of the corresponding cluster and the right, 
# # thinner part showing the batch composition of the full dataset (so you can see enrichment/depletion).


# ## CD4

# In[48]:


if not exists('TCR_clones_CD4'):
    os.makedirs('TCR_clones_CD4')
outfile_prefix_CD4 = 'TCR_clones_CD4/TCR_clones_CD4'

CD4_idx = adata.obs.IFN_Final.isin([
    'CD4_Naive', 'CD4_T_I-IFN', 'Th1', 'Th17', 'Treg', 'Tfh'
])

# # plot by UMAP

# default
min_cluster_size = 2

CD4 = adata[CD4_idx].copy()
all_nbrs_CD4, nndists_gex_CD4, nndists_tcr_CD4 = conga.preprocess.calc_nbrs(
    CD4, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)
nbrs_gex_CD4, nbrs_tcr_CD4 = all_nbrs_CD4[0.1]

conga.plotting.make_graph_vs_graph_logos(
    CD4,
    outfile_prefix_CD4 + '_umap',
    min_cluster_size,
    nbrs_gex_CD4,
    nbrs_tcr_CD4,
    include_full_tcr_cluster_names_in_logo_lines = True,
    nocleanup = True
#     clusters_gex_names = adata.uns['clusters_gex_names'],
#     clusters_tcr_names = adata.uns['clusters_tcr_names']
)

# adata.uns['clusters_tcr_names']

# # or equivalently:
# #
# # conga.plotting.make_graph_vs_graph_logos(
# #     adata,
# #     outfile_prefix,
# #     min_cluster_size,
# #     nbrs_gex,
# #     nbrs_tcr,
# #     gex_header_genes = ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
# #     gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
# #     logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL','GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1','CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
# #     gene_logo_width = 6,
# # )

# # the batch barplot: (from author)
# # Right now the batch colors are arranged in increasing order from bottom to top, 
# # 0 = blue at the bottom, then 1 above that, etc. 
# # The 'tab10' colormap is used if there are fewer than 11 batches, 
# # otherwise the 'tab20' colormap is used: https://matplotlib.org/stable/gallery/color/colormap_reference.html
# # Each batch_key (e.g. 'patient' or 'diagnostic_group' ) gets a bar that is divided left/right with the left, 
# # thicker part showing the batch composition of the corresponding cluster and the right, 
# # thinner part showing the batch composition of the full dataset (so you can see enrichment/depletion).


# ## Stat

# In[49]:


clusters_gex = adata.obs['IFN_Final']
clusters_tcr = adata.obs['TCR_clusters']
good_score_mask = adata.obs.conga_scores < 1
node2cluspair = {
        i:(x,y) for i,(x,y,m) in enumerate(zip(clusters_gex,
                                               clusters_tcr,
                                               good_score_mask)) if m}
node2cluspair_raw = {
        i:(x,y) for i,(x,y) in enumerate(zip(clusters_gex,
                                             clusters_tcr))}
cluspair2nodes = {}
for node, clp in node2cluspair.items():
    cluspair2nodes.setdefault(clp,[]).append(node)
cluspair2nodes_raw = {}
for node, clp in node2cluspair_raw.items():
    cluspair2nodes_raw.setdefault(clp,[]).append(node)
    
bicluster2group = {}
for key,value in cluspair2nodes.items():
    value2 = adata[value,:].obs.groupby('orig.ident').agg({"clone_sizes": "sum"}).to_dict()['clone_sizes']
#     print(adata[value,:].obs.groupby('orig.ident').agg({"clone_sizes": "sum"}))
    bicluster2group[key] = value2
bicluster2group_raw = {}
for key,value in cluspair2nodes_raw.items():
    value2 = adata[value,:].obs.groupby('orig.ident').agg({"clone_sizes": "sum"}).to_dict()['clone_sizes']
#     print(f"{key}, {value2}")
    bicluster2group_raw[key] = value2


# In[50]:


num_clusters_gex = len(np.unique(clusters_gex))+1
num_clusters_tcr = len(np.unique(clusters_tcr))+1
print(f"num_clusters_gex={num_clusters_gex},num_clusters_tcr={num_clusters_tcr}")


# In[51]:


bicluster_count = pd.DataFrame.from_dict(bicluster2group, orient = 'index')
bicluster_count = bicluster_count.applymap(lambda x: 0 if np.isnan(x) else x)
bicluster_count = bicluster_count[['NonCryo2wk','NonCryo1wk','Cryo2wk','Cryo1wk']]
bicluster_count_scale = bicluster_count.apply(lambda x: x/np.sum(x), axis=1)
bicluster_count = bicluster_count.reset_index()
bicluster_count = bicluster_count.sort_values(by = 'level_0', ignore_index=T)
bicluster_count['level_0'] = bicluster_count['level_0'].astype('category')
bicluster_count['level_1'] = bicluster_count['level_1'].astype('category')
bicluster_count_scale = bicluster_count_scale.reset_index()
bicluster_count_scale = bicluster_count_scale.sort_values(by = 'level_0', ignore_index=T)
bicluster_count_scale['level_0'] = bicluster_count_scale['level_0'].astype('category')
bicluster_count_scale['level_1'] = bicluster_count_scale['level_1'].astype('category')


# In[52]:


bicluster_count_raw = pd.DataFrame.from_dict(bicluster2group_raw, orient = 'index')
bicluster_count_raw = bicluster_count_raw.applymap(lambda x: 0 if np.isnan(x) else x)
bicluster_count_raw = bicluster_count_raw[['NonCryo2wk','NonCryo1wk','Cryo2wk','Cryo1wk']]
bicluster_count_scale_raw = bicluster_count_raw.apply(lambda x: x/np.sum(x), axis=1)
bicluster_count_raw = bicluster_count_raw.reset_index()
bicluster_count_raw = bicluster_count_raw.sort_values(by = 'level_0', ignore_index=T)
bicluster_count_raw['level_0'] = bicluster_count_raw['level_0'].astype('category')
bicluster_count_raw['level_1'] = bicluster_count_raw['level_1'].astype('category')
bicluster_count_scale_raw = bicluster_count_scale_raw.reset_index()
bicluster_count_scale_raw = bicluster_count_scale_raw.sort_values(by = 'level_0', ignore_index=T)
bicluster_count_scale_raw['level_0'] = bicluster_count_scale_raw['level_0'].astype('category')
bicluster_count_scale_raw['level_1'] = bicluster_count_scale_raw['level_1'].astype('category')


# In[53]:


bicluster_count


# In[54]:


bicluster_count[bicluster_count.level_0 == 'CD4_Naive'].apply(
    lambda x: [x[1], x[2:].sum()], 
    axis = 1
)


# In[55]:


bicluster_count[bicluster_count.level_0 == 'CD8_Naive'].apply(
    lambda x: [x[1], x[2:].sum()], 
    axis = 1
)


# In[56]:


bicluster_count_raw


# In[57]:


bicluster_count_scale


# In[58]:


bicluster_count_scale_raw


# In[44]:


# save for R to plot

bicluster_count.to_csv("bicluster_count_df.csv")
bicluster_count_scale.to_csv("bicluster_count_scale_df.csv")
bicluster_count_raw.to_csv("bicluster_count_raw_df.csv")
bicluster_count_scale_raw.to_csv("bicluster_count_scale_raw_df.csv")
adata.obs[['TCR_clusters','conga_scores','clone_sizes']].to_csv("bicluster_obs_df.csv")


# In[ ]:





# ## CD8 Naive Tem

# In[59]:


CD8_Naive_top3TCRcluster = ['3_AV12','4_bv13','7_BV1']


# In[60]:


if not exists('TCR_clones_CD8naive'):
    os.makedirs('TCR_clones_CD8naive')
outfile_prefix_CD8 = 'TCR_clones_CD8naive/TCR_clones_CD8'

CD8_idx = adata.obs.IFN_Final.isin([
    'CD8_Naive','CD8_Tem'
]) & adata.obs.TCR_clusters.isin(CD8_Naive_top3TCRcluster)

# # plot by UMAP

# default
min_cluster_size = 25

CD8 = adata[CD8_idx].copy()
all_nbrs_CD8, nndists_gex_CD8, nndists_tcr_CD8 = conga.preprocess.calc_nbrs(
    CD8, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)
nbrs_gex_CD8, nbrs_tcr_CD8 = all_nbrs_CD8[0.1]


conga.plotting.make_graph_vs_graph_logos(
    CD8,
    outfile_prefix_CD8 + '_umap',
    min_cluster_size,
    nbrs_gex_CD8,
    nbrs_tcr_CD8,
    include_full_tcr_cluster_names_in_logo_lines = True,
    nocleanup = True
#     clusters_gex_names = adata.uns['clusters_gex_names'],
#     clusters_tcr_names = adata.uns['clusters_tcr_names']
)

# adata.uns['clusters_tcr_names']

# # or equivalently:
# #
# # conga.plotting.make_graph_vs_graph_logos(
# #     adata,
# #     outfile_prefix,
# #     min_cluster_size,
# #     nbrs_gex,
# #     nbrs_tcr,
# #     gex_header_genes = ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
# #     gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
# #     logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL','GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1','CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
# #     gene_logo_width = 6,
# # )

# # the batch barplot: (from author)
# # Right now the batch colors are arranged in increasing order from bottom to top, 
# # 0 = blue at the bottom, then 1 above that, etc. 
# # The 'tab10' colormap is used if there are fewer than 11 batches, 
# # otherwise the 'tab20' colormap is used: https://matplotlib.org/stable/gallery/color/colormap_reference.html
# # Each batch_key (e.g. 'patient' or 'diagnostic_group' ) gets a bar that is divided left/right with the left, 
# # thicker part showing the batch composition of the corresponding cluster and the right, 
# # thinner part showing the batch composition of the full dataset (so you can see enrichment/depletion).


# In[ ]:





# ## CD4 Naive Th1

# In[61]:


CD4_Naive_top3TCRcluster = ['1_av13','4_bv13','7_BV1']


# In[62]:


if not exists('TCR_clones_CD4naive'):
    os.makedirs('TCR_clones_CD4naive')
outfile_prefix_CD4 = 'TCR_clones_CD4naive/TCR_clones_CD4'

CD4_idx = adata.obs.IFN_Final.isin([
    'CD4_Naive','Th1'
]) & adata.obs.TCR_clusters.isin(CD4_Naive_top3TCRcluster)

# # plot by UMAP

# default
min_cluster_size = 90

CD4 = adata[CD4_idx].copy()
all_nbrs_CD4, nndists_gex_CD4, nndists_tcr_CD4 = conga.preprocess.calc_nbrs(
    CD4, nbr_fracs, also_calc_nndists=True, nbr_frac_for_nndists=nbr_frac_for_nndists)
nbrs_gex_CD4, nbrs_tcr_CD4 = all_nbrs_CD4[0.1]

conga.plotting.make_graph_vs_graph_logos(
    CD4,
    outfile_prefix_CD4 + '_umap',
    min_cluster_size,
    nbrs_gex_CD4,
    nbrs_tcr_CD4,
    include_full_tcr_cluster_names_in_logo_lines = True,
    nocleanup = True
#     clusters_gex_names = adata.uns['clusters_gex_names'],
#     clusters_tcr_names = adata.uns['clusters_tcr_names']
)

# adata.uns['clusters_tcr_names']

# # or equivalently:
# #
# # conga.plotting.make_graph_vs_graph_logos(
# #     adata,
# #     outfile_prefix,
# #     min_cluster_size,
# #     nbrs_gex,
# #     nbrs_tcr,
# #     gex_header_genes = ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
# #     gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
# #     logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL','GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1','CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
# #     gene_logo_width = 6,
# # )

# # the batch barplot: (from author)
# # Right now the batch colors are arranged in increasing order from bottom to top, 
# # 0 = blue at the bottom, then 1 above that, etc. 
# # The 'tab10' colormap is used if there are fewer than 11 batches, 
# # otherwise the 'tab20' colormap is used: https://matplotlib.org/stable/gallery/color/colormap_reference.html
# # Each batch_key (e.g. 'patient' or 'diagnostic_group' ) gets a bar that is divided left/right with the left, 
# # thicker part showing the batch composition of the corresponding cluster and the right, 
# # thinner part showing the batch composition of the full dataset (so you can see enrichment/depletion).


# In[ ]:





# In[49]:


# from PIL import Image               # for looking at some of 
# from IPython.display import display # the output images.

# tag = conga.tags.GRAPH_VS_GRAPH_LOGOS
# pngfile = adata.uns['conga_results'].get(tag,None)
# if pngfile is None:
#     print('no graph-vs-graph clusters of sufficient size found')
# else:
#     help_message = adata.uns['conga_results'][tag+'_help']

#     print(help_message)
#     image = Image.open(pngfile)
# #     display(image)
    
#     r,g,b,a = image.split()
#     im = Image.merge('RGB', (r,g,b))
#     im.save('tmp.pdf','PDF')


# In[50]:


# # # plot by UMAP
# plot_batch = False
# if not plot_batch and 'batch_keys' in adata.uns.keys():
#     adata.uns['batch_keys_old'] = adata.uns['batch_keys']
#     del adata.uns['batch_keys']
# elif plot_batch and 'batch_keys_old' in adata.uns.keys():
#     adata.uns['batch_keys'] = adata.uns['batch_keys_old']
#     del adata.uns['batch_keys_old']

# # default
# nbrs_gex, nbrs_tcr = all_nbrs[0.1]

# min_cluster_size = 100

# conga.plotting.make_graph_vs_graph_logos(
#     adata,
#     outfile_prefix + '_umap2',
#     min_cluster_size,
#     nbrs_gex,
#     nbrs_tcr,
# #     make_batch_bars = True,
#     include_full_tcr_cluster_names_in_logo_lines = True
# #     clusters_gex_names = adata.uns['clusters_gex_names'],
# #     clusters_tcr_names = adata.uns['clusters_tcr_names']
# )

# # adata.uns['clusters_tcr_names']

# # # or equivalently:
# # #
# # # conga.plotting.make_graph_vs_graph_logos(
# # #     adata,
# # #     outfile_prefix,
# # #     min_cluster_size,
# # #     nbrs_gex,
# # #     nbrs_tcr,
# # #     gex_header_genes = ['clone_sizes','CD4','CD8A','CD8B','SELL','GNLY','GZMA','CCL5','ZNF683','IKZF2','PDCD1','KLRB1'],
# # #     gex_header_tcr_score_names = ['imhc', 'cdr3len', 'cd8', 'nndists_tcr'],
# # #     logo_genes = ['CD4','CD8A','CD8B','CCR7','SELL','GNLY','PRF1','GZMA','IL7R','IKZF2','KLRD1','CCL5','ZNF683','KLRB1','NKG7','HLA-DRB1' ],
# # #     gene_logo_width = 6,
# # # )

# # # the batch barplot: (from author)
# # # Right now the batch colors are arranged in increasing order from bottom to top, 
# # # 0 = blue at the bottom, then 1 above that, etc. 
# # # The 'tab10' colormap is used if there are fewer than 11 batches, 
# # # otherwise the 'tab20' colormap is used: https://matplotlib.org/stable/gallery/color/colormap_reference.html
# # # Each batch_key (e.g. 'patient' or 'diagnostic_group' ) gets a bar that is divided left/right with the left, 
# # # thicker part showing the batch composition of the corresponding cluster and the right, 
# # # thinner part showing the batch composition of the full dataset (so you can see enrichment/depletion).


# In[ ]:





# In[51]:


# conga.correlations.run_graph_vs_features(
#     adata, all_nbrs, outfile_prefix=outfile_prefix + '_umap2')


# In[52]:


# conga.plotting.make_graph_vs_features_plots(
#     adata, all_nbrs, outfile_prefix + '_umap2')


# In[ ]:





# In[53]:


# adata.obs['clusters_tcr'] = adata.obs['clusters_tcr'].cat.as_ordered()
# conga.correlations.find_hotspots_wrapper(
#     adata, all_nbrs, None, outfile_prefix + '_umap2')


# In[54]:


# conga.plotting.make_hotspot_plots(
#     adata, all_nbrs, outfile_prefix + '_umap2')


# # Others

# In[55]:


ncells_conga = np.sum(adata.obs.clone_sizes[adata.obs['conga_scores'] < 1])
nclones_conga = np.sum(adata.obs['conga_scores'] < 1)
print(f"conga scores < 1: n_clone={nclones_conga} n_cells={ncells_conga}")


# In[56]:


adata


# In[57]:


adata.uns


# In[58]:


# adata.uns['rank_genes_good_biclusters']['names']


# In[ ]:





# In[59]:


adata_good = adata[adata.obs.conga_scores < 1]
adata_good


# In[60]:


sc.pl.embedding(
    adata_good, 
    basis='X_umap_gex', 
    color = ['IFN_Final','TCR_clusters'], 
    legend_loc='on data', 
    size=15, 
    legend_fontsize=10
)


# In[61]:


adata_good.obs['bicluster'] = adata_good.obs.IFN_Final.astype(str) + "-" + adata_good.obs.TCR_clusters.astype(str)
# sc.pl.embedding(
#     adata_good, 
#     basis='X_umap_gex', 
#     color = ['bicluster'], 
#     legend_loc='on data', 
#     size=15, 
#     legend_fontsize=10
# )
adata_good.obs['bicluster'].value_counts()


# In[62]:


adata_good.obs['bicluster'].value_counts


# In[63]:


adata_good.obs.clone_sizes


# In[ ]:





# In[ ]:





# In[64]:


df_for_plot = bicluster_count_scale[['NonCryo2wk','NonCryo1wk','Cryo2wk','Cryo1wk','level_0']]

plt.figure(dpi=150)
plt.title(f'cell barcode coverage by UTR')

N = 103
ind = np.arange(N)

plt.xlabel('Clusters')
plt.ylabel('Percent')
plt.xticks(ind, df_for_plot['level_0'], rotation = 90)
plt.yticks(np.arange(0, 1.1, 0.1))

B0 = df_for_plot.iloc[:,0]
B1 = df_for_plot.iloc[:,1]
B2 = df_for_plot.iloc[:,2]
C1 = df_for_plot.iloc[:,3]

width = 1  # 设置条形图一个长条的宽度
p1 = plt.bar(ind, B0, width, color='lightcoral') 
p2 = plt.bar(ind, B1, width, bottom=B0,color='grey')  
p3 = plt.bar(ind, B2, width, bottom=B0+B1,color='green')
p4 = plt.bar(ind, C1, width, bottom=B0+B1+B2,color='blue')

plt.legend(
    (p4[0], p3[0], p2[0], p1[0]), 
    ('CA-1wk', 'CA-2wk', 'Non-CA-1wk', 'Non-CA-2wk'), 
    loc = 3, 
    bbox_to_anchor = (1.01, 0.5), borderaxespad = 0
)

plt.show()


# In[65]:


df_for_plot = bicluster_count[['NonCryo2wk','NonCryo1wk','Cryo2wk','Cryo1wk','level_0']]

plt.figure(dpi=150)
plt.title(f'cell barcode coverage by UTR')

N = 103
ind = np.arange(N)

plt.xlabel('Clusters')
plt.ylabel('Percent')
plt.xticks(ind, df_for_plot['level_0'], rotation = 90)
# plt.yticks(np.arange(0, 1, 0.1))

B0 = df_for_plot.iloc[:,0]
B1 = df_for_plot.iloc[:,1]
B2 = df_for_plot.iloc[:,2]
C1 = df_for_plot.iloc[:,3]

width = 1  # 设置条形图一个长条的宽度
p1 = plt.bar(ind, B0, width, color='lightcoral') 
p2 = plt.bar(ind, B1, width, bottom=B0,color='grey')  
p3 = plt.bar(ind, B2, width, bottom=B0+B1,color='green')
p4 = plt.bar(ind, C1, width, bottom=B0+B1+B2,color='blue')

plt.legend(
    (p4[0], p3[0], p2[0], p1[0]), 
    ('CA-1wk', 'CA-2wk', 'Non-CA-1wk', 'Non-CA-2wk'), 
    loc = 3, 
    bbox_to_anchor = (1.01, 0.5), borderaxespad = 0
)

plt.show()


# In[ ]:




