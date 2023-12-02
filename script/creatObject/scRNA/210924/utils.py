from os.path import exists
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import anndata
import scanpy as sc
import seaborn as sns
from typing import Union, Optional
import warnings

def cbind(df1, df2):
    return pd.concat([df1, df2], axis=1)
def rbind(df1, df2):
    return pd.concat([df1, df2], axis=0)

# def get_raw_adata(adata, uns_key = 'raw_adata'):
#     assert uns_key in adata.uns.keys()
#     raw = adata.uns[uns_key]
#     if type(raw) is anndata._core.raw.Raw:
#         raw = raw.to_adata()
#     adata_raw = adata.raw.to_adata()
#     adata_raw.X = raw.X
#     return adata_raw
def get_raw_adata(adata, **unusekwargs):
    assert adata.raw is not None
    adata_raw = adata.raw.to_adata()
    if 'log1p' in adata.uns.keys():
        base = adata.uns['log1p']['base']
        if base is not None:
            adata_raw.X = round(adata_raw.X.expm1().power(np.log(base)))
        else:
            adata_raw.X = round(adata_raw.X.expm1())
    
    return adata_raw

def update_category(
    df: pd.DataFrame, 
    key_add: str, 
    col: str, 
    new_category: pd.Series
):
    if key_add in df.columns:
        print(f"df['{key_add}'] has existed and will be overwritten!")
    df[key_add] = df[col].astype(str)
    df.loc[new_category.index.tolist(), key_add] = new_category
    df[key_add] = df[key_add].astype('category')

def rename_category_not_unique(
    df: pd.DataFrame, 
    col: str, 
    new_category: Union[tuple, list, dict]
):
    assert col in df.columns, 'df did not have this column!'
    categoried_data = df[col]
    if not pd.api.types.is_categorical_dtype(categoried_data):
        categoried_data = categoried_data.astype('category')

    old_category = categoried_data.cat.categories
    assert len(old_category) == len(new_category), \
        "Length of the new category should be equal to the old one!"

    data = categoried_data.astype(old_category.dtype)
    if isinstance(new_category, (list, tuple)):
        data = [ new_category[list(old_category).index(x)] for x in data ]
    elif isinstance(new_category, dict):
        data = [ new_category[x] for x in data ]
    else:
        print("new category should be one of tuple, list, dict.")
        exit(1)
    
    df[col] = data
    df[col] = df[col].astype('category')

def qc(adata, min_genes=300, min_cells=3):
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata.var['mt'] = adata.var_names.str.contains('[Mm][Tt]-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    
def filter(adata, nGenes_max = 6000, nCount_max = 30000, pctMt_max = 15, nGenes_min = 300, nCount_min = 0):
    adata = adata[adata.obs.n_genes_by_counts < nGenes_max, :]
    adata = adata[adata.obs.total_counts < nCount_max, :]
    adata = adata[adata.obs.pct_counts_mt < pctMt_max, :]
    adata = adata[adata.obs.n_genes_by_counts > nGenes_min, :]
    adata = adata[adata.obs.total_counts > nCount_min, :]
    return adata
    
def normalize_variable_scale(adata, ntop = 2500, plot = False, seurat_vst = False):
    # adata.raw = adata
    # adata.uns['raw_adata'] = adata.raw
    if seurat_vst:
        assert 'vst.variable' in adata.var.keys()
        assert 'vst.variance.standardized' in adata.var.keys()
        adata.var['highly_variable'] = adata.var['vst.variable'] > 0
        adata.var['highly_variable_rank'] = adata.var['vst.variance.standardized'][
            adata.var['highly_variable']
            ].rank(
                method = 'first', 
                ascending=False
            )
        if plot:
            warnings.warn("Set plot=False when seurat_vst=True")
    else:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=ntop)
        if plot:
            sc.pl.highly_variable_genes(adata)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata, max_value=10)
    return adata

def reduce(adata, n_neighbors = 20, n_pcs = 50, tsne = True, SEED = 123):
    sc.tl.pca(adata, svd_solver='arpack', random_state=SEED)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=SEED)
    sc.tl.umap(adata, random_state=SEED)
    if tsne:
        sc.tl.tsne(adata, random_state=SEED)
        
def cluster(adata, resolution = 0.5, method = 'leiden', SEED = 123):
    if method == 'leiden':
        fun = sc.tl.leiden
    elif method == 'louvain':
        fun = sc.tl.leiden
    else:
        return
    key = f"{method}.{resolution}"
    fun(adata, resolution = resolution, key_added=key, random_state=SEED)
    print(f"added adata.obs['{key}']")
    
def FindMarkers(adata, groupby, method = 'wilcoxon', plot = False, nPlot = 25):
    assert method in ['t-test', 'wilcoxon']
    sc.tl.rank_genes_groups(adata, groupby=groupby, method=method)
    if plot: sc.pl.rank_genes_groups(adata, n_genes=nPlot, sharey=False)

def count_BiCluster(adata, group1, group2, percent = False):
    assert group1 in adata.obs.keys()
    assert group2 in adata.obs.keys()
    count_dict = dict(adata.obs[[group1, group2]].value_counts())
    count_df = pd.DataFrame(
        0, 
        index=np.unique([x for (x,y) in count_dict.keys()]),
        columns=np.unique([y for (x,y) in count_dict.keys()])
    )
    for key in count_dict.keys():
        (x,y) = key
        count_df.loc[x,y]=count_dict[key]
    if percent:
        count_df = count_df.apply(lambda x: 100 * x / x.sum())
    else:
        count_df['sum'] = count_df.apply(np.sum, axis=1)
    count_df.loc['sum'] = count_df.apply(np.sum)
    
    return count_df

def integrate(adata, method = 'bbknn', batch_key = 'batch', SEED = 123, **kwargs):
    assert method in ['bbknn', 'harmony', 'mnn', 'scanorama']
    adata_int = adata.copy()
    if method == 'bbknn':
        sc.external.pp.bbknn(adata_int, batch_key=batch_key, **kwargs)
        sc.tl.umap(adata_int, random_state=SEED)
    elif method == 'harmony':
        sc.external.pp.harmony_integrate(adata_int, key=batch_key, **kwargs)
        adata_int.obsm['X_pca'],adata_int.obsm['X_pca_raw'] = adata_int.obsm['X_pca_harmony'],adata_int.obsm['X_pca']
        sc.pp.neighbors(adata_int, **adata_int.uns['neighbors']['params'])
        sc.tl.umap(adata_int, random_state=SEED)
    elif method == 'mnn':
        X_raw = adata.X
        adata_int = sc.external.pp.mnn_correct(adata_int, batch_key=batch_key, **kwargs)[0][0]
        X = adata_int.X
        if (X_raw == X).all():
            print("Same data before and after mnn-integration")
            return adata
        sc.tl.pca(adata_int)
        sc.pp.neighbors(adata_int, **adata_int.uns['neighbors']['params'])
        sc.tl.umap(adata_int, random_state=SEED)
    elif method == 'scanorama':
        sc.external.pp.scanorama_integrate(adata_int, key=batch_key, **kwargs)
        adata_int.obsm['X_pca'],adata_int.obsm['X_pca_raw'] = adata_int.obsm['X_scanorama'],adata_int.obsm['X_pca']
        sc.pp.neighbors(adata_int, **adata_int.uns['neighbors']['params'])
        sc.tl.umap(adata_int, random_state=SEED)
    return adata_int