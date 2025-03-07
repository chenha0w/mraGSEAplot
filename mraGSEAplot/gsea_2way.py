from ._regulon_convert import regulon_p2r
import pandas as pd
import rpy2.robjects as robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
import gseapy as gp

base = rpackages.importr('base')

def MR_gsea_2way(MR,sig,regulon, n_thread=4, maxsize=2500, seed=2):

    # sig: signature. It can be pandas Series, dictionary or pd.DataFrame
    # MR: name of MR; string
    # regulon: network; It can be pd.DataFrame or ListVector from R

    #type check and conversion
    if isinstance(sig, dict):
        sig=pd.Series(sig)
    elif isinstance(sig, pd.DataFrame):
        sig=sig.iloc[:,0]
    if not isinstance(sig, pd.Series):
        raise TypeError('The input signature need to be dictionary or pd.Series or pd.DataFrame with 1 column')
    if isinstance(regulon, pd.DataFrame):
        regulon = regulon_p2r(regulon)
    if not isinstance(regulon, robjects.vectors.ListVector):
        raise TypeError('The input regulon has to be either rpy2 listvectors or dataframe')

    # gsea on two sets
    rnk = sig.sort_values(ascending=False)
    gene_sets = {'induced': list(base.names(regulon.rx2(MR).rx2('tfmode')).rx(regulon.rx2(MR).rx2('tfmode').ro>0)),
             'repressed': list(base.names(regulon.rx2(MR).rx2('tfmode')).rx(regulon.rx2(MR).rx2('tfmode').ro<0))}
    pre_res = gp.prerank(rnk=rnk, gene_sets=gene_sets, thread=n_thread, max_size=maxsize, seed=seed) 
    return pre_res.results