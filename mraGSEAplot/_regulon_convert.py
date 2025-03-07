import rpy2.robjects.packages as rpackages
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
import pandas as pd

base = rpackages.importr('base')

def regulon_r2p(regulon):
    # Here regulon is a ListVector with length of number of MRs. Each MR has 2 Lists with tfmode and likelihood of its taragets.
    regulon_df = pd.DataFrame(columns=['MR','target','tfmode','likelihood'])
    for MR in base.names(regulon):
        temp_df = pandas2ri.rpy2py(base.as_data_frame(regulon.rx2[MR]))
        temp_df['MR']=MR
        temp_df = temp_df.reset_index().rename(columns={'index':'target'})
        regulon_df = pd.concat([regulon_df,temp_df],ignore_index=True)
    return regulon_df

def regulon_p2r(regulon_df):
    regulon_r=None
    for MR in regulon_df['MR'].unique():
        tfmode = robjects.FloatVector(regulon_df.loc[regulon_df['MR']==MR,'tfmode'])
        tfmode.names = robjects.StrVector(regulon_df.loc[regulon_df['MR']==MR,'target'])
        likelihood = robjects.FloatVector(regulon_df.loc[regulon_df['MR']==MR,'likelihood'])
        single_net = robjects.ListVector({'tfmode':tfmode,'likelihood':likelihood})
        if regulon_r is None:
            regulon_r = robjects.ListVector({MR:single_net})
        else:
            regulon_r.rx2[MR]=single_net
    return regulon_r