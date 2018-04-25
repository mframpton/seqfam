import os
from seqfam.relatedness import Relatedness
import pandas as pd


#Create the relatedness object.
data_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data"))
wf_file = os.path.join(data_dir,"relatedness",".".join(["king","kinship","ibs"]))
cohort_fam = os.path.join(data_dir,"cohort.fam")
relatedness = Relatedness(wf_file=wf_file,cohort_fam=cohort_fam,bf_file=None)

#Within-family duplicates
wf_duplicate_l = relatedness.find_duplicates(bf_b=False)
print(wf_duplicate_l)
#Between-family duplicates... (Uncomment if you have a bf_file).'''
#bf_duplicate_l = relatedness.find_duplicates(bf_b=True)
#print(bf_duplicate_l)

#Expected versus observed within-family relationships.'''
exp_obs_df = relatedness.get_exp_obs_df()
print(exp_obs_df.loc[(exp_obs_df["EXP_REL"]!=exp_obs_df["OBS_REL"]) & (pd.notnull(exp_obs_df["Kinship"])),:])
