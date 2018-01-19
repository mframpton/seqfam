import os
from cohpy.relatedness import Relatedness
import pandas as pd

'''Create the relatedness object.'''
wf_file = os.path.abspath(os.path.join("..","..","data","relatedness",".".join(["king","kinship","ibs"])))
cohort_tsv = os.path.abspath(os.path.join("..","..","data",".".join(["cohort","tsv"])))
relatedness = Relatedness(bf_file=None,wf_file,cohort_tsv)

'''Within-family duplicates'''
wf_duplicate_l = relatedness.find_duplicates(bf_b=False)
print(wf_duplicate_l)
'''Between-family duplicates... (Uncomment if you have a bf_file).'''
#bf_duplicate_l = relatedness.find_duplicates(bf_b=True)
#print(bf_duplicate_l)

'''Expected versus observed within-family relationships.'''
exp_obs_df = relatedness.get_exp_obs_df()
print(exp_obs_df.ix[(exp_obs_df["EXP_REL"]!=exp_obs_df["OBS_REL"]) & (pd.notnull(exp_obs_df["Kinship"])),:])
