import pandas as pd
import sys
import numpy as np
import os
import itertools
import glob
import re
import relatedness as rf
import log as l
from natsort import natsorted
import argparse
import parent_args as pa


def get_s_p_gp(ind, fam_df):
    '''Get the siblings, parents and grandparents of an individual.'''
    f = fam_df.ix[fam_df["PERSON"]==ind,:].iloc[0]["FATHER"]
    m = fam_df.ix[fam_df["PERSON"]==ind,:].iloc[0]["MOTHER"]
    s_l = fam_df.ix[(fam_df["FATHER"]==f) & (fam_df["MOTHER"]==m) & (fam_df["PERSON"]!=ind),"PERSON"].tolist()
    p_l = filter(lambda x: pd.notnull(x), [f, m])
    [g1f, g1m, g2f, g2m] = [None, None, None, None]
    if pd.notnull(f):
        g1f = fam_df.ix[fam_df["PERSON"]==f,:].iloc[0]["FATHER"]
        g1m = fam_df.ix[fam_df["PERSON"]==f,:].iloc[0]["MOTHER"]
        g2f = fam_df.ix[fam_df["PERSON"]==m,:].iloc[0]["FATHER"]
        g2m = fam_df.ix[fam_df["PERSON"]==m,:].iloc[0]["MOTHER"]
    gp_l = filter(lambda x: pd.notnull(x), [g1f, g1m, g2f, g2m])

    return pd.Series([",".join(s_l),",".join(p_l),",".join(gp_l)], index=["S","P","GP"])


def get_exp_rel(ind_s, s_p_gp_df):
    '''Get the expected relationship between 2 individuals.'''
    ind1, ind2 = ind_s["ID1"], ind_s["ID2"]
    ind1_s_l, ind1_p_l, ind1_gp_l = filter(None,s_p_gp_df.loc[ind1,"S"].split(",")), filter(None,s_p_gp_df.loc[ind1,"P"].split(",")), filter(None,s_p_gp_df.loc[ind1,"GP"].split(","))
    ind2_s_l, ind2_p_l, ind2_gp_l = filter(None,s_p_gp_df.loc[ind2,"S"].split(",")), filter(None,s_p_gp_df.loc[ind2,"P"].split(",")), filter(None,s_p_gp_df.loc[ind2,"GP"].split(","))
    #print ind1, ind2
    #print ind1_s_l, ind1_p_l, ind1_gp_l
    #print ind2_s_l, ind2_p_l, ind2_gp_l
    exp_rel = "4" #i.e. at least a 4th-degree relationship
    if ind2 in ind1_p_l or ind1 in ind2_p_l: #Parent-child
        exp_rel = "1"
    elif len(ind1_p_l) > 0 and set(ind1_p_l) == set(ind2_p_l): #Sibling
        exp_rel = "1"
    elif ind2 in ind1_gp_l or ind1 in ind2_gp_l: #Grandparent-grandchild 
        exp_rel = "2"
    elif (len(ind1_p_l) > 1 and set(ind1_p_l).issubset(set(ind2_gp_l))) or (len(ind2_p_l) > 1 and set(ind2_p_l).issubset(set(ind1_gp_l))): #Uncle/Aunt - nephew/niece
        exp_rel = "2"
    elif len(set(ind1_gp_l).intersection(set(ind2_gp_l))) > 0: #Cousins
        exp_rel = "3"
    elif len(set(ind1_gp_l).intersection(set(ind2_s_l))) > 0 or len(set(ind2_gp_l).intersection(set(ind1_s_l))) > 0: #Great Uncle/Aunt
        exp_rel = "3"
    #print exp_rel
    
    return exp_rel


def get_obs_rel(s, obs_rel_df):
    '''Get the (KING-inferred) observed relationship between 2 individuals.'''
    obs_rel = np.NaN
    df1 = obs_rel_df.ix[(obs_rel_df["FIDID1"]=="_".join(s[["FID1","ID1"]].tolist())) & (obs_rel_df["FIDID2"]=="_".join(s[["FID2","ID2"]].tolist())),:]
    df2 = obs_rel_df.ix[(obs_rel_df["FIDID1"]=="_".join(s[["FID2","ID2"]].tolist())) & (obs_rel_df["FIDID2"]=="_".join(s[["FID1","ID1"]].tolist())),:]
    if df1.empty == False:
        obs_rel = df1.iloc[0]["OBS_REL"]
    elif df2.empty == False:
        obs_rel = df2.iloc[0]["OBS_REL"]
    return obs_rel


def get_errors(ind_eo_df,ind):
    #print ind
    '''Get errors'''
    ind_error_df = ind_eo_df.ix[ind_eo_df["OBS_REL"]!=ind_eo_df["EXP_REL"],:]
    error_str = ""
    if ind_error_df.empty == False:
        ind_error_df["error"] = ind_error_df.ix[:,["ID1","ID2","EXP_REL","OBS_REL"]].apply(axis=1, func=lambda x: "_".join(x))
        error_str = " ".join(ind_error_df["error"].tolist())
        error_str = error_str.replace(ind+"_","")
    return error_str


'''Set paths'''
parser = argparse.ArgumentParser(description="Compare the expected (pedigree) versus observed (KING-inferred) familial relationships.", parents=[pa.get_parent_parser()])
args = parser.parse_args()
pa.add_path_args(args)
log_stream = open(os.path.join(args.data_dir,".".join(["all","rel","familial","exp_v_obs","log","txt"])),'w')

'''Read in the observed (KING-inferred) relationships.'''
kinship_coef_df = rf.read_king_results(log_stream, args.data_dir)
#print kinship_coef_df[(kinship_coef_df["FID1"] == "108") | (kinship_coef_df["FID2"] == "108")].to_string()
#sys.exit()

'''Map the kinship coefficients to observed relations.'''
l.log(log_stream, "Creating OBS_REL column...")
kinship_coef_df["OBS_REL"] = kinship_coef_df["Kinship"].apply(lambda x: "1" if x >= 0.177 else "2" if x >= 0.0844 else "3" if x > 0.0442 else "4")

'''Get the unique individuals and exlude sporadics.'''
kinship_coef_df = rf.remove_sporadics(log_stream, kinship_coef_df)
l.log(log_stream, "# of individuals: " + str(len(pd.unique(kinship_coef_df["FIDID1"].tolist() + kinship_coef_df["FIDID2"].tolist()))))

'''Read in the expected relationships'''
results_df_l = []
fam_tsv_l = natsorted(glob.glob(os.path.join(args.pedigree_dir,"Family_*.tsv")))
#fam_tsv_l = ["Family_108.tsv"] #debugging.
l.log(log_stream, "# pedigrees: " + str(len(fam_tsv_l)))
fam_not_in_kinship_coef_df_l = [] 
for fam_tsv in fam_tsv_l:
    '''Read in pedigree and get a list of unique individuals in the family.'''
    fam = re.search('Family_(.*?)\.tsv', os.path.basename(fam_tsv)).group(1)
    l.log(log_stream, "FAMILY FILE: " + fam)
    if fam == "243" or re.search(r'article|S_|S-|poster|nuclear|old',fam) != None:
        l.log(log_stream, "Ignoring...")
        continue
    fam_df = pd.read_csv(os.path.join(args.pedigree_dir,"Family_"+fam+".tsv"),sep="\t", skiprows=1, header=None, usecols=[0,1,2,3], na_values="0", names=["FAMILY","PERSON","FATHER","MOTHER"], dtype="str")
    fam = fam_df.ix[0,"FAMILY"] #In case family is not the filename.
    l.log(log_stream, "FAMILY NAME: " + fam)
    uniq_ind_l = pd.unique(fam_df["PERSON"].tolist())
    uniq_ind_l = filter(lambda x: x != "0" and x.startswith("p") == False, uniq_ind_l)

    '''Get each individual's parents and grandparents'''
    s_p_gp_df = pd.DataFrame(data={"IND":uniq_ind_l})
    s_p_gp_df[["S","P","GP"]] = s_p_gp_df["IND"].apply(get_s_p_gp, fam_df=fam_df)
    s_p_gp_df.set_index("IND", inplace=True)
    #print s_p_gp_df

    '''Get every pair of individuals and store their expected and observed relationship.'''
    ind1_l, ind2_l = [],[]
    for subset in itertools.combinations(uniq_ind_l, 2):
        ind1_l.append(subset[0])
        ind2_l.append(subset[1])
    eo_df = pd.DataFrame({"FID1":[fam]*len(ind1_l),"ID1":ind1_l,"FID2":[fam]*len(ind2_l),"ID2":ind2_l})
    #print eo_df
    eo_df["EXP_REL"] = eo_df.apply(axis=1, func=get_exp_rel, s_p_gp_df=s_p_gp_df)
    eo_df = eo_df[pd.notnull(eo_df["EXP_REL"])]
    eo_df["OBS_REL"] = eo_df.apply(axis=1, func=get_obs_rel, obs_rel_df=kinship_coef_df)
    #print eo_df
    eo_df = eo_df[pd.notnull(eo_df["OBS_REL"])]
    #print eo_df
    if eo_df.empty:
        l.log(log_stream, "Couldn't find family " + fam + " in kinship results.")
        fam_not_in_kinship_coef_df_l.append(fam)

    '''For each individual, count # times EXP_REL == OBS_REL.'''
    results_dict = {"FID":[],"ID":[],"rel_n":[],"rel_correct_n":[],"correct_pc":[],"errors":[]}
    for ind in pd.unique(eo_df["ID1"].tolist() + eo_df["ID2"].tolist()):
        ind_eo_df = eo_df.ix[(eo_df["ID1"]==ind) | (eo_df["ID2"]==ind),:]
        rel_n = len(ind_eo_df.index)
        rel_correct_n = len(ind_eo_df.ix[ind_eo_df["OBS_REL"]==ind_eo_df["EXP_REL"],:].index)
        correct_pc = rel_correct_n/float(rel_n)
        results_dict["FID"].append(fam)
        results_dict["ID"].append(ind)
        results_dict["rel_n"].append(rel_n)
        results_dict["rel_correct_n"].append(rel_correct_n)
        results_dict["correct_pc"].append(correct_pc)
        results_dict["errors"].append(get_errors(ind_eo_df,ind))
    results_df_l.append(pd.DataFrame(results_dict))
    
    l.log(log_stream, "")

l.log(log_stream, "Merging results for each family and writing them to file...")
results_df = pd.concat(results_df_l)
results_df["sample"] = results_df[["FID","ID"]].apply(axis=1, func=lambda x: "_".join(x.tolist()))
results_df.drop(["FID","ID"], axis=1, inplace=True)
results_df = results_df.reindex(columns=["sample","rel_correct_n","rel_n","correct_pc","errors"])
results_df.set_index("sample", inplace=True)
results_df.sort_values(by="correct_pc", inplace=True)
results_df.to_csv(os.path.join(args.data_dir,".".join(["all","rel","familial","exp_v_obs","csv"])), index=True, index_label="sample")

l.log(log_stream, "Families with pedigrees & not in kinship results (" + str(len(fam_not_in_kinship_coef_df_l)) + "): " + ",".join(fam_not_in_kinship_coef_df_l))
