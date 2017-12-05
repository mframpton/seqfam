import pandas as pd
import os
import numpy as np
import sys
import log as l


def hello_world():
    print "hello_world"


def read_king_results(log_stream, data_dir):  
    '''Read in the observed (KING-inferred) relationships.'''
    infam_kinship_coef_df = pd.read_csv(os.path.join(data_dir,"all.rel.kinship.ibs"),sep="\t", dtype="str")
    infam_kinship_coef_df.rename(columns={"FID":"FID1"}, inplace=True)
    infam_kinship_coef_df["FID2"] = infam_kinship_coef_df["FID1"]
    keep_col_l = ["FID1","ID1","FID2","ID2","Kinship"]
    infam_kinship_coef_df = infam_kinship_coef_df.ix[:,keep_col_l]
    l.log(log_stream, "Dimensions of infam_kinship_coef_df: "  + str(infam_kinship_coef_df.shape))
    kinship_coef_df =  pd.read_csv(os.path.join(data_dir,"all.rel.kinship.ibs0"), sep="\t", usecols=keep_col_l, dtype="str")
    l.log(log_stream, "Dimensions of (between-family) kinship_coef_df: " + str(kinship_coef_df.shape))
    kinship_coef_df = pd.concat([kinship_coef_df,infam_kinship_coef_df])
    l.log(log_stream, "Dimensions of kinship_coef_df: " + str(kinship_coef_df.shape))
    del infam_kinship_coef_df
    kinship_coef_df["Kinship"] = pd.to_numeric(kinship_coef_df["Kinship"], errors="coerce")
    l.log(log_stream, "Creating FIDID1 column...")
    kinship_coef_df["FIDID1"] = kinship_coef_df[["FID1","ID1"]].apply(axis=1, func=lambda x: "_".join(x))
    l.log(log_stream, "Creating FIDID2 column...")
    kinship_coef_df["FIDID2"] = kinship_coef_df[["FID2","ID2"]].apply(axis=1, func=lambda x: "_".join(x))
    
    return kinship_coef_df


def remove_sporadics(log_stream, kinship_coef_df):
    '''Remove the sporadics.'''
    sporadic_l = get_sporadics(kinship_coef_df)
    l.log(log_stream, "# of sporadics: " + str(len(sporadic_l)))
    kinship_coef_df = kinship_coef_df.ix[~(kinship_coef_df["FIDID1"].isin(sporadic_l)) & ~(kinship_coef_df["FIDID2"].isin(sporadic_l)),:]
    l.log(log_stream, "Dimensions of kinship_coef_df: " + str(len(kinship_coef_df.index)))
    
    return kinship_coef_df


def get_sporadics(kinship_coef_df):
    '''Get the list of sporadics.'''    
    [fam_count_s, uniq_ind_l] = get_fam_count_s(kinship_coef_df)
    fam_l = fam_count_s[fam_count_s == 1].index.tolist()
    sporadic_l = [ind for ind in uniq_ind_l if ind.split("_")[0] in fam_l]
    
    return sporadic_l


def get_familials(kinship_coef_df):
    '''Get the list of familials.'''
    [fam_count_s, uniq_ind_l] = get_fam_count_s(kinship_coef_df)
    fam_l = fam_count_s[fam_count_s > 1].index.tolist()
    familial_l = [ind for ind in uniq_ind_l if ind.split("_")[0] in fam_l]

    return familial_l


def get_fam_count_s(kinship_coef_df):
    '''Get a series containing the family sizes.'''
    uniq_ind_l = pd.unique(kinship_coef_df["FIDID1"].tolist() + kinship_coef_df["FIDID2"].tolist())
    print "# of individuals: " + str(len(uniq_ind_l))
    fam_count_s = pd.Series(data=[ind.split("_")[0] for ind in uniq_ind_l]).value_counts()
    
    return [fam_count_s, uniq_ind_l]


def get_fam_score_s(ind, kinship_coef_df, exp_comparison_n):
    ''' '''
    keep_col_l = ["FID1","ID1","FID2","ID2","Kinship","Score"]
    [fid, iid] = ind.split("_")
    df = kinship_coef_df.ix[(kinship_coef_df["FID1"]==fid) & (kinship_coef_df["ID1"]==iid),keep_col_l]
    df2 = kinship_coef_df.ix[(kinship_coef_df["FID2"]==fid) & (kinship_coef_df["ID2"]==iid),keep_col_l]
    if df2.empty == False:
        df2.rename(columns={"FID1":"FID2","ID1":"ID2","FID2":"FID1","ID2":"ID1"},inplace=True)
        df = pd.concat([df,df2])
        del df2
    if df.shape[0] != exp_comparison_n:
        print "WARNING: " + str(df.shape[0]) + " rows."
    '''Aggregate the score by family.'''
    fam_score_s = df[["FID2","Score"]].groupby("FID2").aggregate(np.sum)["Score"]
    fam_score_s.name = ind
    
    return fam_score_s


def add_assigned_fam_cols(fam_score_df): 
    '''Add assigned family columns to a fam_score_df.'''
    fam_max_s = fam_score_df.apply(axis=1, func=lambda x: ",".join(x[(x>0) & (x==np.max(x))].index.tolist()))
    fam_max_s.name = "fam_max"
    fam_all_s = fam_score_df.apply(axis=1, func=lambda x: ",".join(x[x>0].index.tolist()))
    fam_all_s.name = "fam_all"
    fam_score_df = fam_score_df.join(fam_max_s)
    fam_score_df = fam_score_df.join(fam_all_s)
    
return fam_score_df
