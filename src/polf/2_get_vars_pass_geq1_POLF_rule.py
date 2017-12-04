import pandas as pd
import numpy as np
import sys
import os
from natsort import natsorted
import argparse
import sample as s
import log as l
import parent_args as pa


def get_lfam_pass_polf(row_s, family_summary_df, CD_b=False, A_n_min=0, N_n_min=0, N_n_max=np.inf, A_p=None, AN_p_diff_min=None):
	
	A_col = "CDs" if CD_b else "As"
	A_n_col = "CD_n" if CD_b else "A_n"
	lfam_summary_df = family_summary_df.ix[(family_summary_df[A_n_col]>=A_n_min) & (family_summary_df["N_n"]>=N_n_min) & (family_summary_df["N_n"]<=N_n_max),:]
	lfam_pass_l = []
	for lfam in lfam_summary_df.index.tolist():
		A_call_n = (row_s.ix[lfam_summary_df.ix[lfam,A_col]]!="NA").sum()
		N_call_n = (row_s.ix[lfam_summary_df.ix[lfam,"Ns"]]!="NA").sum()
		if A_call_n >= A_n_min and N_call_n >= N_n_min and N_call_n <= N_n_max:
			A_carrier_p = get_carrier_p(row_s,lfam_summary_df.ix[lfam,A_col])
			if A_p != None:
				if A_carrier_p >= A_p:
					lfam_pass_l.append(lfam)
			elif AN_p_diff_min != None:
				N_carrier_p = get_carrier_p(row_s,lfam_summary_df.ix[lfam,"Ns"])
				if A_carrier_p - N_carrier_p >= AN_p_diff_min:
					lfam_pass_l.append(lfam)

	[lfam_pass_str,lfam_pass_n] = get_fam_pass_str_n(lfam_pass_l)

        return pd.Series([lfam_pass_str,lfam_pass_n])


def get_fam_pass_str_n(fam_pass_l):
	
	fam_pass_str, fam_pass_n = np.NaN,0
	if len(fam_pass_l) > 0:
                fam_pass_str = ",".join(fam_pass_l)
                fam_pass_n = len(fam_pass_l)

	return [fam_pass_str,fam_pass_n]


def get_carrier_p(row_s,sample_l):

	if len(sample_l) == 0:
		return 0.0
	geno_value_count_s = row_s.ix[sample_l].value_counts(normalize=True, dropna=False)
	carrier_p = 0.0
	for carrier_geno in ["1","2"]:
		if carrier_geno in geno_value_count_s.index:
			carrier_p += geno_value_count_s.ix[carrier_geno]
	
	return carrier_p


def get_sfam_pass_polf(row_s, family_summary_df):

	sfam_pass_l = []
	for sfam in family_summary_df.index:
		if get_carrier_p(row_s,family_summary_df.ix[sfam,"As"]) > 0.0:
			sfam_pass_l.append(sfam)
	[sfam_pass_str,sfam_pass_n] = get_fam_pass_str_n(sfam_pass_l)

	return pd.Series([sfam_pass_str,sfam_pass_n])


def get_carrier_n(geno_s):
        
	geno_value_count_s = geno_s.value_counts()
        carrier_n = 0
        for carrier_geno in ["1","2"]:
                if carrier_geno in geno_value_count_s.index:
                        carrier_n += geno_value_count_s.ix[carrier_geno]
        return carrier_n


'''Read in UCLEx build argument and set the data directory.'''
parser = argparse.ArgumentParser(description="Merge the results.", parents=[pa.get_parent_parser()])
args = parser.parse_args()
pa.add_path_args(args)
chunk_l = args.chunk.split("/")
log_stream = open(os.path.join(args.data_dir,".".join([str(args.chr)]+chunk_l+["vrare","polf","log","txt"])), 'w')

'''Make family_summary_df.'''
[sample_l, CD_l, UC_l, other_A_l, IBD_l, N_l] = s.get_familial_sample_ll(log_stream, args, uclex_id_b=False, exclude_fam_l=["A","B","LC","SA","SC-1"], remove_dup_b=True)
#[sample_l,IBD_l,N_l,CD_l,UC_l] = s.add_family_79b_to_sample_l(sample_l,IBD_l,N_l,CD_l,UC_l)
family_summary_df = s.get_family_summary_df_wrapper(sample_l)
iso_fam_l = family_summary_df.ix[family_summary_df["A_n"]<=1,:].index.tolist()
iso_l = []
for iso_fam in iso_fam_l:
	iso_l.extend(family_summary_df.ix[iso_fam,"As"])
l.log(log_stream, "# isos: " + str(len(iso_l)))
lfam_l = family_summary_df.ix[(family_summary_df["A_n"]>=3) & (family_summary_df["N_n"]>=2),:].index.tolist() + family_summary_df.ix[(family_summary_df["A_n"]>=4) & (family_summary_df["N_n"]<=1),:].index.tolist()
#print natsorted(lfam_l)
l.log(log_stream, "# lfams: " + str(len(lfam_l)))
sfam_l = natsorted(list(set(family_summary_df.index.tolist()).difference(set(iso_fam_l+lfam_l))))
l.log(log_stream, "# sfams: " + str(len(sfam_l)))
fam_summary_csv = os.path.join(args.data_dir,".".join(["fam","summary","csv"]))
if os.path.exists(fam_summary_csv) == False:
	family_summary_df.ix[natsorted(lfam_l+sfam_l),:].to_csv(fam_summary_csv)

'''Read in variants.'''
l.log(log_stream, "Reading in the rare, damaging, functional variants...")
vrare_geno_df = pd.read_csv(os.path.join(args.data_dir,".".join([args.chr]+chunk_l+["vrare","csv","gz"])), index_col="VARIANT_ID", dtype=object)
col_keep_l = filter(lambda col: col.startswith(("A_","B_","LC_","SA_","SC-1_")) == False, vrare_geno_df.columns.tolist())
vrare_geno_df = vrare_geno_df.ix[:,col_keep_l]
#vrare_geno_df = s.add_family_79b_to_df(vrare_geno_df)
vrare_geno_df[sample_l] = vrare_geno_df[sample_l].fillna("NA")
l.log(log_stream, "Dimensions of vrare_geno_df: " + str(vrare_geno_df.shape))

'''Make columns which indicate for each variant, which lfams they pass a polf rule in.'''
l.log(log_stream, "Making A3N2 A columns...")
vrare_geno_df[["A3N2_As","A3N2_A_n"]] = vrare_geno_df.apply(get_lfam_pass_polf, family_summary_df=family_summary_df, CD_b=False, A_n_min=3, N_n_min=2, AN_p_diff_min=0.5, axis=1)
l.log(log_stream, "Making A3N2 CD columns...")
vrare_geno_df[["A3N2_CDs","A3N2_CD_n"]] = vrare_geno_df.apply(get_lfam_pass_polf, family_summary_df=family_summary_df, CD_b=True, A_n_min=3, N_n_min=2, AN_p_diff_min=0.5, axis=1)
l.log(log_stream, "Making A4N1 A columns...")
vrare_geno_df[["A4N1_As","A4N1_A_n"]] = vrare_geno_df.apply(get_lfam_pass_polf, family_summary_df=family_summary_df, CD_b=False, A_n_min=4, N_n_min=0, N_n_max=1, A_p=1.0, axis=1)
l.log(log_stream, "Making A4N1 CD columns...")
vrare_geno_df[["A4N1_CDs","A4N1_CD_n"]] = vrare_geno_df.apply(get_lfam_pass_polf, family_summary_df=family_summary_df, CD_b=True, A_n_min=4, N_n_min=0, N_n_max=1, A_p=1.0, axis=1)
vrare_geno_df[["A3N2_A_n","A3N2_CD_n","A4N1_A_n","A4N1_CD_n"]] = vrare_geno_df[["A3N2_A_n","A3N2_CD_n","A4N1_A_n","A4N1_CD_n"]].apply(pd.to_numeric, downcast='integer', axis=1)

'''Retain only variants which pass a polf rule in >= 1 lfam.'''
l.log(log_stream, "Retaining variants which pass a polf rule in >= 1 lfam...")
vrare_geno_df = vrare_geno_df.ix[(vrare_geno_df["A3N2_A_n"] > 0) | (vrare_geno_df["A3N2_CD_n"] > 0) | (vrare_geno_df["A4N1_A_n"] > 0) | (vrare_geno_df["A4N1_CD_n"] > 0),:]
l.log(log_stream, "Dimensions of vrare_geno_df: " + str(vrare_geno_df.shape))

if vrare_geno_df.empty == False:
	
	'''Add columns for sfams, isos, carrier_n'''
	l.log(log_stream, "Making sfam columns...")
	vrare_geno_df[["sfams","sfam_n"]] = vrare_geno_df.apply(get_sfam_pass_polf, family_summary_df=family_summary_df.ix[sfam_l,:], axis=1)
	l.log(log_stream, "Making iso_n column...")
	vrare_geno_df["iso_n"] = vrare_geno_df[iso_l].apply(get_carrier_n, axis=1)
	vrare_geno_df[["sfam_n","iso_n"]] = vrare_geno_df[["sfam_n","iso_n"]].apply(pd.to_numeric, downcast='integer', axis=1)
	l.log(log_stream, "Making carrier_n column...")
	vrare_geno_df["carrier_n"] = vrare_geno_df[sample_l].apply(get_carrier_n, axis=1)

	l.log(log_stream, "Moving sample columns to the end.")
	col_l = filter(lambda x: x not in sample_l, vrare_geno_df.columns.tolist()) + sample_l
	vrare_geno_df = vrare_geno_df.reindex(columns=col_l)
	l.log(log_stream, "Dimensions of vrare_geno_df: " + str(vrare_geno_df.shape))
	l.log(log_stream, "Writing vrare_geno_df.")
	vrare_geno_df.to_csv(os.path.join(args.data_dir,".".join([str(args.chr)]+chunk_l+["vrare","polf","csv"])))

