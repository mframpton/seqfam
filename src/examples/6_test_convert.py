import pandas as pd
import sys
import numpy as np
import os
import argparse
import burden as b
import general as g
import log as l
import parent_args as pa


'''Read in arguments, set input directories and create log file.'''
parser = argparse.ArgumentParser(description="Filter variants and prepare files for gene burden tests.", parents=[pa.get_parent_parser()])
args = parser.parse_args()
pa.add_path_args(args)
frq_grp_name_l = args.frq_grp_names.split(",")
frq_grp_range_ll = [[float(frq) for frq in frq_range.split("-")] for frq_range in args.frq_grp_ranges.split(",")]

'''Get the unrelated geno_df'''
[log_stream, geno_df, sample_dict, caco_s] = b.get_unrelated_geno_df(args)
#print geno_df.head()

'''Collapsing: within each gene-transcript, collapse the rare variants.'''
l.log(log_stream, "Collapsing variants...")
l.log(log_stream, str(frq_grp_name_l))
l.log(log_stream, str(frq_grp_range_ll))
drop_col_l = filter(lambda x: x not in ["CHROM","SYMBOL","Gene","Feature","BroadAJcontrols_ALT","gnomad_AF_NFE","EXAC_AF_NFE"] + sample_dict["all"], geno_df.columns.tolist())
geno_df.drop(drop_col_l, axis=1, inplace=True)
geno_df = b.collapse_vars_in_frq_ranges(geno_df=geno_df, sample_l=sample_dict["all"], frq_grp_name_l=frq_grp_name_l, frq_grp_range_ll=frq_grp_range_ll)
#print geno_df.head()
geno_df.reset_index(inplace=True)
geno_df = geno_df.append(caco_s)
geno_df = geno_df.reindex(columns=["SYMBOL","CHROM","Gene","Feature","Collapsed_var","n"] + sample_dict["all"])
#print geno_df.shape
[chunk_num, chunk_denom] = [i for i in args.chunk.split("/")]
geno_df.to_csv(os.path.join(args.data_dir,".".join([str(args.chr),chunk_num,chunk_denom,"collapsed.csv.gz"])), compression="gzip", index=False)
l.log(log_stream, "FINISHED!")

