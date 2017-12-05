import pandas as pd
import sys
import os
import argparse
import parent_args as pa
import log as l
import sample as s


def find_dups_in_king_file(king_file, index_col, sample_l):

    kinship_coef_df = pd.read_csv(king_file, sep="\t", index_col=index_col, dtype=str)
    kinship_coef_df["Kinship"] = pd.to_numeric(kinship_coef_df["Kinship"])
    duplicate_df = kinship_coef_df.ix[kinship_coef_df["Kinship"] > 0.354,:]
    del kinship_coef_df
    duplicate_l = []
    if len(duplicate_df.index) > 0:
            l.log(log_stream,"\n"+duplicate_df.to_string())
            duplicate_l = [convert_tuple_to_ids(duplicate_pair_tuple,sample_l) for duplicate_pair_tuple in duplicate_df.index.tolist()]
    else:
            l.log(log_stream,"0 potential duplicates.")

    return duplicate_l


def convert_tuple_to_ids(duplicate_pair_tuple, sample_l):
    #print duplicate_pair_tuple
    id1,id2 = "",""
    if len(duplicate_pair_tuple) == 3:
        id1,id2 = "_".join([duplicate_pair_tuple[0],duplicate_pair_tuple[1]])+"_","_".join([duplicate_pair_tuple[0],duplicate_pair_tuple[2]])+"_"
    elif len(duplicate_pair_tuple) == 4:
        id1,id2 = "_".join([duplicate_pair_tuple[0],duplicate_pair_tuple[1]])+"_","_".join([duplicate_pair_tuple[2],duplicate_pair_tuple[3]])+"_"
    full_id_l = filter(lambda x: x.startswith(tuple([id1,id2])),sample_l)
    return full_id_l


'''Set paths'''
parser = argparse.ArgumentParser(description="Find duplicates in the King results.", parents=[pa.get_parent_parser()])
args = parser.parse_args()
pa.add_path_args(args)
log_stream = open(os.path.join(args.data_dir,".".join(["all","rel","dups","log","txt"])),'w')

'''Read in the AJ samples we want to keep.'''
[sample_l, CD_l, UC_l, other_A_l, IBD_l, N_l] = s.get_familial_sample_ll(log_stream, args, uclex_id_b=False)

'''Between-family.'''
l.log(log_stream, "Checking for between-family duplicates...")
duplicate_bf_l = find_dups_in_king_file(os.path.join(args.data_dir,".".join(["all","rel","kinship","ibs0"])), ["FID1","ID1","FID2","ID2"], sample_l)
'''Within-family.'''
l.log(log_stream, "Checking for within-family duplicates...")
duplicate_wf_l = find_dups_in_king_file(os.path.join(args.data_dir,".".join(["all","rel","kinship","ibs"])), ["FID","ID1","ID2"], sample_l)

'''Write all duplicates.'''
duplicate_l = duplicate_bf_l + duplicate_wf_l
if len(duplicate_l) > 0:
    duplicate_df = pd.DataFrame(duplicate_l, columns=["sample1","sample2"])
    duplicate_flat_l = [item for sublist in duplicate_l for item in sublist]
    new_fam_id_uclex_dict = s.get_new_family_uclex_id_dict(log_stream, duplicate_flat_l, args)
    duplicate_df["uclex1"] = duplicate_df.apply(lambda row_s: new_fam_id_uclex_dict[row_s.ix["sample1"]], axis=1)
    duplicate_df["uclex2"] = duplicate_df.apply(lambda row_s: new_fam_id_uclex_dict[row_s.ix["sample2"]], axis=1)
    duplicate_path = os.path.join(args.shared_dir,"duplicates.csv")
    l.log(log_stream, "Writing {0}".format(duplicate_path))
duplicate_df.to_csv(duplicate_path,index=False)
