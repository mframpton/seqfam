import pandas as pd
import os
import sys
import numpy as np
import glob
import re
import filter as f
import general as g
import log as l
import sample as s
from bisect import *
from natsort import natsorted, ns


def write_ped(log_stream, geno_df, out_prefix, sample_dict, args, add_unsequenced_fam_members_b=True, write_geno_file=True, write_map_file=True, map_genetic_dist_b=True, map_physical_pos_b=False, write_dat_file=True, write_freq_file=True):

    '''Read in the manifest.'''
    manifest_df = s.get_manifest_df(log_stream, args)
    manifest_df["Sex"].replace({"M":"1","F":"2"}, inplace=True)
    manifest_df["Affection"] = manifest_df["Phenotype"].apply(lambda x: "1" if x == "N" else "2")#Treat IBD as affected.
    manifest_df.drop("Phenotype", axis=1, inplace=True)
    manifest_df.set_index("UCLex", inplace=True)
    manifest_df["old_ID"] = manifest_df.apply(axis=1, func=lambda x: "_".join([x["Family"],x["ID"]])) 
    '''If the analysis includes members of families A or B, then break these large families into sub-families in manifest_df.'''
    l.log(log_stream, "Dimensions of manifest_df: " + str(manifest_df.shape))
    if args.families == "all":
        manifest_df = break_AB_into_subfamilies(log_stream, manifest_df, args.pedigree_dir)
        l.log(log_stream, "Dimensions of manifest_df: " + str(manifest_df.shape))

    '''Convert the UCLex IDs in geno_df to the old IDs.'''
    l.log(log_stream, "Renaming sample IDs in geno_df...")
    uclex_old_id_dict =  manifest_df["old_ID"].to_dict()
    geno_df.rename(columns=uclex_old_id_dict, inplace=True)
    
    '''Write the geno file.'''
    if write_geno_file:
        geno_df.to_csv(os.path.join(args.data_dir,".".join([out_prefix,"csv"])))
    l.log(log_stream, "Dimensions of geno_df: " + str(geno_df.shape))

    '''Retain only the VARIANT identifier and genotype columns.'''
    l.log(log_stream, "Creating ped_df by retaining only the variant identifier and intersection of samples in geno_df and manifest_df...")
    geno_col_l = sorted([col for col in geno_df.columns.tolist() if col in uclex_old_id_dict.values()], key=g.natural_keys)
    ped_df = geno_df[geno_col_l]
    l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))
    '''Retain only SNPs.'''
    l.log(log_stream, "Retaining only the SNPs...")
    ped_df = retain_snps(ped_df)
    l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))
    '''Convert genotypes to ped format.'''
    l.log(log_stream, "Converting the genotypes to ped format...")
    ped_df = convert_genos(ped_df, geno_col_l)
    #f.check_geno_df_sorted(log_stream, ped_df)
    
    '''Transpose so index is samples.'''
    l.log(log_stream, "Transposing so index is samples and columns are variants.")
    ped_df = ped_df.transpose()
    ped_df.reset_index(inplace=True)
    ped_df.rename(columns={"index":"old_ID"}, inplace=True)
    l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))    
    #check_ped_cols_in_order(log_stream, ped_df)

    '''Get sex and affection columns from the manifest_df.'''
    l.log(log_stream, "Merging ped_df and manifest_df...")
    manifest_df.reset_index(inplace=True)
    ped_df = ped_df.merge(manifest_df[["old_ID","Sex","Affection"]], how="left")
    l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))    
    #check_ped_cols_in_order(log_stream, ped_df)

    '''Get father and mother columns from the pedigree_df.'''
    pedigree_df = get_pedigree_df(log_stream, args.pedigree_dir, ped_df)
    ped_df[["Family","Individual"]] = ped_df["old_ID"].str.extract('(?P<Family>.*?)_(?P<Individual>.*)', expand=False)
    ped_df.drop("old_ID", axis=1, inplace=True)
    ped_df[["Father","Mother"]] = ped_df[["Family","Individual"]].apply(axis=1, func=get_parents, pedigree_df=pedigree_df)
    ped_df.set_index(["Family","Individual"], inplace=True)
    l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))
    #check_ped_cols_in_order(log_stream, ped_df)

    if add_unsequenced_fam_members_b == True:
        '''Create a ped_df for the unsequenced individuals using the information in the pedigree_df '''
        l.log(log_stream, "Getting ped_df for unsequenced individuals...")
        ped_unseq_df = get_ped_unseq_df(log_stream, ped_df, pedigree_df)
        l.log(log_stream, "Dimensions of ped_unseq_df: " + str(ped_unseq_df.shape))    
        #check_ped_cols_in_order(log_stream, ped_unseq_df)

        '''Merge the ped_dfs for sequenced and unsequenced individuals.'''
        l.log(log_stream, "Merging ped_df and ped_unseq_df...")
        ped_df_col_l = ped_df.columns.tolist()
        ped_df = pd.concat([ped_df, ped_unseq_df])
        ped_df = ped_df.reindex(columns=ped_df_col_l) #concat does not preserve column order.
        l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))
        #check_ped_cols_in_order(log_stream, ped_df)

    '''Write the ped, map and dat files.'''
    l.log(log_stream, "Indexing ped file by family, individual, father, mother, sex and affection...")
    ped_df.reset_index(inplace=True)
    l.log(log_stream, "# families: " + str(len(pd.unique(ped_df["Family"]))))
    l.log(log_stream, "Sex value counts: ")
    l.log(log_stream, "\n" + ped_df["Sex"].value_counts().to_string())
        l.log(log_stream, "Affection value counts: ")
        l.log(log_stream, "\n" + ped_df["Affection"].value_counts().to_string())
    ped_df = remove_alt_names(ped_df)
    ped_df.set_index(["Family","Individual","Father","Mother","Sex","Affection"], inplace=True)
    if write_map_file:
        ped_df = write_map(log_stream, ped_df, args.data_dir, out_prefix, map_genetic_dist_b, None, map_physical_pos_b)
    if write_dat_file:
        write_dat(log_stream, ped_df, args.data_dir, out_prefix)
    if write_freq_file:
        write_freq(log_stream, ped_df, geno_df, args.data_dir, out_prefix)
    l.log(log_stream, "Dimensions of ped_df: " + str(ped_df.shape))
    ped_df.to_csv(os.path.join(args.data_dir,".".join([out_prefix,"ped"])), header=False, index=True, sep="\t")    


def break_AB_into_subfamilies(log_stream, manifest_df, pedigree_dir):

    l.log(log_stream, "Getting pedigree_df for families in ped_df...")
        fam_tsv_l = glob.glob(os.path.join(pedigree_dir,"Family_*.tsv"))
    fam_tsv_l = filter(lambda x: re.match(r'^Family_(A|B)', os.path.basename(x)) != None, fam_tsv_l)
    pedigree_df = pd.concat([pd.read_csv(fam_tsv_l[i],sep="\t", skiprows=1, header=None, usecols=[0,1], names=["Family","Individual"], dtype="str") for i in xrange(len(fam_tsv_l))])
    pedigree_df.set_index("Individual", inplace=True)
    pedigree_df = pedigree_df[~pedigree_df.index.duplicated(keep='first')]
    manifest_df["old_ID"] = manifest_df["old_ID"].apply(replace_fam_with_subfam, pedigree_df=pedigree_df)
    l.log(log_stream, "Removing members of families A and B who are not in one of the sub-families.")
    manifest_df = manifest_df.ix[manifest_df["old_ID"].str.contains(r'^(A|B)_')==False,]    

    return manifest_df


def replace_fam_with_subfam(fam_id, pedigree_df):

    if re.match(r'(A|B)', fam_id) != None:
        #print fam_id
        fam_id_l = fam_id.split("_")
        #print fam_id_l
        if fam_id_l[1] in pedigree_df.index:
            fam_id_l[0] = pedigree_df.ix[fam_id_l[1],"Family"]
        #else:
            #print "WARNING: DID NOT FIND IN PEDIGREE_DF."
        #print fam_id_l
        return "_".join(fam_id_l)
    else:
        return fam_id


def check_ped_cols_in_order(log_stream, ped_df):

    ped_col_l = ped_df.columns.tolist()
    #print ped_col_l
    ped_col_l = filter(lambda x: re.match("(\d*?)_(\d*?)_(\w*?)_(\w*)",x) != None, ped_col_l)
    #print ped_col_l
    pos_l = [int(re.match("(\d*?)_(\d*?)_(\w*?)_(\w*)",col).group(2)) for col in ped_col_l]
    #print pos_l
    l.log(log_stream, "ped columns are in order: " + str(np.all(np.diff(pos_l) >= 0)))    


def remove_alt_names(ped_df):

    ped_df["Individual"] = ped_df["Individual"].str.split("/").str.get(0)
    ped_df["Father"] = ped_df["Father"].str.split("/").str.get(0)
    ped_df["Mother"] = ped_df["Mother"].str.split("/").str.get(0)
    
    return ped_df


def get_chr_pos_ref_alt(ped_df):

    ped_df.reset_index(inplace=True)
        ped_df[["CHROM","POS","REF","ALT"]] = ped_df["VARIANT_ID"].str.extract('(?P<CHROM>.*?)_(?P<POS>.*?)_(?P<REF>.*?)_(?P<ALT>.*)', expand=False)
    
    return ped_df


def retain_snps(ped_df):

    ped_df = get_chr_pos_ref_alt(ped_df)
    ped_df["SNP"] = ped_df.apply(axis=1, func=lambda x: len(x["REF"])==1 and len(x["ALT"])==1)
    ped_df = ped_df[ped_df["SNP"]==True]
    ped_df.drop(["CHROM","POS","REF","ALT","SNP"], axis=1, inplace=True)
    ped_df.set_index("VARIANT_ID", inplace=True)

    return ped_df    


def convert_genos(ped_df, geno_col_l):

        ped_df = get_chr_pos_ref_alt(ped_df)
        ped_df[geno_col_l] = ped_df.apply(axis=1, func=get_geno_converted_s, geno_col_l=geno_col_l)
        ped_df.drop(["CHROM","POS","REF","ALT"], axis=1, inplace=True)
        ped_df.set_index("VARIANT_ID", inplace=True)

        return ped_df


def get_geno_converted_s(geno_s, geno_col_l):

        ref = geno_s["REF"]
        alt = geno_s["ALT"]
        geno_converted_s = geno_s[geno_col_l].map(lambda x: " ".join([alt,alt]) if x=="2" else " ".join([ref,alt]) if x=="1" else " ".join([ref,ref]) if x=="0" else "0 0")

        return geno_converted_s


def get_ped_unseq_df(log_stream, ped_df, pedigree_df):

    #ped_df.set_index(["Family","Individual"], inplace=True)
    ped_unseq_df = pedigree_df.ix[~pedigree_df.index.isin(ped_df.index),:]
    for var in filter(lambda x: x not in ["Family","Individual","Father","Mother","Sex","Affection"], ped_df.columns.tolist()):
        ped_unseq_df[var] = "0 0"
    
    return ped_unseq_df


def get_pedigree_df(log_stream, pedigree_dir, ped_df):

        l.log(log_stream, "Getting pedigree_df for families in ped_df...")
        fam_tsv_l = glob.glob(os.path.join(pedigree_dir,"Family_*.tsv"))
    '''Only read in the pedigrees for families present in ped_df.'''
    fam_in_ped_df_l = pd.unique(ped_df["old_ID"].str.split("_").str.get(0)).tolist()
    fam_tsv_l = filter(lambda x: re.search('_(.*?)\.tsv$' , os.path.basename(x)).group(1) in fam_in_ped_df_l, fam_tsv_l)
        l.log(log_stream, "# of family pedigree files: " + str(len(fam_tsv_l)))
        pedigree_df = pd.concat([pd.read_csv(fam_tsv_l[i],sep="\t", skiprows=1, header=None, usecols=[0,1,2,3,4,5], names=["Family","Individual","Father","Mother","Sex","Affection"], dtype="str") for i in xrange(len(fam_tsv_l))])
        l.log(log_stream, "# unique families in pedigree_df: " + str(pd.unique(pedigree_df["Family"]).size))
        l.log(log_stream, "Unique families in pedigree_df: " + " ".join(natsorted(pd.unique(pedigree_df["Family"]).tolist())))
    pedigree_df.set_index(["Family","Individual"], inplace=True)
        pedigree_df[["Father","Mother","Sex","Affection"]].fillna("0", inplace=True)
        l.log(log_stream, "Dimensions of pedigree_df: " + str(pedigree_df.shape))
        l.log(log_stream, "Removing duplicate individuals from pedigree_df...")
        pedigree_df = pedigree_df.ix[~pedigree_df.index.duplicated()]
        l.log(log_stream, "Dimensions of pedigree_df: " + str(pedigree_df.shape))

        return pedigree_df


def write_map(log_stream, ped_df, data_dir, out_prefix, map_genetic_dist_b, min_cM_gap, map_physical_pos_b):
    '''Return the ped_df because SNPs may have been pruned from the ped_df.'''

    l.log(log_stream,"Creating map_df from ped_df...")
        map_df = ped_df.columns.to_series().str.split('_', expand=True)
        map_df.columns = ["CHROM","POS","REF","ALT"]
        map_df["ID"] = map_df.apply(axis=1, func=lambda x: "_".join([x["CHROM"],x["POS"],x["REF"],x["ALT"]]))
        map_df["POS"] = map_df["POS"].apply(pd.to_numeric, errors="coerce")
        #print map_df["POS"]
        chr = map_df["CHROM"].iloc[0]
        
    if map_genetic_dist_b == True:
        l.log(log_stream,"Reading in rutgers_df...")
        shared_dir = "/cluster/project9/IBDAJE/exome_analyses/data/March2017/shared"
            rutgers_df = pd.read_csv(os.path.join(shared_dir,"rutgers_map_v3a","RUMapv3a_B137_chr"+chr+".txt"), usecols=["Build37_start_physical_position","Sex_averaged_start_Hald_map_position"], index_col="Build37_start_physical_position", sep=" ")
            l.log(log_stream,"Dimensions of rutgers_df: " + str(rutgers_df.shape))
            l.log(log_stream,"Removing duplicate physical positions from rutgers_df...")
            rutgers_df = rutgers_df[~rutgers_df.index.duplicated()]
            l.log(log_stream,"Dimensions of rutgers_df: " + str(rutgers_df.shape))
            l.log(log_stream, "Getting genetic distances...")
            map_df["Genetic distance"] = map_df["POS"].apply(get_genetic_distance, rutgers_df=rutgers_df)
        '''Prune SNPs'''
        if min_cM_gap != None:
            [map_df, ped_df] = prune_snps(log_stream, ped_df, map_df, min_cM_gap)        
    else:
        l.log(log_stream,"Setting genetic distances to unknown (0).")
        map_df["Genetic distance"] = "0"
    if map_physical_pos_b == True:
        map_df = map_df.reindex(columns=["CHROM","ID","Genetic distance","POS"])
    else:
        map_df = map_df[["CHROM","ID","Genetic distance"]]    

    '''Write the map file.'''
        l.log(log_stream,"Writing map_df...")    
    map_df.to_csv(os.path.join(data_dir,".".join([out_prefix,"map"])), header=False, index=False, sep="\t")

    return ped_df


def prune_snps(log_stream, ped_df, map_df, min_cM_gap):

    '''Prune SNPs which are < min_cM_gap from the previous. '''
        pos_l = map_df["Genetic distance"].tolist()
        #print pos_l
        retained_idx_l = [0]
        for i in xrange(1,len(pos_l)):
            if pos_l[i] - pos_l[retained_idx_l[-1]] >= min_cM_gap:
                    retained_idx_l.append(i)
                          #print retained_idx_l
         retained_var_id_l = map_df.ix[retained_idx_l,"ID"].tolist()
           l.log(log_stream, "Pruning " + str(len(map_df.index) - len(retained_idx_l)) + " variants.")
          map_df = map_df.ix[retained_idx_l]
          ped_df = ped_df[retained_var_id_l]

    return [map_df, ped_df]
    
    
def get_genetic_distance(pos, rutgers_df):
    '''Get the genetic distance from the rutgers map, using linear interpolation where necessary.'''

    #print pos
    if pos in rutgers_df.index:
        #print "Not interpolating..."
        genetic_distance = rutgers_df.ix[pos,"Sex_averaged_start_Hald_map_position"]
        #print genetic_distance
    else:
        #print "Interpolating..."
        rutgers_pos_arr = rutgers_df.index.values
        if pos <= np.min(rutgers_pos_arr):
                        genetic_distance = rutgers_df["Sex_averaged_start_Hald_map_position"].min()
                elif pos >= np.max(rutgers_pos_arr):
                        genetic_distance = rutgers_df["Sex_averaged_start_Hald_map_position"].max()
        else: 
            lower = rutgers_pos_arr[bisect_left(rutgers_pos_arr, pos) - 1]
            above = rutgers_pos_arr[bisect_right(rutgers_pos_arr, pos)]
            #print lower,above
            genetic_distance =  np.interp(pos,[lower,above],[rutgers_df.ix[lower,"Sex_averaged_start_Hald_map_position"],rutgers_df.ix[above,"Sex_averaged_start_Hald_map_position"]])
            #print genetic_distance

    return genetic_distance

                        
def write_dat(log_stream, ped_df, data_dir, out_prefix):
    '''Write a dat (data) file to describe the pedigree (ped) file.'''
    '''Should contain: A CD\n M marker1\n etc.'''    

    l.log(log_stream,"Creating dat_df from ped_df...")
    dat_df = ped_df.columns.to_series().str.split('_', expand=True)    
    dat_df.index.name = "VARIANT_ID"
    dat_df.columns = ["CHROM","POS","REF","ALT"]
    dat_df["ID"] = dat_df.apply(axis=1, func=lambda x: "_".join([x["CHROM"],x["POS"],x["REF"],x["ALT"]]))
    dat_df.drop(["CHROM","POS","REF","ALT"], axis=1, inplace=True)
    dat_df["Type"] = "M"
    dat_df.reset_index(inplace=True)
    dat_df.drop("VARIANT_ID", axis=1, inplace=True)
    #print dat_df
    dat_df.loc[-1] = ["IBD","A"]
    dat_df.sort_index(inplace=True)
    dat_df = dat_df.reindex(columns=["Type","ID"])
    l.log(log_stream,"Writing dat_df...")
    dat_df.to_csv(os.path.join(data_dir,".".join([out_prefix,"dat"])), header=False, index=False, sep="\t")


def write_freq(log_stream, ped_df, geno_df, data_dir, out_prefix):

    #print geno_df.columns.tolist()
    freq_stream = open(os.path.join(data_dir,out_prefix+".freq"), 'w')
    var_id_l = ped_df.columns.tolist()
    for var_id in var_id_l:
        freq_stream.write(" ".join(["M",var_id]) + "\n")
        [ref, alt] = var_id.split("_")[2:]  
        alt_frq, ref_frq = geno_df.ix[var_id,"BroadAJcontrols_ALT"], 1-geno_df.ix[var_id,"BroadAJcontrols_ALT"]
        if alt_frq == 0.0:
            alt_frq, ref_frq = geno_df.ix[var_id,"gnomad_AF_NFE"], 1-geno_df.ix[var_id,"gnomad_AF_NFE"]
        if alt_frq == 0.0:
            alt_frq, ref_frq = geno_df.ix[var_id,"EXAC_NFE"], 1-geno_df.ix[var_id,"EXAC_NFE"]
        if alt_frq == 0.0:
            alt_frq, ref_frq = geno_df.ix[var_id,"Kaviar"], 1-geno_df.ix[var_id,"Kaviar"]
        if alt_frq == 0.0:
            alt_frq, ref_frq = 0.00001, 0.99999
        alt_frq = str(alt_frq)
        ref_frq = str(ref_frq)
        freq_stream.write(" ".join(["A",ref,ref_frq]) + "\n")
        freq_stream.write(" ".join(["A",alt,alt_frq]) + "\n")
    freq_stream.close()


def get_parents(id_s, pedigree_df):
    '''Find the parents in pedigree_df, but if these parents are not in ped_df, treat them as unknown.'''

    parent_s = pd.Series(["0","0"], index=["Father","Mother"])
    if tuple(id_s.tolist()) in pedigree_df.index:    
        parent_s = pedigree_df.ix[tuple(id_s.tolist()),["Father","Mother"]]

    return parent_s


def merge_ped_dat_map(in_dir, out_prefix):

    '''Create merged ped file.'''
    to_merge_l = sorted(filter(lambda x: "all" not in x, glob.glob(os.path.join(in_dir,"*.ped"))), key=g.natural_keys)
    #print to_merge_l
    ped_df = pd.concat([pd.read_csv(ped_tsv, sep="\t", header=None, index_col=[0,1,2,3,4,5]) for ped_tsv in to_merge_l], axis=1)
    ped_df.to_csv(os.path.join(in_dir, out_prefix + ".ped"), sep="\t", index=True, header=None)
    del ped_df

    '''Create merged dat file.'''
    to_merge_l = sorted(filter(lambda x: "all" not in x, glob.glob(os.path.join(in_dir,"*.dat"))), key=g.natural_keys)
    #print to_merge_l
    dat_df = pd.concat([pd.read_csv(ped_tsv, sep="\t", header=None) for ped_tsv in to_merge_l])
    dat_df.drop_duplicates(inplace=True)
    dat_df.to_csv(os.path.join(in_dir, out_prefix + ".dat"), sep="\t", index=False, header=None)
    del dat_df

    '''Create merged map file.'''
    to_merge_l = sorted(filter(lambda x: "all" not in x, glob.glob(os.path.join(in_dir,"*.map"))),key=g.natural_keys)
    #print to_merge_l
        map_df = pd.concat([pd.read_csv(map_tsv, sep="\t", header=None) for map_tsv in to_merge_l])
        map_df.to_csv(os.path.join(in_dir, out_prefix + ".map"), sep="\t", index=False, header=None)