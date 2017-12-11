import pandas as pd
import numpy as np
import sys
import os
from cohpy.gene_burden import CMC


'''Read in the genotype file.'''
data_dir = os.path.abspath(os.path.join("..","..","data","gene_burden"))
geno_df = pd.read_csv(os.path.join(data_dir,".".join(["16","3","5","csv"])), dtype=str)

'''Get the samples.'''
sample_l = geno_df.columns.tolist()
sample_l = sample_l[sample_l.index("hardy_wein_p_all")+1:]
sample_dict = {}
sample_dict["all"] = sample_l
case_l = filter(lambda x: x.startswith(("Levine_","BGI_Feb2017","BGI_Sept2015")), sample_l)
control_l = filter(lambda x: x not in case_l, sample_l)
sample_dict["case"] = case_l
sample_dict["control"] = control_l

'''Make the gene column '''
geno_df[["CHROM","POS","REF","ALT"]] = geno_df["VARIANT_ID"].str.extract('(?P<CHROM>.*?)_(?P<POS>.*?)_(?P<REF>.*?)_(?P<ALT>.*)', expand=False)
geno_df["gene_transcript"] = geno_df[["CHROM","Gene","Feature","SYMBOL"]].apply(lambda row_s: "_".join(row_s.tolist()), axis=1)
geno_df.drop(["CHROM","POS","REF","ALT"], inplace=True, axis=1)
geno_df.set_index("VARIANT_ID", inplace=True)

'''Do the gene burden.'''
cmc = CMC(sample_dict=sample_dict, gene_col="gene_transcript", frq_col_l=["BroadAJcontrols_ALT","gnomad_AF_NFE","EXAC_NFE"])
#cmc_result_s = cmc.do_multivariate_tests(geno_df)
#print cmc_result_s.sort_values()

'''Read in the ancestry covariates.'''
covar_df = pd.read_csv(os.path.join(data_dir,"all.ancestry_results.csv"))
covar_df = covar_df.ix[:,['UCLex', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5']]
covar_df.set_index("UCLex", inplace=True)
covar_df = covar_df.ix[sample_dict["all"],:]
covar_df = covar_df.transpose()

cmc_result_df = cmc.do_multivariate_tests(geno_df, covar_df)
print cmc_result_df.sort_values(by="llr_cov_p")

