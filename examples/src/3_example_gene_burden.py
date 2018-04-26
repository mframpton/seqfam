import os, sys
import pandas as pd
from seqfam.gene_burden import CMC


#Set paths.
data_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data","gene_burden"))
genotypes_path = os.path.join(data_dir,".".join(["variant_genotypes","csv"]))
samples_path = os.path.join(data_dir,".".join(["samples","csv"]))
covariates_path = os.path.join(data_dir,".".join(["covariates","csv"]))
results_path = os.path.join(data_dir,".".join(["cmc","results","csv"]))

#Read the samples into a Series.
sample_s = pd.read_csv(samples_path, dtype=str, index_col="Sample ID")
sample_s["Affection"] = sample_s["Affection"].astype(int)
sample_s = sample_s[sample_s != 0]

#Read the variant annotations + genotypes into a DataFrame.
variant_col,gene_col = "VARIANT_ID","Gene"
pop_frq_col_l = ["database1_AF","database2_AF","database3_AF"]
geno_df = pd.read_csv(genotypes_path, dtype=str, usecols=[variant_col,gene_col] + pop_frq_col_l + sample_s.index.tolist(), index_col=variant_col)
geno_df.loc[:,pop_frq_col_l] = geno_df.loc[:,pop_frq_col_l].apply(pd.to_numeric, axis=1)
geno_df.loc[:,sample_s.index] = geno_df.loc[:,sample_s.index].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
geno_df.loc[:,sample_s.index] = geno_df.loc[:,sample_s.index].fillna(0)

#Read the covariates into a DataFrame.
covar_df = None if covariates_path == None else pd.read_csv(covariates_path, index_col=0)

#Do gene burden (CMC) tests.
cmc = CMC()
pop_frq_cat_dict = {"rare":0.01,"mod_rare":0.05}
geno_df = cmc.assign_variants_to_pop_frq_cats(geno_df, pop_frq_col_l, pop_frq_cat_dict)
cmc_result_df = cmc.do_multivariate_tests(sample_s, geno_df, group_col=gene_col, agg_col="pop_frq_cat", agg_val_l=list(pop_frq_cat_dict.keys()), covar_df=covar_df, results_path=results_path)
#Uncomment the following 2 lines to run with the tests without the covariates.
#cmc_result_df = cmc.do_multivariate_tests(sample_s, geno_df, group_col=gene_col, agg_col="pop_frq_cat", agg_val_l=list(pop_frq_cat_dict.keys()), covar_df=None, results_path=results_path)
