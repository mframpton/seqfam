import os
import pandas as pd
from cohpy.gene_burden import CMC
import sys


'''Set paths.'''
data_dir = os.path.abspath(os.path.join("..","..","data","gene_burden"))
genotypes_path = os.path.join(data_dir,".".join(["variant_genotypes.csv"]))
samples_path = os.path.join(data_dir,".".join(["samples","csv"]))
covariates_path = os.path.join(data_dir,".".join(["covariates","csv"]))
results_path = os.path.join(data_dir,".".join(["cmc","results","csv"]))

'''Read the samples into a Series'''
sample_s = pd.read_csv(samples_path, dtype=str, index_col="Sample ID")
sample_s["Affection"] = sample_s["Affection"].astype(int)
sample_s = sample_s[sample_s != 0]

'''Read the variant annotations + genotypes into a DataFrame '''
variant_col,gene_col = "VARIANT_ID","gene_transcript"
pop_frq_col_l = ["BroadAJcontrols_ALT","gnomad_AF_NFE","EXAC_NFE"]
geno_df = pd.read_csv(genotypes_path, dtype=str, usecols=[variant_col,gene_col] + pop_frq_col_l + sample_s.index.tolist(), index_col=variant_col)
geno_df[pop_frq_col_l] = geno_df[pop_frq_col_l].apply(pd.to_numeric, axis=1)
geno_df[sample_s.index] = geno_df[sample_s.index].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
geno_df[sample_s.index].fillna(0, inplace=True)

'''Read the covariates into a DataFrame.'''
covar_df = None if covariates_path == None else pd.read_csv(covariates_path, index_col=0)

'''Do gene burden (CMC) tests.'''
cmc = CMC()
geno_df = cmc.assign_variants_to_pop_frq_cats(geno_df, pop_frq_col_l, {"rare":0.01,"mod_rare":0.05})
cmc_result_df = cmc.do_multivariate_tests(sample_s, geno_df, group_col=gene_col, agg_col="pop_frq_cat", agg_val_l=["rare","mod_rare"], covar_df=covar_df, results_path=results_path)
print(cmc_result_df)
