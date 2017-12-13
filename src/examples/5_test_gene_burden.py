import os
from cohpy.gene_burden import CMC


'''Read in the genotype file.'''
data_dir = os.path.abspath(os.path.join("..","..","data","gene_burden"))
genotypes_path = os.path.join(data_dir,".".join(["16","3","5","formatted","csv"]))
samples_path = os.path.join(data_dir,".".join(["samples","csv"]))
covariates_path = os.path.join(data_dir,".".join(["covariates","csv"]))
results_path = os.path.join(data_dir,".".join(["cmc","results","csv"]))

'''Do gene burden (CMC) tests.'''
cmc = CMC()
cmc_result_df = cmc.do_multivariate_tests(samples_path, genotypes_path, covariates_path=covariates_path, results_path=results_path)
print cmc_result_df
