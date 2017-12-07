import pandas as pd
import numpy as np
import sys
import os
from gene_burden.gene_burden import CMC


'''Read in the genotype file.'''
data_dir = os.path.abspath(os.path.join("..","..","data","gene_burden"))
#print data_dir
geno_df = pd.read_csv(os.path.join(data_dir,".".join(["16","3","5","csv"])), dtype=str)
#print geno_df.ix[geno_df["SYMBOL"]=="NOD2",["Gene","Feature","SYMBOL"]]
#sys.exit()

'''Get the samples.'''
sample_l = geno_df.columns.tolist()
sample_l = sample_l[sample_l.index("hardy_wein_p_all")+1:]
#print sample_l
#print len(sample_l)
sample_dict = {}
sample_dict["all"] = sample_l
#Need to work out which samples are cases and which controls.
print sample_l
case_l = filter(lambda x: x.startswith(("Levine_","BGI_Feb2017","BGI_Sept2015")), sample_l)
control_l = filter(lambda x: x not in case_l, sample_l)

'''Make the gene column '''
geno_df[["CHROM","POS","REF","ALT"]] = geno_df["VARIANT_ID"].str.extract('(?P<CHROM>.*?)_(?P<POS>.*?)_(?P<REF>.*?)_(?P<ALT>.*)', expand=False)
geno_df["gene_transcript"] = geno_df[["CHROM","Gene","Feature","SYMBOL"]].apply(lambda row_s: "_".join(row_s.tolist()), axis=1)
geno_df.drop(["CHROM","POS","REF","ALT"], inplace=True, axis=1)

'''Do the gene burden.'''
cmc = CMC(sample_dict=sample_dict, gene_col="gene_transcript", frq_col_l=["BroadAJcontrols_ALT","gnomad_AF_NFE","EXAC_NFE"])
geno_collapsed_df = cmc.get_geno_collapsed_df(geno_df)
print geno_collapsed_df
print len(geno_collapsed_df.index)


#from scipy.stats import chisqprob
#import statsmodels.api as sm
#from statsmodels.discrete.discrete_model import Logit
import statsmodels.discrete.discrete_model as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

#def likelihood_ratio(llmin, llmax):
    #return(2*(llmax-llmin))


#Try NOD2
print geno_collapsed_df.ix["16_ENSG00000167207_ENST00000300589_NOD2"]
print geno_collapsed_df.index

#print geno_collapsed_df[]
#print geno_collapsed_df[]

print [1 if sample in case_l else 0 for sample in sample_l]
y = np.array([1 if sample in case_l else 0 for sample in sample_l])
print y.shape
X = np.transpose(geno_collapsed_df.ix["16_ENSG00000167207_ENST00000300589_NOD2",sample_l].values)
print X
print X.shape

logit_model=sm.Logit(y,X)
result=logit_model.fit()
print(result.summary())
print result.params
print result.llr_pvalue

#LR = likelihood_ratio(L1,L2)
#p = chisqprob(LR, 1) # L2 has 1 DoF more than L1
#print 'p: %.30f' % p 


#def do_multivariate_test(geno_collapsed_):

#result_df = geno_collapsed_df.apply(do_multivariate_test, axis=1)

#geno_collapsed_df[]

def do_multivariate_test(geno_collapsed_gene_df, sample_l, y):
    
    print geno_collapsed_gene_df.name
    X = np.transpose(geno_collapsed_gene_df[sample_l].values)
    print np.all(X==0)
    logit_model=sm.Logit(y,X)
    result=logit_model.fit()
    return result.llr_pvalue

result_s = geno_collapsed_df.groupby(level=['gene_transcript']).apply(do_multivariate_test, sample_l=sample_l, y=y)
print result_s
