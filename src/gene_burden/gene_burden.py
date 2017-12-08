import pandas as pd
import numpy as np
import statsmodels.discrete.discrete_model as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)


class CMC(object):
    
    def __init__(self, sample_dict, gene_col, frq_col_l, frq_grp_name_l=["rare","mod_rare"], frq_grp_range_ll=[[0.0,0.01],[0.01,0.05]]):
        
        self.sample_dict = sample_dict
        self.gene_col = gene_col
        self.frq_col_l = frq_col_l
        self.frq_grp_name_l=frq_grp_name_l
        self.frq_grp_range_ll=frq_grp_range_ll
    
    
    def do_multivariate_tests(self, geno_df):
       
        geno_collapsed_df = self.get_geno_collapsed_df(geno_df)
        y = np.array([1 if sample in self.sample_dict["case"] else 0 for sample in self.sample_dict["all"]])
        result_s = geno_collapsed_df[self.sample_dict["all"]].groupby(level=[self.gene_col]).apply(self.do_multivariate_test, y=y)
        return result_s
    
       
    def do_multivariate_test(self, geno_collapsed_gene_df, y):
       
        #print geno_collapsed_gene_df.name
        X = np.transpose(geno_collapsed_gene_df.values)
        #print np.all(X==0)
        logit_model=sm.Logit(y,X)
        result=logit_model.fit(method='bfgs')
        return result.llr_pvalue
       
       
    def get_geno_collapsed_df(self, geno_df):
        
        geno_df[self.frq_col_l] = geno_df[self.frq_col_l].apply(pd.to_numeric, axis=1)
        geno_collapsed_df = self.collapse_vars_in_frq_ranges(geno_df=geno_df)
        return geno_collapsed_df
    
        
    def collapse_vars_in_frq_ranges(self, geno_df):

        geno_df["collapsed_var"] = geno_df.apply(axis=1, func=self.assign_collapse_group)
        collapse_group_n_df = geno_df.groupby(self.gene_col)["collapsed_var"].value_counts().to_frame("n")
        geno_collapsed_df = geno_df.groupby([self.gene_col] + ["collapsed_var"]).apply(func=self.collapse_vars_in_samples)
        geno_collapsed_df = geno_collapsed_df.merge(collapse_group_n_df, left_index=True, right_index=True)
        return geno_collapsed_df


    def assign_collapse_group(self, geno_s):
        '''Return each variant's collapse group. If the variant is not within a collapse group, just return the variant ID.'''
        
        var_frq = filter(lambda x: np.isnan(x) == False, geno_s.ix[self.frq_col_l].tolist()+[0.0])[0]
        for i in xrange(len(self.frq_grp_name_l)):
            if var_frq >= self.frq_grp_range_ll[i][0] and var_frq <= self.frq_grp_range_ll[i][1]:
                return self.frq_grp_name_l[i]
        return geno_s.name


    def collapse_vars_in_samples(self, gene_trans_group_df):   
        
        return gene_trans_group_df.ix[:,self.sample_dict["all"]].apply(axis=0, func=self.collapse_vars_in_1_sample)


    def collapse_vars_in_1_sample(self, sample_geno_s):
        
        geno_count_s = sample_geno_s.value_counts()
        if "1" in geno_count_s.index and geno_count_s.ix["1"] > 0:
                return 1
        elif "2" in geno_count_s.index and geno_count_s.ix["2"] > 0:
                return 1
        else:
                return 0
