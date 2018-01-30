import pandas as pd
import numpy as np
import statsmodels.discrete.discrete_model as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
import sys
from cohpy.misc import Logger
import operator
from collections import OrderedDict


class CMC(object):
    
    def __init__(self):
    
        self.logger = Logger()
        self.sample_s = None
        self.geno_df = None
        self.group_col = None
        self.agg_col = None
        self.agg_val_l = None
    
    
    def do_multivariate_tests(self, sample_s, geno_df, group_col, agg_col, agg_val_l, covar_df=None, results_path=None):
        '''Main method for doing multivariate tests.
        
        Args:
            sample_s (Series): contains samples and affection status.
            geno_df (DataFrame): contains variants annotations (the group and aggregate columns) and genotypes.
            group_col (str): column by which to group variants (1 test per group e.g. gene/functional unit).
            agg_col (str): column by which to aggregate variants e.g. allele population frequency.
            agg_val_l (list of strs): names of the aggregated categories.
            covar_df (DataFrame): contains covariates.
            results_path (str): path to results file.
        
        Returns:
            result_df (DataFrame): multivariate test results.
        '''
        
        #Set object attributes.
        self.group_col = group_col
        self.agg_col = agg_col
        self.agg_val_l = agg_val_l
        #Read in samples, genotypes and covariates.
        self.logger.log("Reading in samples and annotated genotypes...")
        self.sample_s = sample_s
        self.sample_s["Affection"] = self.sample_s["Affection"].astype(int)
        self.sample_s = self.sample_s[self.sample_s != 0]
        self.geno_df = geno_df
        geno_df[self.sample_s.index] = geno_df[self.sample_s.index].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
        geno_df[self.sample_s.index].fillna(0, inplace=True)
        self.covar_df = covar_df
        self.logger.log("# variants: {0}".format(len(geno_df.index)))
        self.logger.log("# genes to test: {0}".format(len(pd.unique(geno_df[self.group_col]))))
        #Aggregate the genotypes.
        self.logger.log("Aggregating genotypes...")
        geno_agg_df = self.aggregate_by_agg_col(geno_df)
        #Do the multivariate tests.
        self.logger.log("Doing multivariate tests...")
        result_df = geno_agg_df[self.sample_s.index].groupby(level=[self.group_col]).apply(self.do_multivariate_test, y=self.sample_s.values-1, covar_df=covar_df)
        #Merge the results with the collapsed variant counts.
        self.logger.log("Merge the results with the population frequency variant category counts...")
        agg_cat_count_df = self.get_agg_cat_count_df(geno_agg_df)
        result_df = agg_cat_count_df.join(result_df)
        if covar_df is None:
            result_df.sort_values(by="llr_p", inplace=True)
        else:
            result_df.sort_values(by="llr_cov_p", inplace=True)
        self.logger.log("Write results.")
        result_df.to_csv(results_path, index=True)
        return result_df
    
    
    def aggregate_by_agg_col(self, geno_df):
        '''Aggregate genotypes within variant population frequency categories.
        
        Args:
            geno_df (DataFrame): contains the sample genotypes, variant identifier, gene identifier and population frequencies.
            
        Returns:
            geno_agg_df (DataFrame): contains the genotypes aggregated within variant population frequency categories.
        '''

        self.logger.log("For each gene, count the # of variants in each aggregation category.")
        agg_cat_n_s = geno_df.groupby(self.group_col)[self.agg_col].value_counts()
        agg_cat_n_s.name = "n"
        self.logger.log("In each group, aggregate sample genotypes by self.agg_col.")
        geno_agg_df = geno_df.groupby([self.group_col,self.agg_col])[self.sample_s.index].apply(lambda geno_col: geno_col.any() > 0).astype(int)
        geno_agg_df = geno_agg_df.join(agg_cat_n_s)
        return geno_agg_df
    
    
    def assign_variants_to_pop_frq_cats(self, geno_df, pop_frq_col_l, pop_frq_cat_dict):
        '''Assign variants to allele population frequency categories.
        
        Args:
            geno_df (DataFrame): contains variant identified and population frequencies.
            pop_frq_col_l (list of str): contains the names of the population allele frequency columns in descending order of preference.
            pop_frq_cat_dict (dict of (str,float)): mapping of frequency category name to exclusive upper bound. 
        
        Returns:
            geno_df (DataFrame): the same DataFrame which is given as input with an extra column called "pop_frq_cat".
        '''
        
        self.logger.log("Assign variants to a population frequency category.")
        pop_frq_cat_dict = OrderedDict(sorted(pop_frq_cat_dict.items(),key=operator.itemgetter(1)))
        pop_frq_cat_l = list(pop_frq_cat_dict.keys())
        pop_frq_bin_arr = np.array(list(pop_frq_cat_dict.values()))
        pop_frq_s = geno_df.apply(lambda row_s: list(filter(lambda x: np.isnan(x) == False, row_s.ix[pop_frq_col_l].tolist()+[0.0]))[0],axis=1)
        pop_frq_idx_arr = np.digitize(pop_frq_s.values,pop_frq_bin_arr)
        geno_df = geno_df.join(pd.Series(data=pop_frq_idx_arr, index=pop_frq_s.index, name="pop_frq_cat_idx"))
        geno_df["pop_frq_cat"] = geno_df.apply(lambda row_s: row_s.name if row_s["pop_frq_cat_idx"] == pop_frq_bin_arr.size else pop_frq_cat_l[row_s["pop_frq_cat_idx"]], axis=1)
        geno_df.drop("pop_frq_cat_idx", inplace=True, axis=1)
        return geno_df
        
    
    def do_multivariate_test(self, geno_agg_gene_df, y, covar_df=None):
        '''Do a multivariate test for 1 gene.
        
        Args:
            geno_agg_gene_df (DataFrame): contains the genotypes for 1 gene aggregated in by the aggregation column e.g. by population allele frequency.
            y (numpy.ndarray): affection status where 1=affected and 0=unaffected. 
            covar_df (DataFrame): contains the covariates. 
            
        Returns:
            test_result_s (Series): contains the multivariate test results.
        '''
        
        geno_agg_gene_df.index = geno_agg_gene_df.index.droplevel() #Drop the group name so it won't be in the agg cat variable names.
        logit_result = self.fit_logit_model(geno_agg_gene_df,y)
        test_result_l = [("llr_p",logit_result.llr_pvalue)]
        if covar_df is not None:
            '''Model with only covariates'''
            logit_result_h0 = self.fit_logit_model(covar_df,y)
            '''Model with covariates plus aggregated variant independent variables.'''
            logit_result_h1 = self.fit_logit_model(pd.concat([covar_df, geno_agg_gene_df]),y)
            llr_cov_p = stats.chisqprob(2*(logit_result_h1.llf - logit_result_h0.llf), logit_result_h1.df_model - logit_result_h0.df_model)
            test_result_l.append(("llr_cov_p",llr_cov_p))
            h1_coef_l = [(agg_cat+"_c",logit_result_h1.params.loc[agg_cat]) if agg_cat in logit_result_h1.params.index else (agg_cat+"_c",np.NaN) for agg_cat in self.agg_val_l]
            h1_pval_l = [(agg_cat+"_p",logit_result_h1.pvalues.loc[agg_cat]) if agg_cat in logit_result_h1.pvalues.index else (agg_cat+"_p",np.NaN) for agg_cat in self.agg_val_l]
            h1_coef_pval_l = [item for sublist in zip(h1_coef_l,h1_pval_l) for item in sublist]
            test_result_l.extend(h1_coef_pval_l)
        test_result_s = pd.Series(OrderedDict(test_result_l))
        return test_result_s
    
    
    def fit_logit_model(self, X_df, y):
        '''Fit a logit model.
        
        Args:
            X_df: the independent variables (covariates and or aggregated genotypes) for 1 gene. 
            y (numpy.ndarray): sample affection status where 1=affected and 0=unaffected. 
        Returns:
            logit_result (statsmodels.discrete.discrete_model.BinaryResultsWrapper): contains results from fitting logit regression model.
        '''
        
        logit_model = sm.Logit(y,X_df.transpose())
        logit_result = logit_model.fit(method='bfgs')
        return logit_result
    
    
    def get_agg_cat_count_df(self, geno_agg_df):
        '''For each gene, get the number of variants in each population frequency category.
        
        Args:
            geno_agg_df (DataFrame): contains the genotypes aggregated by the aggregation column in each gene.
            
        Returns:
            agg_cat_count_df (DataFrame): contains the number of variants in each aggregtion category in each group.
        '''
        
        agg_cat_count_df = geno_agg_df['n'].reset_index()
        agg_cat_count_df[self.agg_col] = agg_cat_count_df[self.agg_col].apply(lambda x: x if x in self.agg_val_l else "unagg")
        agg_cat_count_df = pd.pivot_table(data=agg_cat_count_df, values='n', index=self.group_col, columns=self.agg_col)
        agg_cat_count_df = agg_cat_count_df.reindex(columns=self.agg_val_l+["unagg"])
        agg_cat_count_df.fillna(0, inplace=True)
        agg_cat_count_df = agg_cat_count_df.apply(pd.to_numeric, downcast='integer', axis=1)
        return agg_cat_count_df
        