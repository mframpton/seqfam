import pandas as pd
import numpy as np
import statsmodels.discrete.discrete_model as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
import sys
from seqfam.misc import Logger
import operator
from collections import OrderedDict


class CMC(object):
    '''Implements the combined multivariate and collapsing (CMC) burden test, where the multivariate test is a log-likelihood ratio test.'''
    
    def __init__(self):
    
        self.logger = Logger()
        self.sample_s = None
        self.geno_df = None
        self.group_col = None
        self.agg_col = None
        self.agg_cat_l = None
    
    
    def do_multivariate_tests(self, sample_s, geno_df, group_col, agg_col, agg_cat_l, covar_df=None, results_path=None):
        '''Main method for doing multivariate tests.
        
        Args:
            | sample_s (Series): index is the sample names, and values are 1/2 (unaffected/affected).
            | geno_df (DataFrame): index is the variant ID and columns must include (1) group_col (see below); (2) agg_col (see below); (3) sample genotypes (# copies of alternate allele 0-2).
            | group_col (str): column in geno_df identifying the groups of variants to be tested for association with the phenotype e.g. gene, pathway or any other entity.
            | agg_col (str): column containing the aggregation categories (e.g. population allele frequency range) - variants in the same group and aggregation category will be aggregated (collapsed into a dichotomous variable 0/1).
            | agg_cat_l (list of strs): specifies which aggregation categories to present results for (i.e. coefficient and p-value in the alternative hypothesis logit model).
            | covar_df (DataFrame): index is the covariate names and columns are the samples.
            | results_path (str): path to results file.
        
        Returns:
            result_df (DataFrame): multivariate test results, where the index is the variant group (e.g. gene) and columns are (1) # variants in each aggregated (& unaggregated) category; (2) proportion of (un)affecteds carrying a variant in each aggregated (& unaggregated) category; (3) llr_p (and llr_cov_p), the log-likelihood ratio p-value (after controlling for covariates); (4) coefficient & p-value for each independent variable.
        '''
        
        #Set object attributes.
        self.group_col = group_col
        self.agg_col = agg_col
        self.agg_cat_l = agg_cat_l
        #Read in samples, genotypes and covariates.
        self.logger.log("Reading in samples and annotated genotypes...")
        self.sample_s = sample_s
        self.sample_s["Affection"] = self.sample_s["Affection"].astype(int)
        self.sample_s = self.sample_s[self.sample_s != 0]
        self.geno_df = geno_df
        geno_df[self.sample_s.index] = geno_df[self.sample_s.index].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
        geno_df.loc[:,self.sample_s.index].fillna(0, inplace=True)
        self.covar_df = covar_df
        self.logger.log("# variants: {0}".format(len(geno_df.index)))
        self.logger.log("# genes to test: {0}".format(len(pd.unique(geno_df[self.group_col]))))
        #Aggregate the genotypes.
        self.logger.log("Aggregating genotypes...")
        geno_agg_df = self.aggregate_by_agg_col(geno_df)
        self.logger.log("Retain groups which have >= 1 non-zero genotype.")
        all_zero_s = geno_agg_df[self.sample_s.index].groupby(level=group_col).apply(lambda df: df.values.sum()==0)
        if all_zero_s.size > 0:
            self.logger.log("Dropped because they have all zero genotypes: {0}".format(",".join(all_zero_s[all_zero_s==True].index.tolist())))
            geno_agg_df = geno_agg_df.ix[all_zero_s[all_zero_s==False].index,:]
        #Do the multivariate tests.
        self.logger.log("Doing multivariate tests...")
        result_df = geno_agg_df[self.sample_s.index].groupby(level=[self.group_col]).apply(self.do_multivariate_test, y=self.sample_s.values-1, covar_df=covar_df)
        #Merge the results with the collapsed variant counts.
        self.logger.log("Merge the results with the population frequency variant category counts...")
        agg_cat_count_df = self.get_agg_cat_count_df(geno_agg_df)
        agg_cat_count_df = agg_cat_count_df.join(self.get_agg_cat_prop_by_affection(geno_agg_df))
        result_df = agg_cat_count_df.join(result_df)
        if covar_df is None:
            result_df.sort_values(by="llr_p", inplace=True)
        else:
            result_df.sort_values(by="llr_cov_p", inplace=True)
        self.logger.log("Write results.")
        result_df.to_csv(results_path, index=True)
        self.logger.log("\n"+result_df.to_string())
        return result_df
    
    
    def aggregate_by_agg_col(self, geno_df):
        '''Aggregate genotypes within variant population frequency categories.
        
        Args:
            geno_df (DataFrame): index is the variant ID and columns must include (1) group_col (see below); (2) agg_col (see below); (3) sample genotypes (# copies of alternate allele 0-2).
            
        Returns:
            geno_agg_df (DataFrame): index is the gene & variant aggregation category, and columns are the # of variants for each sample.
        '''

        self.logger.log("For each gene, count the # of variants in each aggregation category.")
        agg_cat_n_s = geno_df.groupby(self.group_col)[self.agg_col].value_counts()
        agg_cat_n_s.name = "n"
        self.logger.log("In each group, aggregate sample genotypes by self.agg_col.")
        geno_agg_df = geno_df.groupby([self.group_col,self.agg_col])[self.sample_s.index].apply(lambda geno_col: geno_col.any() > 0).astype(int)
        geno_agg_df = geno_agg_df.join(agg_cat_n_s)
        return geno_agg_df
    
    
    def assign_variants_to_pop_frq_cats(self, geno_df, pop_frq_col_l, pop_frq_cat_dict):
        '''Assign variants to allele population frequency range categories.
        
        Args:
            | geno_df (DataFrame): index is the variant ID and columns must include the population allele frequency columns listed in the pop_frq_col_l parameter (see below)
            | pop_frq_col_l (list of str): contains the names of the population allele frequency columns in descending order of preference.
            | pop_frq_cat_dict (dict of (str,float)): mapping of frequency category name to exclusive upper bound. 
        
        Returns:
            geno_df (DataFrame): the inputted geno_df DataFrame with an extra column for the variant aggregation category ("pop_frq_cat").
        '''
        
        self.logger.log("Assign variants to a population frequency category.")
        pop_frq_cat_dict = OrderedDict(sorted(pop_frq_cat_dict.items(),key=operator.itemgetter(1)))
        pop_frq_cat_l = list(pop_frq_cat_dict.keys())
        pop_frq_bin_arr = np.array(list(pop_frq_cat_dict.values()))
        pop_frq_s = geno_df.apply(lambda row_s: list(filter(lambda x: pd.notnull(x), row_s.loc[pop_frq_col_l].tolist()+[0.0]))[0],axis=1)
        pop_frq_idx_arr = np.digitize(pop_frq_s.values,pop_frq_bin_arr)
        geno_df = geno_df.join(pd.Series(data=pop_frq_idx_arr, index=pop_frq_s.index, name="pop_frq_cat_idx"))
        geno_df["pop_frq_cat"] = geno_df.apply(lambda row_s: row_s.name if row_s["pop_frq_cat_idx"] == pop_frq_bin_arr.size else pop_frq_cat_l[row_s["pop_frq_cat_idx"]], axis=1)
        geno_df.drop("pop_frq_cat_idx", inplace=True, axis=1)
        return geno_df
        
    
    def do_multivariate_test(self, geno_agg_gene_df, y, covar_df=None):
        '''Do a multivariate test for 1 gene.
        
        Args:
            | geno_agg_gene_df (DataFrame): index is the gene & variant aggregation category, and columns are the # of variants for each sample.
            | y (numpy.ndarray): values are 0/1 (unaffected/affected). 
            | covar_df (DataFrame): index is the covariate names and columns are samples.
            
        Returns:
            test_result_s (Series): multivariate test results for 1 gene/functional unit containing (1) llr_p (and llr_cov_p), the log-likelihood ratio p-value (after controlling for covariates); (2) coefficient & p-value for each independent variable.
        '''
        
        geno_agg_gene_df.index = geno_agg_gene_df.index.droplevel() #Drop the group name so it won't be in the agg cat variable names.
        logit_result = self.fit_logit_model(geno_agg_gene_df,y)
        test_result_l = [("llr_p",logit_result.llr_pvalue)]
        coef_pval_l = None
        if covar_df is not None:
            #Model with only covariates
            logit_result_h0 = self.fit_logit_model(covar_df,y)
            #Model with covariates plus aggregated variant independent variables
            logit_result_h1 = self.fit_logit_model(pd.concat([covar_df, geno_agg_gene_df]),y)
            llr_cov_p = stats.chisqprob(2*(logit_result_h1.llf - logit_result_h0.llf), logit_result_h1.df_model - logit_result_h0.df_model)
            test_result_l.append(("llr_cov_p",llr_cov_p))
            coef_pval_l = self.get_coef_pval_l(logit_result_h1, covar_b=True)
        else:
            coef_pval_l = self.get_coef_pval_l(logit_result, covar_b=False)
        test_result_l.extend(coef_pval_l)
        test_result_s = pd.Series(OrderedDict(test_result_l))
        return test_result_s
    
    
    def get_coef_pval_l(self, logit_result, covar_b=False):
        '''Get the coefficients and corresponding p-values of the independent variables in the logit model.
        
        Args:
            | logit_result (statsmodels.discrete.discrete_model.BinaryResultsWrapper): contains results from fitting logit regression model.
        
        Returns:
            coef_pval_l (list): list of coefficients and corresponding p-values.
        '''
        
        ind_var_l = self.agg_cat_l if covar_b == False else self.agg_cat_l + self.covar_df.index.tolist()
        coef_l = [(agg_cat+"_c",logit_result.params.loc[agg_cat]) if agg_cat in logit_result.params.index else (agg_cat+"_c",np.NaN) for agg_cat in ind_var_l]
        pval_l = [(agg_cat+"_p",logit_result.pvalues.loc[agg_cat]) if agg_cat in logit_result.pvalues.index else (agg_cat+"_p",np.NaN) for agg_cat in ind_var_l]
        coef_pval_l = [item for sublist in zip(coef_l,pval_l) for item in sublist]
        return coef_pval_l
    
    
    def fit_logit_model(self, X_df, y):
        '''Fit a logit model.
        
        Args:
            | X_df: the independent variables (covariates and or aggregated genotypes) for 1 gene. 
            | y (numpy.ndarray): values are 0/1 (unaffected/affected). 
        
        Returns:
            logit_result (statsmodels.discrete.discrete_model.BinaryResultsWrapper): contains results from fitting logit regression model.
        '''
        
        logit_model = sm.Logit(y,X_df.transpose())
        logit_result = logit_model.fit(method='bfgs')
        return logit_result
    
    
    def get_agg_cat_count_df(self, geno_agg_df):
        '''For each group in group_col (e.g. gene), get the number of variants in each variant aggregation category (in agg_col) e.g. population allele frequency range.
        
        Args:
            geno_agg_df (DataFrame): index is the gene & variant aggregation category, and columns are the # of variants for each sample.
            
        Returns:
            agg_cat_count_df (DataFrame): index is the gene and columns are # variants in each aggregated category and # unaggregated.
        '''
        
        agg_cat_count_df = geno_agg_df['n'].reset_index()
        agg_cat_count_df[self.agg_col] = agg_cat_count_df[self.agg_col].apply(lambda x: x if x in self.agg_cat_l else "unagg")
        
        def sum_unagg_n(agg_cat_count_df):
            '''If the user has not listed all of the aggregation categories in the agg_cat_l parameter of the do_multivariate_tests method,
            then agg_cat_count_df may contain multiple rows per gene for "unagg", and so these must be summed.
            
            Args:
                agg_cat_count_df (DataFrame): contains the group_col, agg_col and a count (n) column.
                
            Returns:
                agg_cat_count_df (DataFrame): contains the group_col, agg_col and a count (n) column.
            '''
            
            unagg_count_df = agg_cat_count_df.loc[agg_cat_count_df["pop_frq_cat"]=="unagg",:].groupby("Gene")["n"].sum().to_frame()
            unagg_count_df.reset_index(inplace=True)
            unagg_count_df["pop_frq_cat"] = "unagg"
            agg_cat_count_df = pd.concat([agg_cat_count_df.loc[agg_cat_count_df["pop_frq_cat"]!="unagg",:],unagg_count_df])
            return agg_cat_count_df
        
        agg_cat_count_df = sum_unagg_n(agg_cat_count_df)
        agg_cat_count_df = pd.pivot_table(data=agg_cat_count_df, values='n', index=self.group_col, columns=self.agg_col)
        agg_cat_count_df = agg_cat_count_df.reindex(columns=self.agg_cat_l+["unagg"])
        agg_cat_count_df.fillna(0, inplace=True)
        agg_cat_count_df = agg_cat_count_df.apply(pd.to_numeric, downcast='integer', axis=1)
        return agg_cat_count_df
    
    
    def get_agg_cat_prop_by_affection(self, geno_agg_df):
        '''For each group in group_col (e.g. gene), get the proportion of (un)affecteds who are carriers in each variant aggregation category (in agg_col) e.g. population allele frequency range.
              
        Args:
            geno_agg_df (DataFrame): index is the group_col & aggregation category, and columns are the # of variants for each sample.
        
        Returns:
            agg_cat_prop_df(DataFrame): index is the group_col, and columns indicate the proportion of (un)affected carriers in each variant aggregation category, plus unaggregated variants.
        '''
        
        geno_agg_df = geno_agg_df.reset_index()
        geno_agg_df[self.agg_col] = geno_agg_df[self.agg_col].apply(lambda x: x if x in self.agg_cat_l else "unagg")
        
        def get_prop_s(geno_agg_df, sample_l, name):
            '''
            Args:
                | geno_agg_df (DataFrame): index is the group_col & aggregation category, and columns are the # of variants for each sample.
                | sample_l (list of strs): list of sample names (columns in geno_agg_df).
                | name (str): name to give to returned Series. 
            
            Returns:
                prop_s (Series): index is the group_col & aggregation category, and values are the proportion of (un)affecteds who are carriers.
            '''
            
            prop_s = geno_agg_df.groupby([self.group_col,self.agg_col]).apply(lambda group: group[sample_l].values.mean())
            prop_s.name = name
            return prop_s
        
        ca_l,co_l = self.sample_s[self.sample_s["Affection"]==2].index.tolist(), self.sample_s[self.sample_s["Affection"]==1].index.tolist()
        prop_df = pd.concat([get_prop_s(geno_agg_df,ca_l,"aff_p"),get_prop_s(geno_agg_df,co_l,"unaff_p")], axis=1)
        prop_df.reset_index(inplace=True)

        def pivot_prop_df(prop_df, affection):
            '''
            Args:
                | prop_df (DataFrame): contains the group_col and agg_col columns, plus columns for the proportion of (un)affected carriers. 
                | affection (str): affection 
            
            Returns:
                agg_cat_prop_df (DataFrame): index is the group_col, and columns are the proportion of (un)affecteds who are carriers of each variant aggregation category in agg_cat_l, plus unaggregated variants.
            '''
            
            agg_cat_prop_df = pd.pivot_table(data=prop_df, values="{0}_p".format(affection), index=self.group_col, columns=self.agg_col)
            column_rename_dict = dict(zip(agg_cat_prop_df.columns.tolist(),["{0}_{1}".format(col,"{0}_p".format(affection)) for col in agg_cat_prop_df.columns.tolist()]))
            agg_cat_prop_df.rename(columns=column_rename_dict, inplace=True)
            agg_cat_prop_df = agg_cat_prop_df.reindex(columns=["{0}_{1}_p".format(agg_cat,affection) for agg_cat in self.agg_cat_l + ["unagg"]])
            return agg_cat_prop_df

        agg_cat_prop_df = (pivot_prop_df(prop_df,"aff")).merge(pivot_prop_df(prop_df,"unaff"), left_index=True, right_index=True)
        agg_cat_prop_df = agg_cat_prop_df.reindex(columns=["{0}_{1}_p".format(agg_cat,affection) for agg_cat in self.agg_cat_l + ["unagg"] for affection in ["aff","unaff"]])
        
        return agg_cat_prop_df 
    