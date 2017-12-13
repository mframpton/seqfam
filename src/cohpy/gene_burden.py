import pandas as pd
import numpy as np
import statsmodels.discrete.discrete_model as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)
import sys
from misc import Logger
import operator
from collections import OrderedDict


class CMC(object):
    
    def __init__(self):
    
        self.logger = Logger()
        self.sample_s = None
        self.variant_col = None
        self.gene_col = None
        self.pop_frq_col_l = None
        self.pop_frq_cat_l = None
        self.pop_frq_bin_arr = None
    
    
    def do_multivariate_tests(self, samples_path, genotypes_path, variant_col="VARIANT_ID", gene_col="gene_transcript",
                              pop_frq_col_l=["BroadAJcontrols_ALT","gnomad_AF_NFE","EXAC_NFE"], pop_frq_cat_dict={"rare":0.01,"mod_rare":0.05},
                              covariates_path=None, results_path=None):
        '''Main method for doing multivariate tests.
        
        Args:
            samples_path (str): path to samples file.
            genotypes_path (str): path to genotypes file.
            variant_col (str): unique variant ID column in samples file.
            gene_col (str): gene ID column.
            pop_frq_col_l (list of strs): ordered list of population frequency columns.
            pop_frq_grp_dict (dict): mapping from variant population frequency categories to their upper bounds.
            covariates_path (str): path to covariates file.
            results_path (str): path to results file. 
        
        Returns:
            result_df (DataFrame): multivariate test results.
        '''
        
        #Set object attributes.
        self.variant_col = variant_col
        self.gene_col = gene_col
        self.pop_frq_col_l = pop_frq_col_l
        pop_frq_cat_dict = OrderedDict(sorted(pop_frq_cat_dict.items(),key=operator.itemgetter(1)))
        self.pop_frq_cat_l = pop_frq_cat_dict.keys()
        self.pop_frq_bin_arr = np.array(pop_frq_cat_dict.values())
        #Read in samples, genotypes and covariates.
        self.logger.log("Reading in samples and annotated genotypes...")
        self.sample_s = pd.read_csv(samples_path, dtype=str, index_col="Sample ID")
        self.sample_s["Affection"] = self.sample_s["Affection"].astype(int)
        self.sample_s = self.sample_s[self.sample_s != 0]
        geno_df = pd.read_csv(genotypes_path, dtype=str, usecols=[self.variant_col,self.gene_col] + self.pop_frq_col_l + self.sample_s.index.tolist(), index_col=self.variant_col)
        geno_df[self.pop_frq_col_l] = geno_df[self.pop_frq_col_l].apply(pd.to_numeric, axis=1)
        geno_df[self.sample_s.index.tolist()] = geno_df[self.sample_s.index.tolist()].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
        geno_df[self.sample_s.index.tolist()].fillna(0, inplace=True)
        covar_df = None if covariates_path == None else pd.read_csv(covariates_path, index_col=0)
        self.logger.log("# variants: {0}".format(len(geno_df.index)))
        self.logger.log("# genes to test: {0}".format(len(pd.unique(geno_df[self.gene_col]))))
        #Aggregate the genotypes.
        self.logger.log("Aggregating genotypes...")
        geno_aggregated_df = self.aggregate_in_var_pop_frq_cats(geno_df=geno_df)
        #Do the multivariate tests.
        self.logger.log("Doing multivariate tests...")
        result_df = geno_aggregated_df[self.sample_s.index.tolist()].groupby(level=[self.gene_col]).apply(self.do_multivariate_test, y=self.sample_s.values-1, covar_df=covar_df)
        #Merge the results with the collapsed variant counts.
        self.logger.log("Merge the results with the population frequency variant category counts...")
        pop_frq_cat_count_df = self.get_pop_frq_cat_count_df(geno_aggregated_df)
        result_df = pop_frq_cat_count_df.join(result_df)
        result_df.sort_values(by="llr_cov_p", inplace=True)
        self.logger.log("Write results.")
        result_df.to_csv(results_path, index=True)
        return result_df
    
    
    def aggregate_in_var_pop_frq_cats(self, geno_df):
        '''Aggregate genotypes within variant population frequency categories.
        
        Args:
            geno_df (DataFrame): contains the sample genotypes, variant identifier, gene identifier and population frequencies.
            
        Returns:
            geno_aggregated_df (DataFrame): contains the genotypes aggregated within variant population frequency categories.
        '''

        self.logger.log("Assign variants to a population frequency group.")
        pop_frq_s = geno_df.apply(lambda row_s: filter(lambda x: np.isnan(x) == False, row_s.ix[self.pop_frq_col_l].tolist()+[0.0])[0],axis=1)
        pop_frq_idx_arr = np.digitize(pop_frq_s.values,self.pop_frq_bin_arr)
        geno_df = geno_df.join(pd.Series(data=pop_frq_idx_arr, index=pop_frq_s.index, name="pop_frq_cat_idx"))
        geno_df["pop_frq_cat"] = geno_df.apply(lambda row_s: row_s.name if row_s["pop_frq_cat_idx"] == len(self.pop_frq_cat_l) else self.pop_frq_cat_l[row_s["pop_frq_cat_idx"]], axis=1)
        geno_df.drop("pop_frq_cat_idx", inplace=True, axis=1)
        self.logger.log("For each gene, count the # of variants in each population frequency category.")
        pop_frq_cat_n_s = geno_df.groupby(self.gene_col)["pop_frq_cat"].value_counts()
        pop_frq_cat_n_s.name = "n"
        self.logger.log("In each gene, aggregate sample genotypes by variant population frequency group.")
        geno_aggregated_df = geno_df.groupby([self.gene_col] + ["pop_frq_cat"])[self.sample_s.index].apply(lambda geno_col: geno_col.any() > 0).astype(int)
        #print geno_aggregated_df[self.sample_s.index.tolist()[:2]].to_string()
        #sys.exit()
        geno_aggregated_df = geno_aggregated_df.join(pop_frq_cat_n_s)
        return geno_aggregated_df
    
    
    def do_multivariate_test(self, geno_aggregated_gene_df, y, covar_df=None):
        '''Do a multivariate test for 1 gene.
        
        Args:
            geno_aggregated_gene_df (DataFrame): contains the genotypes for 1 gene aggregated in variant population frequency categories.
            y (numpy.ndarray): affection status where 1=affected and 0=unaffected. 
            covar_df (DataFrame): contains the covariates. 
            
        Returns:
            test_result_s (Series): contains the multivariate test results.
        '''
        
        return_data_l,return_index_l = [],[]
        
        [degf,llf,llr_p] = self.fit_logit_model(geno_aggregated_gene_df,y)
        return_data_l.append(llr_p)
        return_index_l.append("llr_p")
        
        if covar_df is not None:
            '''Model with only covariates'''
            [h0_degf,h0_llf,h0_llr_p] = self.fit_logit_model(covar_df,y)
            '''Model with covariates plus aggregated variant independent variables.'''
            [h1_degf,h1_llf,h1_llr_p] = self.fit_logit_model(pd.concat([geno_aggregated_gene_df, covar_df]),y)
            #print h0_llf,h1_llf,h0_degf,h1_degf,2*(h1_llf-h0_llf),h1_degf-h0_degf
            llr_cov_p = stats.chisqprob(2*(h1_llf-h0_llf),h1_degf-h0_degf)
            return_data_l.append(llr_cov_p)
            return_index_l.append("llr_cov_p")
            
        return pd.Series(data=return_data_l, index=return_index_l) 

    
    def fit_logit_model(self, X_df, y):
        '''Fit a logit model.
        
        Args:
            X_df: the independent variables (covariates and or aggregated genotypes) for 1 gene. 
            y (numpy.ndarray): sample affection status where 1=affected and 0=unaffected. 
        Returns:
            result_l (list of int and 2 floats): degrees of freedom, log-likelihood and log-likelihood ratio p-value. 
        '''
        
        X = np.transpose(X_df.values)
        logit_model = sm.Logit(y,X)
        result = logit_model.fit(method='bfgs')
        degf = result.df_model
        llf = result.llf
        llr_p = result.llr_pvalue
        return [degf,llf,llr_p]
    
    
    def get_pop_frq_cat_count_df(self, geno_aggregated_df):
        '''For each gene, get the number of variants in each population frequency category.
        
        Args:
            geno_aggregated_df (DataFrame): contains the genotypes aggregated in variant population frequency categories in each gene.
            
        Returns:
            pop_frq_cat_count_df (DataFrame): contains the number of variants in each population frequency category in each gene.
        '''
        
        pop_frq_cat_count_df = geno_aggregated_df['n'].reset_index()
        pop_frq_cat_count_df["pop_frq_cat"] = pop_frq_cat_count_df["pop_frq_cat"].apply(lambda x: x if x in self.pop_frq_cat_l else "unagg")
        pop_frq_cat_count_df = pd.pivot_table(data=pop_frq_cat_count_df, values='n', index=self.gene_col, columns='pop_frq_cat')
        pop_frq_cat_count_df = pop_frq_cat_count_df.reindex(columns=self.pop_frq_cat_l+["unagg"])
        pop_frq_cat_count_df.fillna(0, inplace=True)
        pop_frq_cat_count_df = pop_frq_cat_count_df.apply(pd.to_numeric, downcast='integer', axis=1)
        return pop_frq_cat_count_df
        