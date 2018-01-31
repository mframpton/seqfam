import sys
import pandas as pd
import numpy as np
import itertools
from collections import OrderedDict
from cohpy.misc import Logger


class Relatedness(object):
    '''Provides methods to find duplicates either between or within families, and to identify pairs of individuals within a family whose observed
    relationship (KING kinship coefficient) is different than expected given the pedigree.'''
    
    def __init__(self,bf_file=None,wf_file=None,cohort_tsv=None,kinship_coef_thresh_dict={"0":0.354,"1":0.177,"2":0.0844,"3":0.0442}):
        
        self.logger = Logger()
        self.bf_file = bf_file
        self.wf_file = wf_file
        self.kinship_coef_thresh_dict = kinship_coef_thresh_dict
        self.cohort_tsv = cohort_tsv
        self.bf_df = None
        self.wf_df = None
        self.cohort_df = None

        
    def find_duplicates(self, bf_b=True):
        '''Find either between-family or within-family duplicates using the corresponding KING kinship coefficient output.
        
        Args:
            bf_b (boolean): whether to look for duplicates from different families (True) or duplicates within the same family, defaults to True.
        
        Returns:
            duplicate_l (list of strs): pairs of identified duplicates.'''
        
        self.logger.log("Checking for {between_or_within} duplicates...".format(between_or_within="between-family" if bf_b == True else "within-family"))
        kinship_coef_file = self.bf_file if bf_b == True else self.wf_file
        index_col = ["FID1","ID1","FID2","ID2"] if bf_b == True else ["FID","ID1","ID2"]
        kinship_coef_df = pd.read_csv(kinship_coef_file, sep="\t", index_col=index_col, dtype=str)
        kinship_coef_df["Kinship"] = pd.to_numeric(kinship_coef_df["Kinship"])
        duplicate_df = kinship_coef_df.ix[kinship_coef_df["Kinship"] >= self.kinship_coef_thresh_dict["0"],:]
        duplicate_l = []
        if ~duplicate_df.empty:
            duplicate_l = [self.convert_tuple_to_ids(duplicate_pair_tuple) for duplicate_pair_tuple in duplicate_df.index.tolist()]
        return duplicate_l
        
        
    def convert_tuple_to_ids(self, duplicate_pair_tuple):
        '''Convert a tuple containing the IDs of 2 duplicate individuals into a string. 
        
        Args:
            duplicate_pair_tuple (tuple): tuple containing IDs of 2 duplicate individuals.
        
        Returns:
            duplicate_pair_str (str): string containing IDs of the 2 duplicate individuals.'''
        
        #print duplicate_pair_tuple
        id1,id2 = "",""
        if len(duplicate_pair_tuple) == 3:
            id1,id2 = "_".join([duplicate_pair_tuple[0],duplicate_pair_tuple[1]]),"_".join([duplicate_pair_tuple[0],duplicate_pair_tuple[2]])
        elif len(duplicate_pair_tuple) == 4:
            id1,id2 = "_".join([duplicate_pair_tuple[0],duplicate_pair_tuple[1]]),"_".join([duplicate_pair_tuple[2],duplicate_pair_tuple[3]])
        duplicate_pair_str = ",".join([id1,id2])
        return duplicate_pair_str


    def get_exp_obs_df(self):
        '''Make a DataFrame containing the expected and observed relationships for individuals in the same family.
        
        Returns:
            exp_obs_df (DataFrame): contains columns for the expected and observed relationships.'''
        
        self.logger.log("Checking expected versus observed within-family relationships...")
        if self.wf_df == None:
            self.wf_df = pd.read_csv(self.wf_file, sep="\t", dtype=str)
        if self.cohort_df == None:
            self.cohort_df = pd.read_csv(self.cohort_tsv, sep="\t", dtype=str)

        self.logger.log("Cohort contains {0} individuals and {1} families.".format(len(self.cohort_df.drop_duplicates(subset=["PERSON","FAMILY"]).index),
                                                                                   pd.unique(self.cohort_df["FAMILY"]).size))
        fam_ind_seq_dict = self.wf_df.groupby("FID").apply(lambda fam_df: pd.unique(fam_df["ID1"].append(fam_df["ID2"])).tolist()).to_dict()
        self.logger.log("# individuals with >=1 kinship coefficient: {0}".format(sum([len(fam_ind_seq_dict[fam]) for fam in fam_ind_seq_dict])))
        self.logger.log("Removing families with no kinship coefficients.") 
        self.cohort_df = self.cohort_df.ix[self.cohort_df["FAMILY"].isin(fam_ind_seq_dict),:]
        self.logger.log("Cohort contains {0} individuals and {1} families.".format(len(self.cohort_df.drop_duplicates(subset=["PERSON","FAMILY"]).index),
                                                                                   pd.unique(self.cohort_df["FAMILY"]).size))
        self.logger.log("Get the expected relationships between members of each family.")
        exp_obs_df = self.cohort_df.groupby(by="FAMILY", sort=False).apply(self.get_exp_rel_df, fam_ind_seq_dict)
        self.logger.log("Obtained expected relationship for {0} pairs across {1} families.".format(len(exp_obs_df.index),
                                                                                                   len(exp_obs_df.index.get_level_values(0).unique())))
        
        self.logger.log("Get and merge the observed relationships between members of each family.")
        exp_obs_df["Kinship"] = exp_obs_df.apply(func=self.get_kinship_coef, axis=1)
        self.logger.log("Retrieved kinship coefficients for {0} pairs across {1} families.".format(exp_obs_df["Kinship"].notnull().sum(),
                                                                                                   len(exp_obs_df.ix[exp_obs_df["Kinship"].notnull(),:].index.get_level_values(0).unique())))
        exp_obs_df["OBS_REL"] = exp_obs_df["Kinship"].apply(lambda x: np.NaN if pd.isnull(x) else "0" if x >= self.kinship_coef_thresh_dict["0"] else "1" if x >= self.kinship_coef_thresh_dict["1"] else "2" if x >= self.kinship_coef_thresh_dict["2"] else "3" if x >= self.kinship_coef_thresh_dict["3"] else "4")
        exp_obs_df.reset_index(inplace=True)
        exp_obs_df.drop("level_1", inplace=True, axis=1)
        exp_obs_df.set_index(["FAMILY","ID1","ID2"], inplace=True)
        return exp_obs_df


    def get_exp_rel_df(self, fam_df, fam_ind_seq_dict):
        '''Make a dataframe containing the expected relationships between individuals in a family.
        
        Args:
            fam_df (DataFrame): contains the pedigree information for a family.
        
        Returns:
            exp_rel_df (DataFrame): contains the expected relationships between the family members.'''
    
        #print(fam_df.name)        
        fam_df.replace({"0":np.NaN}, inplace=True) #Tidy fam_df
        
        #Get each individual's parents and grandparents.
        fam_uniq_ind_l = pd.unique(fam_df["PERSON"]).tolist()
        fam_uniq_ind_l = list(filter(lambda x: x in fam_ind_seq_dict[fam_df.name], fam_uniq_ind_l)) #Only need relations for samples with kinship coefs.
        relations_df = pd.DataFrame(data={"IND":fam_uniq_ind_l})
        relations_df[["siblings","parents","grandparents"]] = relations_df["IND"].apply(self.get_relations_df, fam_df=fam_df)
        relations_df.set_index("IND", inplace=True)
    
        #Get every pair of individuals and store their expected and observed relationship.
        ind1_l, ind2_l = [],[]
        for subset in itertools.combinations(fam_uniq_ind_l,2):
            ind1_l.append(subset[0])
            ind2_l.append(subset[1])
        exp_rel_df = pd.DataFrame(OrderedDict([("ID1",ind1_l),("ID2",ind2_l)]))
        exp_rel_df["EXP_REL"] = exp_rel_df.apply(axis=1, func=self.get_exp_rel, relations_df=relations_df)
        return exp_rel_df


    def get_relations_df(self, ind, fam_df):
        '''For 1 individual, get a series containing lists of their siblings, parents and grandparents.
        
        Args:
            ind (str): individual ID 
            fam_df (DataFrame): contains the pedigree information for a family.
        
        Returns:
            relations_s (Series): contains lists of the individual's siblings, parents and grandparents.'''
    
        #print "IND: {0}".format(ind)
        father,mother = fam_df.ix[fam_df["PERSON"]==ind,:].iloc[0]["FATHER"],fam_df.ix[fam_df["PERSON"]==ind,:].iloc[0]["MOTHER"]
        #print "PARENTS: {0},{1}".format(father,mother)
        sibling_l = fam_df.ix[(fam_df["FATHER"]==father) & (fam_df["MOTHER"]==mother) & (fam_df["PERSON"]!=ind),"PERSON"].tolist() if pd.notnull(father) and pd.notnull(mother) else []
        #print "SIBLINGS: {0}".format(",".join(sibling_l))
        parental_gf = fam_df.ix[fam_df["PERSON"]==father,:].iloc[0]["FATHER"] if pd.notnull(father) else np.NaN
        #print "PARENTAL GF: {0}".format(parental_gf)
        parental_gm = fam_df.ix[fam_df["PERSON"]==father,:].iloc[0]["MOTHER"] if pd.notnull(father) else np.NaN
        #print "PARENTAL GM: {0}".format(parental_gm)
        maternal_gf = fam_df.ix[fam_df["PERSON"]==mother,:].iloc[0]["FATHER"] if pd.notnull(mother) else np.NaN
        #print "PARENTAL GF: {0}".format(maternal_gf)
        maternal_gm = fam_df.ix[fam_df["PERSON"]==mother,:].iloc[0]["MOTHER"] if pd.notnull(mother) else np.NaN
        #print "PARENTAL GM: {0}".format(maternal_gm)
        parent_l = list(filter(lambda x: pd.notnull(x), [father,mother]))
        grandparent_l = list(filter(lambda x: pd.notnull(x), [parental_gf,parental_gm,maternal_gf,maternal_gm]))
        
        relations_s = pd.Series([sibling_l,parent_l,grandparent_l],
                         index=["siblings","parents","grandparents"])
        return relations_s


    def get_exp_rel(self, ind_pair_s, relations_df):
        '''Get the expected degree of relationship (1-4) between a pair of individuals.
        
        Args:
            ind_pair_s (Series): the IDs of a pair of individuals.
            
        Returns:
            exp_rel (str): the expected degree of relationship.'''
    
        ind1, ind2 = ind_pair_s["ID1"], ind_pair_s["ID2"]
        ind1_sibling_l, ind2_sibling_l = relations_df.loc[ind1,"siblings"], relations_df.loc[ind2,"siblings"]
        ind1_parent_l, ind2_parent_l = relations_df.loc[ind1,"parents"], relations_df.loc[ind2,"parents"]
        ind1_grandparent_l, ind2_grandparent_l = relations_df.loc[ind1,"grandparents"], relations_df.loc[ind2,"grandparents"]
    
        exp_rel = "4" #i.e. at least a 4th-degree relationship.
        if ind2 in ind1_parent_l or ind1 in ind2_parent_l: #Parent-child
            exp_rel = "1"
        elif len(ind1_parent_l) > 0 and set(ind1_parent_l) == set(ind2_parent_l): #Sibling
            exp_rel = "1"
        elif ind2 in ind1_grandparent_l or ind1 in ind2_grandparent_l: #Grandparent-grandchild 
            exp_rel = "2"
        elif (len(ind1_parent_l) > 1 and set(ind1_parent_l).issubset(set(ind2_grandparent_l))) or (len(ind2_parent_l) > 1 and set(ind2_parent_l).issubset(set(ind1_grandparent_l))): #Uncle/Aunt - nephew/niece
            exp_rel = "2"
        elif len(set(ind1_grandparent_l).intersection(set(ind2_grandparent_l))) > 0: #Cousins
            exp_rel = "3"
        elif len(set(ind1_grandparent_l).intersection(set(ind2_sibling_l))) > 0 or len(set(ind2_grandparent_l).intersection(set(ind1_sibling_l))) > 0: #Great Uncle/Aunt
            exp_rel = "3"
        return exp_rel
    
    
    def get_kinship_coef(self, s):
        '''Retrieve the kinship coefficient from a Series representing a row in a KING output file.
        
        Args:
            s (Series): represents a row in a KING output file.
            
        Returns:
            kinship_coef (str): kinship coefficient in the row.'''

        #print s
        kinship_coef = np.NaN
        df = self.wf_df.ix[(self.wf_df["FID"]==s.name[0]) & (self.wf_df["ID1"]==s.ix["ID1"]) & (self.wf_df["ID2"]==s.ix["ID2"]),:]  
        if df.empty == False:
            #print df
            kinship_coef = float(df.iloc[0]["Kinship"])
        else:
            df = self.wf_df.ix[(self.wf_df["FID"]==s.name[0]) & (self.wf_df["ID2"]==s.ix["ID1"]) & (self.wf_df["ID1"]==s.ix["ID2"]),:]
            if df.empty == False:
                kinship_coef = float(df.iloc[0]["Kinship"])
        return kinship_coef
    