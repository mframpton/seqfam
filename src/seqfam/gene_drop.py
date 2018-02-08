import numpy as np
import pandas as pd
from collections import OrderedDict
import time
from seqfam.misc import Logger
import sys
from natsort import natsorted
from itertools import groupby


class Node(object):
    '''Represents an individual (node) in a family tree.'''

    def __init__(self, id, parent_l, spouse, children_l):
        self.id = id
        self.parent_l = parent_l
        self.spouse = spouse
        self.children_l = children_l
        self.genotype = None
    
    
    def get_summary_str(self):
        '''Get a summary string for the node.
        
        Returns:
            summary_str: the object attributes as a string.'''
        
        summary_str_l = ["ID: {0}".format(self.id),
                         "Parents: {0}".format(",".join([parent_node.id for parent_node in self.parent_l])),
                         "Spouse: {0}".format(self.spouse.id if self.spouse != None else ""),
                         "Children: {0}".format(",".join([child_node.id for child_node in self.children_l])),
                         "Genotype: {0}".format(self.genotype)]
        summary_str = "; ".join(summary_str_l)
        return summary_str


class NodeGenerator(object):
    '''Generates nodes to represent individuals in a pedigree.'''

    def convert_ped_df_to_node_l(self, ped_df):
        '''Create the nodes.
        
        Args:
            ped_df (DataFrame): 
            
        Returns:
            node_l (list of Nodes): the Node objects represent the individuals in the pedigree DataFrame.'''

        #Create the nodes.
        name_node_dict = OrderedDict(ped_df.apply(lambda row_s: (row_s["PERSON"], Node(row_s["PERSON"],[],None,[])), axis=1).tolist())
        #Now set the relationships between them.
        ped_df.apply(self.set_relationships, name_node_dict=name_node_dict, axis=1)
        node_l = name_node_dict.values()
        return node_l
        

    def set_relationships(self, ped_row_s, name_node_dict):
        '''Set the relationships between the nodes.
        
        Args:
            ped_row_s (Series): row from a pedigree dataframe, representing 1 individual.
            name_node_dict (dictionary): maps from node name to node.'''
        
        [person,father,mother] = ped_row_s[["PERSON","FATHER","MOTHER"]].tolist()
        if father != "0":
            name_node_dict[person].parent_l = [name_node_dict[father],name_node_dict[mother]]
            name_node_dict[father].children_l.append(name_node_dict[person])
            name_node_dict[mother].children_l.append(name_node_dict[person])
            if name_node_dict[father].spouse == None:
                name_node_dict[father].spouse = name_node_dict[mother]
            if name_node_dict[mother].spouse == None:
                name_node_dict[mother].spouse = name_node_dict[father]
            

class FamilyTree(object):
    '''Represents a family tree with node objects and has methods to perform gene-dropping.'''

    def __init__(self, logger, id, node_l):
        self.logger = logger
        self.id = id
        self.node_l = node_l
        self.founder_l = list(filter(lambda x: not x.parent_l, self.node_l))


    def gene_drop(self, pop_af, ind_gtyped_l):
        '''Main method for performing gene dropping for 1 variant in this family tree.
        
        Args:
            pop_af (float): population allele frequency of the variant.
            ind_gtyped_l (list of strs): IDs of individuals in this family who have a genotype for the variant. 
        Returns:
            carrier_allele_count (int): the number of (genotyped) carriers in the family after dropping the gene.'''
        
        #Set founder genotypes using the population allele frq & select 1 per spousal pair from which to start gene dropping.
        founder_no_spouse_l = []
        for founder in self.founder_l:
            founder.genotype = np.random.binomial(2,pop_af)
            if founder.spouse not in founder_no_spouse_l:
                founder_no_spouse_l.append(founder)
        #Perform gene dropping from each of the selected founders (using a depth-first search). 
        for founder_node in founder_no_spouse_l:
            self.gene_drop_dfs(founder_node)
        #Get the carrier allele count and reset the genotypes.
        carrier_allele_count = sum((node.genotype for node in self.node_l if node.id in ind_gtyped_l))
        for node in self.node_l:
            node.genotype = None
        return carrier_allele_count


    def gene_drop_dfs(self, start_node):
        '''Perform gene dropping in the family tree starting from the specified start node and using a depth-first traversal.
        
        Args:
            start_node (Node): node from which to start the gene dropping.
            
        Returns:
            visited (set): nodes visited during the depth-first traversal.'''
        
        visited, stack = set(), [start_node]
        while stack:
            node = stack.pop(0)
            if node not in visited:
                #print "Visiting " + node.id
                if node.genotype == None:
                    self.set_offspring_genotype(node)
                #print node.genotype
                visited.add(node)
                # new nodes are added to the start of stack
                #if node.spouse != None:
                    #stack = [node.spouse] + stack
                if len(node.children_l) > 0:
                    stack = node.children_l + stack
        #print ""
        return visited
    

    def set_offspring_genotype(self, node):
        '''Set a node's genoytpe based on its parents' genotypes. If the parent is heterozygous then the probability of the variant allele being passed to the offspring is 0.5.
        
        Args:
            node (Node): the node whose genotype will be set.'''
        
        get_allele_passed_to_offspring = lambda parent_genotype: 0 if parent_genotype == 0 else 1 if parent_genotype == 2 else np.random.randint(0,2) if parent_genotype == 1 else None 
        if node.parent_l[0].genotype == None or node.parent_l[1].genotype == None:
            node.genotype = None
        else:
            node.genotype = get_allele_passed_to_offspring(node.parent_l[0].genotype) + get_allele_passed_to_offspring(node.parent_l[1].genotype)


    def log_all_genotypes(self):
        '''Log all of the node genotypes.'''

        self.logger.log(" ".join(["FAMILY {}:".format(self.id)] + ["{0}:{1}".format(node.id,node.genotype) for node in self.node_l]))


    def log_all_info(self):
        '''Log all of the information about each node in the family tree.'''
        
        summary_str_l = ["FAMILY: {0}".format(self.id)]
        summary_str_l += [node.get_summary_str() + "; Founder: {0}".format(node in self.founder_l) for node in self.node_l]
        self.logger.log("\n".join(summary_str_l))
            

class Cohort(object):
    '''Represents a cohort of familial individuals as a list of FamilyTree objects.'''
    
    def __init__(self, cohort_tsv):
        self.logger = Logger()
        self.node_generator = NodeGenerator()
        cohort_df = pd.read_csv(cohort_tsv, sep="\t", dtype=str)
        self.logger.log("Making family tree objects...")
        self.fam_tree_l = cohort_df.groupby("FAMILY").apply(self.make_fam_tree).tolist()
        self.fam_tree_l = natsorted(self.fam_tree_l, key=lambda fam_tree: fam_tree.id)
    
    
    def make_fam_tree(self, ped_df):
        '''For a family, make a list of Nodes and from these, a Family Tree.
        
        Args:
            ped_df (pandas.core.frame.DataFrame): contains the pedigree information for the family.
        
        Returns:
            family_tree (FamilyTree obj): represents the family.'''
        
        family_id = ped_df.iloc[0]["FAMILY"]
        node_l = self.node_generator.convert_ped_df_to_node_l(ped_df)
        return FamilyTree(self.logger, family_id, node_l)
    
    
    def get_all_sample_l(self):
        '''Get a list of all of the samples in the cohort.
        
        Returns:
            all_sample_l (list): list of sample IDs of the form <FAMILY_ID>_<INDIVIDUAL_ID>.'''
        
        all_sample_l = ["{0}_{1}".format(fam_tree.id,node.id) for fam_tree in self.fam_tree_l for node in fam_tree.node_l]
        return all_sample_l
    
    
    def gene_drop(self, pop_af, cohort_af, sample_genotyped_l, gene_drop_n):
        '''Perform gene dropping across the cohort and return the proportion of iterations in which the simulated allele frequency is less than or equal to the cohort frequency.
        
        Args:
            pop_af (float): population allele frequency.
            cohort_af (float): cohort allele frequency.
            sample_genotyped (list of strs): the list of samples genotyped for this variant from which cohort af was calculated.
            gene_drop_n (int): number of iterations to perform.
        
        Returns:
            cohort_enriched_p (float): proportion of iterations in which the simulated allele frequency is less than or equal to the cohort frequency.'''
        
        ind_gtyped_grps = groupby(sample_genotyped_l, lambda x: x.split("_")[0])
        fam_ind_gtyped_dict = OrderedDict([(fam,[id.split("_")[1] for id in list(ind_gtyped_grp)]) for fam, ind_gtyped_grp in ind_gtyped_grps])
        self.logger.log("# families with >=1 genotyped sample: {0}".format(len(fam_ind_gtyped_dict)))
        self.logger.log("Start gene drop with pop_af={0}, cohort_af={1}, # genotype calls={2} & gene_drop_n={3}.".format(pop_af,cohort_af,len(sample_genotyped_l),gene_drop_n))
        
        t0 = time.time()
        total_allele_count = len(sample_genotyped_l)*2
        lt_thresh_n = 0
        for n in range(gene_drop_n):
            carrier_allele_count = sum(fam_tree.gene_drop(pop_af,fam_ind_gtyped_dict[fam_tree.id]) for fam_tree in self.fam_tree_l if fam_tree.id in fam_ind_gtyped_dict)
            gene_drop_af = carrier_allele_count/float(total_allele_count)
            if cohort_af <= gene_drop_af:
                lt_thresh_n += 1        
        cohort_enriched_p = lt_thresh_n/float(gene_drop_n)
        self.logger.log("p-value={0}".format(cohort_enriched_p))
        t1 = time.time()
        self.logger.log("Processing time: {0:.2f} secs\n".format(t1-t0))
        return cohort_enriched_p
