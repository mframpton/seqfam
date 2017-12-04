import numpy as np
import pandas as pd
from collections import OrderedDict
from datetime import datetime


class Node(object):

    def __init__(self, id, parent_l, spouse, children_l, sequenced_b):
        self.id = id
        self.parent_l = parent_l
        self.spouse = spouse
        self.children_l = children_l
        self.genotype = None
        self.sequenced_b = sequenced_b


class NodeGenerator(object):

    def convert_ped_df_to_node_l(self, ped_df):
        '''Create the nodes.'''
        name_node_dict = OrderedDict()
        for i in ped_df.index.tolist():
            person,sequenced = ped_df.ix[i,"PERSON"], ped_df.ix[i,"SEQUENCED"]
            name_node_dict[person] = Node(person,[],None,[],sequenced)
        '''Set the relationships.'''
        for i in ped_df.index.tolist():
            person = ped_df.ix[i,"PERSON"]
            if ped_df.ix[i,"FATHER"] != "0":
                father,mother = ped_df.ix[i,"FATHER"],ped_df.ix[i,"MOTHER"]
                name_node_dict[person].parent_l = [name_node_dict[father],name_node_dict[mother]]
                name_node_dict[father].children_l.append(name_node_dict[person])
                name_node_dict[mother].children_l.append(name_node_dict[person])
                if name_node_dict[father].spouse == None:
                    name_node_dict[father].spouse = name_node_dict[mother]
                if name_node_dict[mother].spouse == None:
                    name_node_dict[mother].spouse = name_node_dict[father]

        return name_node_dict.values()


class FamilyTree(object):

    def __init__(self, node_l):
        self.node_l = node_l
        self.founder_l = filter(lambda x: not x.parent_l, self.node_l)
        self.sequenced_n = len([node for node in node_l if node.sequenced_b == True])

    def gene_drop(self, pop_af):
        '''Set the founder genotypes.'''
        self.set_founder_genotype(pop_af)
        '''Get 1 founder per spousal pair.'''
        founder_no_spouse_l = []
        for founder in self.founder_l:
            if founder.spouse not in founder_no_spouse_l:
                founder_no_spouse_l.append(founder)
        '''Gene drop.'''
        for founder_node in founder_no_spouse_l:
            self.gene_drop_dfs(founder_node)

    def set_founder_genotype(self, pop_af):
        for founder in self.founder_l:
            founder.genotype = np.random.binomial(2,pop_af)

    def gene_drop_dfs(self, start_node):
        visited, stack = set(), [start_node]
        while stack:
            node = stack.pop(0)
            if node not in visited:
                #print "Visiting " + node.id
                if node.genotype == None:
                    node.genotype = self.get_genotype_from_parents(node)
                #print node.genotype
                visited.add(node)
                # new nodes are added to the start of stack
                #if node.spouse != None:
                    #stack = [node.spouse] + stack
                if len(node.children_l) > 0:
                    stack = node.children_l + stack
        #print ""
        return visited

    def get_genotype_from_parents(self, node):
        if node.parent_l[0].genotype == None or node.parent_l[1].genotype == None:
            return None
        else:
            return self.simulate_segregation(node.parent_l[0].genotype) + self.simulate_segregation(node.parent_l[1].genotype)

    def simulate_segregation(self, genotype):
        if genotype == 0:
            return 0
        elif genotype == 2:
            return 1
        elif genotype == 1:
            return np.random.randint(0,2)

    def print_all_genotypes(self, sequenced_only_b=False):
        if sequenced_only_b == False:
                print " ".join([node.id + ":" + str(node.genotype) for node in self.node_l])
        else:
                print " ".join([node.id + ":" + str(node.genotype) for node in self.node_l if node.sequenced_b == True])

    def print_all_info(self):
        for node in self.node_l:
            summary_str_l = ["ID: " + node.id,
                             "Founder: " + str(node in self.founder_l),
                             "Parents: " + ",".join([parent_node.id for parent_node in node.parent_l]),
                             "Spouse: " + node.spouse.id if node.spouse != None else "Spouse: ",
                             "Children: " + ",".join([child_node.id for child_node in node.children_l]),
                             "Genotype: " + str(node.genotype),
                             "Sequenced: " + str(node.sequenced_b)]
            print "; ".join(summary_str_l)
        print ""

    def get_carrier_allele_count(self):
        return sum((node.genotype for node in self.node_l if node.sequenced_b == True))

    def reset_genotypes(self):
        for node in self.node_l:
            node.genotype = None
            

class Cohort(object):
    
    def __init__(self, cohort_tsv):
        self.logger = Logger()
        self.node_generator = NodeGenerator()
        cohort_df = pd.read_csv(cohort_tsv, sep="\t", dtype=str)
        cohort_df["SEQUENCED"] = cohort_df["SEQUENCED"].apply(lambda x: eval(x))
        self.logger.log("Making family tree objects...")
        self.fam_tree_l = cohort_df.groupby("FAMILY").apply(self.make_fam_tree).tolist()
    
    def make_fam_tree(self, ped_df):
        '''
        
        Args:
            ped_df (pandas.core.frame.DataFrame): contains the pedigree information for the family.
        Returns:
            family_tree (FamilyTree obj): represents the family.
        '''
        node_l = self.node_generator.convert_ped_df_to_node_l(ped_df)
        return FamilyTree(node_l)
    
    def gene_drop(self, pop_af, cohort_af, gene_drop_n):
        '''
        Perform gene dropping across the cohort and return the proportion of iterations in which the simulated allele frequency is less than or equal to the cohort frequency.
        
        Args:
            pop_af (float): population allele frequency.
            cohort_af (float): cohort allele frequency.
            gene_drop_n (int): number of iterations to perform.
        
        Returns:
            cohort_enriched_p (float): proportion of iterations in which the simulated allele frequency is less than or equal to the cohort frequency.
        '''
        
        self.logger.log("Start gene drop with pop_af={0}, cohort_af={1} & gene_drop_n={2}.".format(pop_af,cohort_af,gene_drop_n))
        lt_thresh_n = 0
        for n in xrange(gene_drop_n):
            total_allele_count,carrier_allele_count = 0,0
            for fam_tree in self.fam_tree_l:
                fam_tree.gene_drop(pop_af)
                carrier_allele_count += fam_tree.get_carrier_allele_count()
                total_allele_count += fam_tree.sequenced_n*2
                fam_tree.reset_genotypes()
            gene_drop_af = carrier_allele_count/float(total_allele_count)
            if cohort_af <= gene_drop_af:
                lt_thresh_n += 1
        cohort_enriched_p = lt_thresh_n/float(gene_drop_n)
        self.logger.log("p-value={0}".format(cohort_enriched_p))
        
        return cohort_enriched_p


class Logger(object):
    
    def log(self, txt): 
        '''Prints a time-stamped text string.
        
        Args:
            txt (str): text string to print.
        '''
        
        timestamped_log_str = datetime.now().strftime('%Y-%m-%d %H:%M:%S') + ": " + txt
        print timestamped_log_str
