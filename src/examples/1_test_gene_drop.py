import os
import sys
from gene_drop.gene_drop import Cohort


data_dir = os.path.join("..","..","data","gene_drop")
cohort_tsv = os.path.join(data_dir,"cohort.tsv")
cohort = Cohort(cohort_tsv)
#print help(cohort)
#print help(cohort.gene_drop)

p = cohort.gene_drop(0.025, 0.025, 1000)
p = cohort.gene_drop(0.025, 0.03, 1000)
p = cohort.gene_drop(0.025, 0.035, 1000)
p = cohort.gene_drop(0.025, 0.04, 1000)
