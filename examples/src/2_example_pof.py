import pandas as pd
from seqfam.pof import Family,Pof


#Create some example families.    
family_1 = Family("1","A3N2",["1_1","1_2","1_3"],["1_4","1_5"],A_n_min=3,N_n_min=2,AN_p_diff_min=0.5)
family_1.log_info()
family_2 = Family("2","A4N1",["2_10","2_11","2_12","2_13","2_14"],["2_15"],A_n_min=4,N_n_min=1,A_p_min=1.0)
family_2.log_info()

#Create some genotypes.
variant_genotypes_s = pd.Series(data=["1","1","1","1","0","1","1","1","1","1","1"],
                                index=["1_1","1_2","1_3","1_4","1_5","2_10","2_11","2_12","2_13","2_14","2_15"])

#Create a Pof object and check whether this variant passes in the example families.
family_l = [family_1,family_2]
pof = Pof(family_l)
family_pass_l = pof.get_family_pass_name_l(variant_genotypes_s)
print(family_pass_l)
