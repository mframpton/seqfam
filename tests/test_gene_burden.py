import unittest
from seqfam.gene_burden import CMC
import pandas as pd
import numpy as np

class Test(unittest.TestCase):


    def setUp(self):
        self.cmc = CMC()
        self.geno_df = pd.DataFrame({"Variant":["var{0}".format(i) for i in range(1,7)],
                                     "Gene":["gene1","gene1","gene1","gene2","gene2","gene2"],
                                     "DB1_AF":[0.01,0.004,0.05,np.NaN,np.NaN,np.NaN],
                                     "DB2_AF":[0.04,0.11,0.005,0.001,np.NaN,0.046],
                                     "DB3_AF":[0.005,0.11,0.009,0.001,np.NaN,np.NaN],
                                     "sample1":[1,0,1,0,1,0],
                                     "sample2":[0,1,0,1,0,1]})
        self.geno_df.set_index("Variant", inplace=True)
        self.cmc.group_col = "Gene"
        self.cmc.agg_col = "pop_frq_cat"
        self.cmc.sample_s = pd.Series(data=[1,2],index=["sample1","sample2"])

    def tearDown(self):
        del self.geno_df

    def test_assign_variants_to_pop_frq_cats(self):
        expected_s = pd.Series(data=["mod_rare","rare","var3","rare","rare","mod_rare"],
                               index=["var{0}".format(i) for i in range(1,7)])
        result_s = self.cmc.assign_variants_to_pop_frq_cats(self.geno_df,["DB1_AF","DB2_AF","DB3_AF"],{"rare":0.01,"mod_rare":0.05})["pop_frq_cat"]
        self.assertEqual(expected_s.tolist(),result_s.tolist())
        
    def test_aggregate_by_agg_col(self):
        self.geno_df = self.cmc.assign_variants_to_pop_frq_cats(self.geno_df,["DB1_AF","DB2_AF","DB3_AF"],{"rare":0.01,"mod_rare":0.05})
        self.geno_df = self.cmc.aggregate_by_agg_col(self.geno_df)
        self.assertEqual([1,0,1,0,1],self.geno_df["sample1"].tolist())
        self.assertEqual([0,1,0,1,1],self.geno_df["sample2"].tolist())
        self.assertEqual([1,1,1,1,2],self.geno_df["n"].tolist())

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()