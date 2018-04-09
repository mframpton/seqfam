from natsort import natsorted
from seqfam.misc import Logger


class Family(object):
	'''Represents a Family with various attributes: name, category (with respect to number of affected and unaffected members), IDs of affected and unaffected
	members, plus attributes relating to conditions which must be satisfied for a variant to "pass" in the family. These are the minimum number of affected members
	who are carriers, the minimum number of unaffecteds who are carriers, the minimum proportion of affecteds who are carriers, and the minimum difference in the
	proportion of affecteds and unaffecteds who are carriers. 
	'''

	def __init__(self, name, category, A_l, N_l, A_n_min=0, N_n_min=0, A_p_min=None, AN_p_diff_min=None):
		self.logger = Logger()
		self.name = name
		self.category = category
		self.A_l = natsorted(A_l)
		self.N_l = natsorted(N_l)
		self.A_n_min = A_n_min
		self.N_n_min = N_n_min
		self.A_p_min = A_p_min
		self.AN_p_diff_min = AN_p_diff_min


	def log_info(self):
		'''Log the object attributes.'''
		
		info_str_l = ["\nFamily name:\t{0}".format(self.name), "Category:\t{0}".format(self.category), 
					"Affecteds:\t{0}".format(",".join(self.A_l)),"Unaffecteds:\t{0}".format(",".join(self.N_l)),
					"A_p_min:\t{0}".format(self.A_p_min), "AN_p_diff_min:\t{0}\n".format(self.AN_p_diff_min)]
		info_str = "\n".join(info_str_l)
		self.logger.log(info_str)
		

	def pass_po(self, variant_genotypes_s, no_call="NA", carrier_call=["1","2"]):
		'''Check whether a variant passes in the family.
		
		Args:
			| variant_genotypes_s (Series): the genotypes of family members for the variant of interest.
			| no_call (str): how a no-call is represented in the genotype data.
			| carrier_call (list of strs): genotypes which correspond to carrying the variant.
		
		Returns:
			boolean: whether the variant passes.
		'''
		
		#Check the A_n_min and N_n_min checks are satisfied.
		if (variant_genotypes_s[self.A_l] != no_call).sum() < self.A_n_min or (variant_genotypes_s[self.N_l] != no_call).sum() < self.N_n_min:
			return False
		A_p = None
		if self.A_p_min != None or self.AN_p_diff_min != None:
			A_geno_count_s = variant_genotypes_s.ix[self.A_l].value_counts(normalize=True, dropna=False)
			A_p = A_geno_count_s.ix[carrier_call].sum()
		if self.A_p_min != None: 
			if A_p < self.A_p_min:
				return False
		if self.AN_p_diff_min != None:
			N_geno_count_s = variant_genotypes_s.ix[self.N_l].value_counts(normalize=True, dropna=False)
			N_p = N_geno_count_s.ix[carrier_call].sum()
			if A_p - N_p < self.AN_p_diff_min:
				return False
			
		return True


class Pof(object):
	'''Pattern of (variant) occurrence: stores a list of family objects, and has a function to check if the genotypes for a variant of interest satisfy the
	sepcified pattern of occurence criteria in any of these families.'''

	def __init__(self, family_l):
		self.family_l = family_l

	def get_family_pass_name_l(self, variant_genotypes_s):
		'''Checks whether the genotypes for a variant of interest satisfy the specified pattern of occurrence criteria (pass) in any of the supplied families.
		
		Args:
			variant_genotypes_s (Series): contains the genotypes for the variant of interest for all individuals in the families contained in the family_l attribute.
		Returns:
			pass_l (list of Family objects): the list of families in which the variant of interest passes.
		'''
		
		pass_l = [family.name for family in self.family_l if family.pass_po(variant_genotypes_s)]
		return pass_l
