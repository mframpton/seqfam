from natsort import natsorted


class Family(object):

	def __init__(self, name, category, A_l, N_l, A_n_min=0, N_n_min=0, A_p_min=None, AN_p_diff_min=None):
		self.name = name
		self.category = category
		self.A_l = natsorted(A_l)
		self.N_l = natsorted(N_l)
		self.A_n_min = A_n_min
		self.N_n_min = N_n_min
		self.A_p_min = A_p_min
		self.AN_p_diff_min = AN_p_diff_min


	def print_info(self):
		print "Family name:\t{0}".format(self.name)
		print "Category:\t{0}".format(self.category)
		print "Affecteds:\t{0}".format(",".join(self.A_l))
		print "Unaffecteds:\t{0}".format(",".join(self.N_l))
		print "A_p_min:\t{0}".format(self.A_p_min)
		print "AN_p_diff_min:\t{0}\n".format(self.AN_p_diff_min)


	def pass_po(self, geno_s, no_call="NA", carrier_call=["1","2"]):
		
		#Check the A_n_min and N_n_min checks are satisfied.
		if (geno_s[self.A_l] != no_call).sum() < self.A_n_min or (geno_s[self.N_l] != no_call).sum() < self.N_n_min:
			return False
		A_p = None
		if self.A_p_min != None or self.AN_p_diff_min != None:
			A_geno_count_s = geno_s.ix[self.A_l].value_counts(normalize=True, dropna=False)
			A_p = A_geno_count_s.ix[carrier_call].sum()
		if self.A_p_min != None: 
			if A_p < self.A_p_min:
				return False
		if self.AN_p_diff_min != None:
			N_geno_count_s = geno_s.ix[self.N_l].value_counts(normalize=True, dropna=False)
			N_p = N_geno_count_s.ix[carrier_call].sum()
			if A_p - N_p < self.AN_p_diff_min:
				return False
			
		return True


class Pof(object):

	def __init__(self, family_l):
		self.family_l = family_l

	def get_family_pass_name_l(self, geno_s):
		pass_l = [family.name for family in self.family_l if family.pass_po(geno_s)]
		return pass_l
