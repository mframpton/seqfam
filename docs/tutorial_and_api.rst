******
seqfam
******

Introduction
############

The *seqfam* package is primarily designed for analysing next generation sequencing (NGS) DNA data from families with known pedigree information in order to identify rare variants that are potentially causal of a disease/trait of interest.
It uses the popular and versatile Pandas library, and can be straightforwardly integrated into existing analysis code/pipelines.
*Seqfam* can be used to verify pedigree information, to perform Monte Carlo gene dropping, to undertake regression-based gene burden testing, and to identify variants which segregate by affection status in families via user-defined pattern of occurrence rules.
Additionally, it can generate scripts for running analyses in a *MapReduce pattern* on a computer cluster, something which is usually desirable in NGS data analysis and indeed *big data* analysis in general.

Requirements and installation
#############################

*Seqfam* is compatible with Windows, Mac OX X and Linux operating systems. It is coded using Python 3.6 but can also be run by Python 2.7. It requires the following packages:

* pandas==0.20.3
* scipy==0.19.1
* natsort==5.1.1
* numpy==1.13.3
* setuptools==38.4.0
* statsmodels==0.8.0

Having cloned the repository, run the command :code:`python setup.py install`.

Tutorial
########

This section describes the functionality and methods employed by seqfam’s 5 modules, which are:

    1. *gene_drop.py*: Monte Carlo gene dropping;

    2. *pof.py*: variant pattern of occurrence in families;

    3. *gene_burden.py*: regression-based gene burden testing;

    4. *relatedness.py*: identification of duplicates and verification of ascertained pedigree information via kinship coefficients;

    5. *sge.py*: Sun Grid Engine (SGE) array job creation.

Figure 1 provides a visual representation of modules 1–4.

.. figure:: module_flowchart.png
    :align: center
    :alt: alternate text
    :figclass: align-center

    Panel A represents the Cohort.gene_drop method in the *gene_drop.py* module which performs Monte Carlo gene dropping.
    On a single iteration, for each family the algorithm seeds founder genotypes based on the variant population allele frequency and then gene drops via depth-first traversals.
    Having done this for all families, a simulated cohort allele frequency (AF) is calculated and following many iterations (e.g. 10,000), a p-value, the proportion of iterations where cohort AF < simulated cohort AF, is outputted.
    Panel B represents the Pof.get_family_pass_name_l method in the *pof.py* module.
    Prior to calling the method, each family is assigned a variant pattern of occurrence in family (POF) rule.
    The method then takes a variant’s genotypes and returns families whose POF rule is passed.
    Panel C represents the CMC.do_multivariate_tests method in the *gene_burden.py* module.
    This method takes sample affection status and variant genotypes across multiple genes, plus optionally covariates such as ancestry PCA coordinates.
    For each gene, the method aggregates the variants by allele frequency, constructs null and alternative hypothesis logit models which may include the covariates, and then performs a log-likelihood ratio test.
    Panel D represents the Relatedness.get_exp_obs_df method in the *relatedness.py* module.
    For input, this takes pedigree information and kinship coefficients from KING for each within-family sample pair.
    It maps these data to expected and observed degrees of relationship respectively, returning a data frame.

The repository contains additional scripts in src/examples which demonstrate the functionality of the modules on example data, including files in the data directory.
The scripts are *1_example_gene_drop.py*, *2_example_pof.py*, *3_example_gene_burden.py*, *4_example_relatedness.py*, and *5_example_sge.py*.
The reader can also refer to Table 1 for a summary of the main user functions of the 5 seqfam modules, which includes their input/output.
Data in the example data files are derived from the whole exome sequencing of a large cohort of over 200 families with inflammatory bowel disease (unpublished study1).

.. table:: Table 1. Summary of main user functions in seqfam modules

   +-------------+-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | Module      | Method (Class)                | Description                                                                         | Input                                                                                                        | Output                                                           |
   +=============+===============================+=====================================================================================+==============================================================================================================+==================================================================+
   | gene_drop   | gene_drop (Cohort)            | Monte Carlo gene dropping                                                           | Cohort tsv file (pedigree info), variant population AF, cohort AF, list of samples genotyped, # interactions | p-value                                                          |
   +-------------+-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | pof         | get_family_pass_name_l        | Variant POF with respect to affected (& unaffected) members                         | Variant POF rule & genotypes                                                                                 | List of families whose POF rules is passed by variant            |
   +-------------+-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | gene_burden | do_multivariate_tests (CMC)   | Regression-based gene burden testing                                                | Files containing samples, genotypes & covariates files; output path                                          | Data frame and csv file containing burden test results           |
   +-------------+-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   |             | find_duplicates (Relatedness) | Identify duplicates from kinship coefficient                                        | King file                                                                                                    | List of duplicates                                               |
   | relatedness +-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   |             | get_exp_obs_df (Relatedness)  | Map pedigrees & kinship coefficients to expected & observed degrees of relationship | Cohort tsv, KING within-family sample pair kinship coefficient file                                          | Data frame of expected & observed degrees of relationship        |
   +-------------+-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | sge         | make_map_reduce_jobs (SGE)    | Make computer cluster array job scripts.                                            | Filename prefix, lists of map tasks, map tasks to execute and reduce tasks.                                  | Scripts required to run array job including master submit script |
   +-------------+-------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+

gene_drop
=========

For a rare variant to be considered potentially causal of a particular trait/disease based on in silico analysis, it must satisfy various criteria, such as being biologically plausible and predicted to be pathogenic.
The user can run analyses with the *gene_drop.py* and *pof.py* modules to acquire additional evidence.
Given the structure of the families, Monte Carlo gene dropping can assess whether a variant is enriched in the cohort relative to the general population, and assuming the trait/disease is more prevalent in the cohort, such enrichment supports causality.
The user can use the *pof.py* module to identify variants which are carried by most or all affected members of a family, or even which segregate between affected and unaffected members.
The authors are unaware of existing software packages for performing these analyses in familial cohorts, so *gene_drop.py* and *pof.py* fulfil this need.
The *gene_drop.py* module can be considered complementary to the *RVsharing* *R* package :cite:`Bureau2014` which calculates the probability of multiple affected relatives sharing a rare variant under the assumption of no disease association or linkage.

By default, for each variant of interest, the *gene_drop.py* module performs 10,000 iterations of gene dropping in the familial cohort.
In each iteration it gene drops in each family once, seeding the founder genotypes based on the population allele frequency.
It then calculates the resulting simulated cohort allele frequency from samples specified by the user i.e. those which were sequenced.
After completing all iterations of gene dropping, *gene_drop.py* outputs the proportion of iterations in which the true cohort allele frequency is less than or equal to the simulated cohort allele frequency: a low proportion, e.g. < 5%, is evidence of enrichment.

The module gene drops in a family in the following way.
First, it assigns a genotype (number of copies of the mutant allele) to each founder using a Binomial distribution where the number of trials is 2 and the probability of success in each trial is the population allele frequency.
Hence the founders are assumed to be unrelated.
It then performs a depth-first traversal starting from each founder (1 per spousal pair), and for heterozygotes, uses a random number generator to determine which parental allele to pass onto the child.
Thus, every individual in the family is assigned a genotype.

The only input file for *1_example_gene_drop.py* is cohort.tsv, which contains the pedigree information for an example familial cohort.
This cohort has 3,608 samples from 251 families, and the complexity of the families, calculated as *2n-f* where n and f are the number of non-founders and founders respectively (Abecasis et al., 2002), has median 9 and range 0–103.
The cohort.csv file is in fam file format (Purcell et al., 2007), meaning it has 1 row per individual and 6 columns for family ID, person ID, father, mother, sex and affection.

The *1_example_gene_drop.py* script first creates a Cohort object from *cohort.tsv*, which stores each family tree, then calls the *gene_drop* method with the following arguments: *pop_af* and *cohort_af* are the allele frequency of a particular variant in the general population and cohort respectively, *sample_genotyped_l* is the list of cohort samples with a genotype, and *gene_drop_n* is the number of iterations of gene dropping to perform.
Hence the samples in *sample_genotyped_l* are used by the user to calculate *cohort_af*, and by the method to calculate the simulated cohort allele frequencies.
The method returns a p-value.
The script calls the gene_drop method with ascending values for *cohort_af*, and so descending p-values are returned.

.. code-block:: python
       
   import os
   from seqfam.gene_drop import Cohort

   #Create cohort object from cohort.tsv file.
   data_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data"))
   cohort_tsv = os.path.join(data_dir,"cohort.tsv")
   cohort = Cohort(cohort_tsv)

   #Get a list of genotyped individuals.'''
   cohort_sample_l = cohort.get_all_sample_l()
   sample_genotyped_l = list(filter(lambda sample: sample.find("p") == -1,cohort_sample_l))

   #Do gene dropping.'''
   for cohort_af in [0.025,0.03,0.035,0.04]:
       p = cohort.gene_drop(0.025, cohort_af, sample_genotyped_l, 1000)

pof
===

Like *gene_drop.py*, the *pof.py* module can provide additional evidence for whether particular rare variants are causal of a particular trait/disease. It is intended for identifying variants which are carried by most or all affected members of a family, or even which segregate between affected and unaffected members.

For each family, the user uses the *pof.py* module to define a variant pattern of occurrence rule and check whether any supplied variants pass. The rule can specify a minimum value for the proportion of affected members (*As*) who are carriers (*A_carrier_p*), and/or a minimum difference between the proportion of *As* and unaffected members (*Ns*) who are carriers (*AN_carrier_diff*). Constraints for the number of genotyped *As* and *Ns* can also be added, (*A_n_min* and *N_n_min* respectively).

The example script *2_example_pof.py* provides an illustrative example.
It first creates a couple of *Family* objects to represent 2 families and their variant pattern of occurrence rule.

.. code-block:: python
   
   import pandas as pd
   from seqfam.pof import Family,Pof

   ...

   family_1 = Family("1","A3N2",["1_1","1_2","1_3"],["1_4","1_5"],A_n_min=3,N_n_min=2,AN_carrier_diff=0.5)
   family_2 = Family("2","A4N1",["2_10","2_11","2_12","2_13","2_14"],["2_15"],A_n_min=4,N_n_min=1,A_carrier_p=1.0)

Family 1 is specified as having 3 affected and 2 unaffected members, and its pattern of occurrence rule requires *AN_carrier_diff* to be 0.5.
Family 2 has 4 affecteds and 1 unaffected, and a rule requiring all the affecteds to be carriers.
The rule in both families requires all members to be genotyped (see the *A_n_min* and *N_n_min* parameters).

The script next makes the genotypes for a hypothetical variant in a Pandas Series called *variant_genotypes_s*.
Finally, it creates a *Pof* object from the 2 *Family* objects, and then calls the *get_family_pass_name_l* method to obtain a list of the families whose pattern of occurrence rule is passed by this variant.

.. code-block:: python

   family_l = [family_1,family_2]
   pof = Pof(family_l)
   family_pass_l = pof.get_family_pass_name_l(variant_genotypes_s)
   print(family_pass_l)

gene_burden
===========

This module implements the Combined Multivariate and Collapsing (CMC) burden test :cite:`Li2008` for detecting rare causal variants, where the multivariate test is a log-likelihood ratio test.
The user can supply covariates to control for potential confounders such as divergent ancestry.
This burden test should be applied to unrelated samples, and hence is of no use for cohorts containing few families.
However, for cohorts containing a relatively large number of families, a sufficient number of unrelated cases can be extracted and combined with a separate set of unrelated controls.
Burden tests aggregate rare variants in a gene or functional unit into a single score (:cite:`Li2008`; :cite:`Madsen2009`; :cite:`Morris2010`; :cite:`Price2010`, and are one broad class of statistical methods which combine the effects of rare variants in order to increase power over single marker approaches.
Sequence kernel association testing (SKAT) :cite:`Wu2011` is another widely-used sub-category of such methods.
In general, burden testing is more powerful than SKAT when a large proportion of variants are causal and are all deleterious/protective.

To use the *gene_burden.py* module, the user must first read the various required data into Pandas data frames.
These data include variant annotations by which to group (e.g. gene/functional unit) and aggregate (e.g. population allele frequency), and the genotypes, affection status and covariates for the unrelated samples.
The user can specify multiple categories in which to aggregate variants (e.g. into population allele frequency ranges of 0–1% and 1–5%), and variants outside these categories (e.g. more common variants) remain unaggregated.
An aggregated variant category takes the value 0 or 1.
For each variant group, having aggregated the variants, *gene_burden.py* will perform a multivariate test, which is a log-likelihood ratio test based on Wilk’s theorem:

.. math::

   \chi^2 = 2(ll_{h0} - ll_{h1}); df = df_{h1} - df_{h0}

where *ll* is log-likelihood, *h1* is the alternative hypothesis, *h0* is the null hypothesis and *df* is degrees of freedom.
Specifically, it is a log-likelihood ratio test on null and alternative hypothesis logit models where the dependent variable is derived from affection status, the variant variables (aggregated and/or unaggregated) are independent variables in the alternative model and the covariates are independent variables in both.
The logit models are fitted using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm.

Since aggregation by population allele frequency is common, the module provides a function to assign variants to allele frequency ranges based on multiple columns (derived from different population databases such as gnomAD :cite:`Lek2016`.
The user specifies a preference order in case the variant is absent from the most preferred database(s).

The *3_example_gene_burden.py* script uses the *gene_burden.py* module to perform CMC tests on example data: it performs 1 CMC test per gene where variants in the population allele frequency ranges of 0-1% and 1-5% are aggregated, and any variants ≥ 5% remain unaggregated.
The input files are in the data/gene_burden directory: samples.csv, genotypes.csv and optionally covariates.csv.
The samples.csv file contains the samples’ ID and affection status where 2 indicates a case and 1 a control, and genotypes.csv contains 1 row per variant with columns for the sample genotypes, and for variant grouping and aggregation e.g. gene and population allele frequency.
A sample’s genotype is the number of alternate alleles which it carries (0-2).
The covariates.csv file contains the covariates to control for, which in this case are ancestry PCA coordinates.

The script first reads samples.csv into a Pandas Series, and genotypes.csv and covariates.csv into Pandas DataFrames.
These data frames are indexed by variant ID and covariate name respectively.
Having created a CMC object, the script calls its assign_vars_to_pop_frq_cats method in order to map the variants to the desired population allele frequency ranges.
Multiple population allele frequency columns (databases) are used here, ordered by descending preference.
The mapping is stored in a new column in the genotypes data frame.
Finally, the script calls the do_multivariate_tests method to perform the CMC tests, specifying the gene column for grouping the variants, and the new allele frequency range column for aggregation.
The results are written to a CSV (comma-separated values) file and returned in a data frame.
They include the number of variants in each aggregation category, the number of unaggregated variants (“unagg” column), the log-likelihood ratio test p-value with/without covariates (“llr_p” and “llr_cov_p”), and the coefficient/p-value for each aggregated variant variable (“_c” and “_p”).

.. code-block:: python

   import os
   import pandas as pd
   from seqfam.gene_burden import CMC
   import sys


   #Set paths.
   data_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data","gene_burden"))
   genotypes_path = os.path.join(data_dir,".".join(["variant_genotypes","csv"]))
   samples_path = os.path.join(data_dir,".".join(["samples","csv"]))
   covariates_path = os.path.join(data_dir,".".join(["covariates","csv"]))
   results_path = os.path.join(data_dir,".".join(["cmc","results","csv"]))

   #Read the samples into a Series.
   sample_s = pd.read_csv(samples_path, dtype=str, index_col="Sample ID")
   sample_s["Affection"] = sample_s["Affection"].astype(int)
   sample_s = sample_s[sample_s != 0]

   #Read the variant annotations + genotypes into a DataFrame '''
   variant_col,gene_col = "VARIANT_ID","Gene"
   pop_frq_col_l = ["database1_AF","database2_AF","database3_AF"]
   geno_df = pd.read_csv(genotypes_path, dtype=str, usecols=[variant_col,gene_col] + pop_frq_col_l + sample_s.index.tolist(), index_col=variant_col)
   geno_df.loc[:,pop_frq_col_l] = geno_df.loc[:,pop_frq_col_l].apply(pd.to_numeric, axis=1)
   geno_df.loc[:,sample_s.index] = geno_df.loc[:,sample_s.index].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
   geno_df.loc[:,sample_s.index] = geno_df.loc[:,sample_s.index].fillna(0)

   #Read the covariates into a DataFrame.
   covar_df = None if covariates_path == None else pd.read_csv(covariates_path, index_col=0)

   #Do gene burden (CMC) tests.
   cmc = CMC()
   geno_df = cmc.assign_variants_to_pop_frq_cats(geno_df, pop_frq_col_l, {"rare":0.01,"mod_rare":0.05})
   cmc_result_df = cmc.do_multivariate_tests(sample_s, geno_df, group_col=gene_col, agg_col="pop_frq_cat", agg_val_l=["rare","mod_rare"], covar_df=covar_df, results_path=results_path)
   print(cmc_result_df)

relatedness
===========

The potential for genetic discovery in DNA sequencing data is reduced when samples are mislabelled.
Hence, necessary quality control steps include identifying duplicates, and in the case of familial samples, verifying the ascertained familial relationships described in the pedigrees.
The *relatedness.py* module facilitates these quality control steps and is used in conjunction with KING software :cite:`Manichaikul2010`.
Given genotypes for relatively common variants, KING can efficiently calculate a kinship coefficient for each sample pair.
The *relatedness.py* module can then map each kinship coefficient to a degree of relationship and check it corresponds with the pedigree.
KING is often already part of NGS analysis pipelines, so incorporating *relatedness.py* is straightforward.
Peddy :cite:`Pedersen2017` is an alternative which does not require KING.

As input, the *relatedness.py* module requires a pedigree information file and a kinship coefficients file from KING, either containing within or between-family sample pairs.
The example script *4_example_relatedness.py* uses data/cohort.tsv for the former and data/relatedness/king.kinship.ibs (within-family sample pairs) for the latter.
It creates a Relatedness object with these paths, and then calls the object's find_duplicates and get_exp_obs_df methods.

.. code-block:: python
   
   from seqfam.relatedness import Relatedness
   import pandas as pd
   ...
   relatedness = Relatedness(wf_file=wf_file,cohort_tsv=cohort_tsv,bf_file=None)

   #Within-family duplicates.
   wf_duplicate_l = relatedness.find_duplicates(bf_b=False)
   print(wf_duplicate_l)
   #Between-family duplicates... (Uncomment if you have a bf_file).
   #bf_duplicate_l = relatedness.find_duplicates(bf_b=True)
   #print(bf_duplicate_l)

   #Expected versus observed within-family relationships.
   exp_obs_df = relatedness.get_exp_obs_df()

The find_duplicates method returns a list of any duplicate samples, and the get_exp_obs_df method returns a Pandas DataFrame containing the expected and observed degree of relationship for each within-family sample pair.

>>> print(wf_duplicate_l)
['171b_1448,171b_1449']
>>> print(exp_obs_df)
                 EXP_REL  Kinship OBS_REL
FAMILY ID1  ID2                          
1      44   47         1   0.2584       1
2      6    20         4   0.0390       4
            21         4   0.0688       3
       20   21         4   0.0051       4
3      501  838        1   0.2052       1
            844        2   0.1081       2
       838  844        1   0.2572       1
...

The user can modify the mapping from kinship coefficient to relationship degree, but by default it is as specified in KING documentation: > 0.354 for duplicate samples/monozygotic twins, 0.177–0.354 for 1st degree relatives, 0.0884–0.177 for 2nd degree relatives, 0.0442–0.0884 for 3rd degree relatives, and < 0.0442 for unrelated. The final line of the script prints the sample pairs which have a different expected and observed degree of relationship.

.. code-block:: python
   
   print(exp_obs_df.loc[(exp_obs_df["EXP_REL"]!=exp_obs_df["OBS_REL"]) & (pd.notnull(exp_obs_df["Kinship"])),:])

sge
===

The sge.py module has general utility in analysing NGS data, and indeed any big data on computer clusters. Many NGS data analyses can be cast as "embarassingly parallel problems" and hence executed more efficiently on a computer cluster via a "MapReduce pattern": the overall task is decomposed into independent sub-tasks ("map" tasks) which run in parallel, and on their completion, a "reduce" action merges/filters/summarises the results. For example, gene burden testing across the whole exome can be decomposed into independent sub-tasks by splitting the exome into sub-units e.g. chromosomes. Sun Grid Engine (SGE) is a widely used batch-queueing system, and analyses can be performed in a MapReduce pattern on SGE via array jobs. Given a list of map tasks and the reduce task(s), the sge.py module can create the scripts for submitting and running an array job.

The script 5_example_sge.py provides an example. 
It first makes lists of map tasks (map_task_l) and map tasks to execute (map_task_exec_l) via the custom get_map_task_l method (see the script), and then a reduce task string (reduce_tasks).
While map_task_l contains all map tasks, map_task_exec_l contains the subset which have not yet completed successfully and hence need to run.

.. code-block:: python

   from seqfam.sge import SGE   
   ...
   print("Making map and reduce tasks...")
   chr_l = [str(chrom) for chrom in range(1,23)] + ["X","Y"]
   [map_task_l, map_task_exec_l] = get_map_task_l(chr_l)
   reduce_tasks = "\n".join(["python 2_merge_results.py","python 3_summarise_results.py"])

Next, the script creates an SGE object which stores the directory where job scripts will be written (the variable script_dir which here has the value data/sge). Finally it calls the object's make_map_reduce_jobs method with the following arguments: a job script name prefix (here "test"), map_task_l, map_task_exec_l and reduce_tasks.

.. code-block:: python
   
   sge = SGE(script_dir)
   sge.make_map_reduce_jobs("test", map_task_l, reduce_tasks, map_task_exec_l)

This writes the job scripts, and were they for a real array job (they are not), the user could submit it to the job scheduler by running the master executable submit script data/sge/submit_map_reduce.sh. The generated file test.map_task_exec.txt specifies which map tasks to run (map_tasks_exec_l).

References
==========

.. bibliography:: references.bib

API reference
#############

gene_drop
=========

.. automodule:: gene_drop
   :members:

gene_burden
===========

.. automodule:: gene_burden
   :members:
   
pof
===

.. automodule:: pof
   :members:

relatedness
===========

.. automodule:: relatedness
   :members:

sge
===

.. automodule:: sge
   :members:

misc
====

.. automodule:: misc
   :members:
 