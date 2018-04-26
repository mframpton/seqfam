******
seqfam
******

Introduction
############

The :py:mod:`seqfam` package is primarily designed for analysing next generation sequencing (NGS) DNA data from families with known pedigree information in order to identify rare variants that are potentially causal of a disease/trait of interest.
It uses the popular and versatile `Pandas <https://pandas.pydata.org/>`_ library, and can be straightforwardly integrated into existing analysis code/pipelines.
:py:mod:`Seqfam` can be used to verify pedigree information, to perform Monte Carlo gene dropping, to undertake regression-based gene burden testing, and to identify variants which segregate by affection status in families via user-defined pattern of occurrence rules.
Additionally, it can generate scripts for running analyses in a *MapReduce pattern* on a computer cluster, something which is usually desirable in NGS data analysis and indeed *big data* analysis in general.

Requirements and installation
#############################

:py:mod:`Seqfam` is compatible with Windows, Mac OX X and Linux operating systems.
It is coded using Python 3.6 but can also be run by Python 2.7.
It requires the following packages (listed in *requirements.txt*):

* pandas==0.20.3
* scipy==0.19.1
* natsort==5.1.1
* numpy==1.13.3
* setuptools==38.4.0
* statsmodels==0.8.0

Run the following commands to clone and install from GitHub.

.. code-block:: bash
   
   $ git clone https://github.com/mframpton/seqfam
   $ cd seqfam
   $ pip install -r requirements.txt
   $ python setup.py install

Tutorial
########

This section describes the functionality and methods employed by :py:mod:`seqfam’s` 5 modules, which are:

    1. :py:mod:`gene_drop`: Monte Carlo gene dropping;

    2. :py:mod:`pof`: variant pattern of occurrence in families;

    3. :py:mod:`gene_burden`: regression-based gene burden testing;

    4. :py:mod:`relatedness`: identification of duplicates and verification of ascertained pedigree information via kinship coefficients;

    5. :py:mod:`sge`: Sun Grid Engine (SGE) array job creation.

Figure 1 provides a visual representation of modules 1–4.

.. figure:: module_flowchart.png
    :align: center
    :alt: alternate text
    :figclass: align-center

    Panel A represents the :py:meth:`Cohort.gene_drop` method in the :py:mod:`gene_drop` module which performs Monte Carlo gene dropping.
    On a single iteration, for each family the algorithm seeds founder genotypes based on the variant population allele frequency (AF) and then gene drops via depth-first traversals.
    Having done this for all families, a simulated cohort AF is calculated and following many iterations (e.g. 10,000), a p-value, the proportion of iterations where cohort AF < simulated cohort AF, is outputted.
    Panel B represents the :py:meth:`Pof.get_family_pass_name_l` method in the :py:mod:`pof` module.
    Prior to calling the method, each family is assigned a variant pattern of occurrence in family (POF) rule.
    The method then takes a variant’s genotypes and returns families whose POF rule is passed.
    Panel C represents the :py:meth:`CMC.do_multivariate_tests` method in the :py:mod:`gene_burden` module.
    This method takes sample affection status and variant genotypes across multiple genes, plus optionally covariates such as ancestry PCA coordinates.
    For each gene, the method aggregates the variants by AF, constructs *h0* and *h1* logit models which may include the covariates, and then performs a log-likelihood ratio test.
    Panel D represents the :py:meth:`Relatedness.get_exp_obs_df` method in the :py:mod:`relatedness` module.
    For input, this takes pedigree information and kinship coefficients from `KING <http://people.virginia.edu/~wc9c/KING/>`_ for each within-family sample pair.
    It maps these data to expected and observed degrees of relationship respectively, returning a :py:mod:`DataFrame`.

The repository contains additional scripts in :file:`src/examples` which demonstrate the functionality of the modules on example data, including files in the :file:`data` directory.
The scripts are :file:`1_example_gene_drop.py`, :file:`2_example_pof.py`, :file:`3_example_gene_burden.py`, :file:`4_example_relatedness.py`, and :file:`5_example_sge.py`.
The reader can also refer to Table 1 for a summary of the main user functions of the 5 :py:mod:`seqfam` modules.
Data in the example data files are derived from the whole exome sequencing of a large cohort of over 200 families with inflammatory bowel disease.

.. table:: Table 1. Summary of main user functions in :py:mod:`seqfam` modules

   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | Method                                             | Description                                                                         | Input                                                                                                        | Output                                                           |
   +====================================================+=====================================================================================+==============================================================================================================+==================================================================+
   | :py:meth:`gene_drop.Cohort.gene_drop`              | Monte Carlo gene dropping                                                           | Cohort fam file (pedigree info), variant population AF, cohort AF, list of samples genotyped, # interactions | p-value                                                          |
   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | :py:meth:`pof.Pof.get_family_pass_name_l`          | Variant POF with respect to affected (& unaffected) members                         | Variant POF rule & genotypes                                                                                 | List of families whose POF rules is passed by variant            |
   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | :py:meth:`gene_burden.CMC.do_multivariate_tests`   | Regression-based gene burden testing                                                | Files containing samples, genotypes & covariates files; output path                                          | Data frame and csv file containing burden test results           |
   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | :py:meth:`relatedness.Relatedness.find_duplicates` | Identify duplicates from kinship coefficient                                        | KING sample pair kinship coefficient file                                                                    | List of duplicates                                               |
   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | :py:meth:`relatedness.Relatedness.get_exp_obs_df`  | Map pedigrees & kinship coefficients to expected & observed degrees of relationship | Cohort fam, KING within-family sample pair kinship coefficient file                                          | Data frame of expected & observed degrees of relationship        |
   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+
   | :py:meth:`sge.SGE.make_map_reduce_jobs`            | Make computer cluster array job scripts.                                            | Filename prefix, lists of map tasks, map tasks to execute and reduce tasks.                                  | Scripts required to run array job including master submit script |
   +----------------------------------------------------+-------------------------------------------------------------------------------------+--------------------------------------------------------------------------------------------------------------+------------------------------------------------------------------+

gene_drop
=========

For a rare variant to be considered potentially causal of a particular trait/disease based on in silico analysis, it should satisfy various criteria, such as being biologically plausible and predicted to be pathogenic.
The :py:mod:`gene_drop` module can be used to further assess candidate variants via Monte Carlo gene dropping.

Given the structure of the families, Monte Carlo gene dropping can indicate whether a variant is enriched in the cohort relative to the general population, and assuming the trait/disease is more prevalent in the cohort, such enrichment supports causality.
The :py:mod:`gene_drop` module can be considered complementary to the `RVsharing <https://cran.r-project.org/web/packages/RVsharing/index.html>`_ *R* package :cite:`Bureau2014` which calculates the probability of multiple affected relatives sharing a rare variant under the assumption of no disease association or linkage.

The module requires a pedigree file in `fam format <https://www.cog-genomics.org/plink2/formats#fam>`_ as input.
The example script :file:`1_example_gene_drop.py` shows how to use :py:mod:`gene_drop` with the pedigrees in :file:`data/cohort.fam`.
It first creates a :py:obj:`gene_drop.Cohort` object from :file:`cohort.fam` which stores all of the pedigrees as trees.

.. code-block:: python

   from seqfam.gene_drop import Cohort
   ...
   cohort_fam = os.path.join(data_dir,"cohort.fam")
   cohort = Cohort(cohort_fam)

For a hypothetical variant of interest, the script then specifies:

1. allele frequency in the general population (*pop_af*) is 0.025;
2. the subset of samples which have genotypes (*sample_genotyped_l*).

Now the gene dropping can be performed via the :py:meth:`gene_drop.Cohort.gene_drop` method.
The script uses the method to assess whether increasing cohort allele frequencies (*cohort_af*) indicate enrichment relative to the general population.
For each *cohort_af*, the method returns an enrichment p-value (*p*), and so as *cohort_af* increases, *p* decreases.

.. code-block:: python
   
   pop_af = 0.025
   for cohort_af in [0.025,0.03,0.035,0.04]:
          p = cohort.gene_drop(pop_af, cohort_af, sample_genotyped_l, 1000)

The method gene drops in a family in the following way.
First, it assigns a genotype (number of copies of the mutant allele) to each founder using a Binomial distribution where the number of trials is 2 and the probability of success in each trial is *pop_af*.
Hence the founders are assumed to be unrelated.
It then performs a depth-first traversal starting from each founder (1 per spousal pair), and for heterozygotes, uses a random number generator to determine which parental allele to pass onto the child.
Thus, every individual in the family is assigned a genotype.

By default, for each variant of interest, the method performs 10,000 iterations of gene dropping in the familial cohort.
In each iteration it gene drops in each family once and then calculates the resulting simulated cohort allele frequency from the genotyped samples (*sample_genotyped_l*).
After completing all iterations of gene dropping, the method returns the p-value (*p*), which is the proportion of iterations where *cohort_af* is less than or equal to the simulated cohort allele frequency.
A low proportion, e.g. < 5%, can be taken as evidence of enrichment.

pof
===

Like :py:mod:`gene_drop`, the :py:mod:`pof` module can provide additional evidence for whether particular rare variants are causal of a particular trait/disease.
It is intended for identifying variants which are carried by most or all affected members of a family (*As*), or even which segregate between *As* and unaffected members (*Ns*).

For each family, the user uses the :py:mod:`pof` module to define a variant pattern of occurrence (POF) rule and check whether any supplied variants pass.
The rule can specify a minimum value for the proportion *As* who are carriers (*A_carrier_p*), and/or a minimum difference between the proportion of *As* and *Ns* who are carriers (*AN_carrier_diff*).
Constraints for the number of genotyped *As* and *Ns* can also be added, (*A_n_min* and *N_n_min* respectively).

The example script :file:`2_example_pof.py` provides an illustrative example.
It first creates a couple of :py:obj:`pof.Family` objects to represent 2 families and their POF rule.

.. code-block:: python
   
   import pandas as pd
   from seqfam.pof import Family,Pof

   ...

   family_1 = Family("1","A3N2",["1_1","1_2","1_3"],["1_4","1_5"],A_n_min=3,N_n_min=2,AN_carrier_diff=0.5)
   family_2 = Family("2","A4N1",["2_10","2_11","2_12","2_13","2_14"],["2_15"],A_n_min=4,N_n_min=1,A_carrier_p=1.0)

Family 1 is specified as having 3 *As* and 2 *Ns*, and its pattern of occurrence rule requires *AN_carrier_diff* to be 0.5.
Family 2 has 4 *As* and 1 *N*, and a rule requiring all the *As* to be carriers.
The rule in both families requires all members to be genotyped (see the *A_n_min* and *N_n_min* parameters).

The script next makes the genotypes for a hypothetical variant in a Pandas :py:mod:`Series` called *variant_genotypes_s*.
Finally, it creates a :py:obj:`pof.Pof` object from the 2 :py:obj:`pof.Family` objects, and then calls the :py:meth:`pof.Pof.get_family_pass_name_l` method to obtain a list of the families whose POF rule is passed by this variant.

.. code-block:: python

   family_l = [family_1,family_2]
   pof = Pof(family_l)
   family_pass_l = pof.get_family_pass_name_l(variant_genotypes_s)
   print(family_pass_l)

gene_burden
===========

The :py:mod:`gene_burden.py` module implements the Combined Multivariate and Collapsing (CMC) burden test :cite:`Li2008` for detecting rare causal variants, where the multivariate test is a log-likelihood ratio test.
The user can supply covariates to control for potential confounders such as divergent ancestry.
This burden test should be applied to unrelated samples, and hence is of no use for cohorts containing few families.
However, for cohorts containing a relatively large number of families, a sufficient number of unrelated cases can be extracted and pooled with a separate set of unrelated controls.
Burden tests aggregate rare variants in a gene or functional unit into a single score (:cite:`Li2008`; :cite:`Madsen2009`; :cite:`Morris2010`; :cite:`Price2010`, and are one broad class of statistical methods which combine the effects of rare variants in order to increase power over single marker approaches.
Sequence kernel association testing (SKAT) :cite:`Wu2011` is another widely-used sub-category of such methods.
In general, burden testing is more powerful than SKAT when a large proportion of variants are causal and are all deleterious/protective.

The :file:`3_example_gene_burden.py` script shows how to use the :py:mod:`gene_burden` module to perform CMC tests which control for covariates.
Here, we say that variants are *grouped* by the tested units (gene / other functional unit), and within the groups, they are *aggregated*, usually within population allele frequency (PAF) ranges.
Aggregation means that within each aggregation category (e.g. PAF < 1%), an individual sample is given the value 1 if it carries any variants, otherwise 0.
The example script performs 1 CMC test per gene (i.e. it groups variants by gene), where variants are aggregated within 2 PAF ranges: PAF < 1% and 1% <= PAF < 5% (any variants with PAF >= 5% remain unaggregated). 

The input files are in the :file:`data/gene_burden` directory: :file:`samples.csv`, :file:`genotypes.csv` and :file:`covariates.csv`.
The :file:`samples.csv` file contains the samples’ ID and affection status where 2 indicates a case and 1 a control.
The :file:`genotypes.csv` file can be created by combining genotypes from a `VCF file <http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/>`_ with variant annotations (e.g. from the `Variant Effect Predictor <https://www.ensembl.org/info/docs/tools/vep/index.html>`_).
It contains 1 row per variant with columns for the sample genotypes (the number of alternate alleles), plus columns for variant grouping and aggregation e.g. gene and PAF.
The :file:`covariates.csv` file contains the covariates to control for, which in this instance are ancestry Principal Components Analysis (PCA) coordinates.

The script first reads :file:`samples.csv` into a Pandas :py:mod:`Series`, and :file:`genotypes.csv` and :file:`covariates.csv` into Pandas :py:mod:`DataFrames`.
These :py:mod:`DataFrames` are indexed by variant ID and covariate name respectively.

.. code-block:: python

   import pandas as pd
   from seqfam.gene_burden import CMC
   ...
   #Read the samples into a Series.
   sample_s = pd.read_csv(samples_path, dtype=str, index_col="Sample ID")
   sample_s["Affection"] = sample_s["Affection"].astype(int)
   sample_s = sample_s[sample_s != 0]

   #Read the variant annotations + genotypes into a DataFrame.
   variant_col,gene_col = "VARIANT_ID","Gene"
   pop_frq_col_l = ["database1_AF","database2_AF","database3_AF"]
   geno_df = pd.read_csv(genotypes_path, dtype=str, usecols=[variant_col,gene_col] + pop_frq_col_l + sample_s.index.tolist(), index_col=variant_col)
   geno_df.loc[:,pop_frq_col_l] = geno_df.loc[:,pop_frq_col_l].apply(pd.to_numeric, axis=1)
   geno_df.loc[:,sample_s.index] = geno_df.loc[:,sample_s.index].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
   geno_df.loc[:,sample_s.index] = geno_df.loc[:,sample_s.index].fillna(0)

   #Read the covariates into a DataFrame.
   covar_df = None if covariates_path == None else pd.read_csv(covariates_path, index_col=0)

Having created a :py:obj:`gene_burden.CMC` object, the script calls its :py:func:`gene_burden.CMC.assign_variants_to_pop_frq_cats` method in order to map the variants to the desired PAF range categories.
Multiple PAF columns (databases) are used here, ordered by descending preference.
The mapping is stored in a new *pop_frq_cat* column in the genotypes :py:mod:`DataFrame`.

.. code-block:: python

   cmc = CMC()
   geno_df = cmc.assign_variants_to_pop_frq_cats(geno_df, pop_frq_col_l, {"rare":0.01,"mod_rare":0.05})

Finally, the script calls the :py:meth:`gene_burden.CMC.do_multivariate_tests` method to perform the CMC tests, specifying the gene column for grouping the variants, and the new *pop_frq_cat* column for aggregation.

.. code-block:: python

   agg_col = "pop_frq_cat"
   cmc_result_df = cmc.do_multivariate_tests(sample_s, geno_df, group_col=gene_col, agg_col="pop_frq_cat", agg_val_l=["rare","mod_rare"], covar_df=covar_df, results_path=results_path)

For each gene, this method performs a multivariate test, which is a log-likelihood ratio test based on Wilk’s theorem:

.. math::

   \chi^2 = 2(ll_{h0} - ll_{h1}); df = df_{h1} - df_{h0}

where *ll* is log-likelihood, *h1* is the alternative hypothesis, *h0* is the null hypothesis and *df* is degrees of freedom.
Specifically, it is a log-likelihood ratio test on null and alternative hypothesis logit models where the dependent variable is derived from affection status, the variant variables (aggregated and/or unaggregated) are independent variables in the alternative model and the covariates are independent variables in both.
The logit models are fitted using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm.

The results are written to a CSV (comma-separated values) file (*results_path*) and also returned in a :py:mod:`DataFrame` (*cmc_result_df*).
They include the number of variants in each aggregation category (*rare*, *mod_rare*), the number of unaggregated variants (*unagg*), the proportion of affecteds and unaffecteds which have value 1 for each variant variable (*rare_aff_p* ... *unagg_unaff_p*), the log-likelihood ratio test p-value with/without covariates (*llr_p* and *llr_cov_p*), and the coefficient/p-value for each aggregated variant variable and covariate in the *h1* logit model (*rare_c*, *rare_p* ... *PC5_c*, *PC5_p*).

>>> print(cmc_result_df.head().to_string())
         rare  mod_rare  unagg  rare_aff_p  rare_unaff_p  mod_rare_aff_p  mod_rare_unaff_p  unagg_aff_p  unagg_unaff_p         llr_p  llr_cov_p     rare_c    rare_p  mod_rare_c  mod_rare_p     PC1_c         PC1_p     PC2_c         PC2_p     PC3_c         PC3_p     PC4_c         PC4_p     PC5_c     PC5_p
Gene                                                                                                                                                                                                                                                                                                            
GENE_81    17         4      0    0.083333      0.041667        0.340909          0.090909          NaN            NaN  8.776693e-11   0.000075   0.260649  0.607068    1.317017    0.000041  0.350426  2.719782e-10  0.742543  9.178604e-11 -0.286762  1.989026e-10 -0.933229  2.921982e-26 -0.196124  0.100419
GENE_10     1         0      0    0.000000      0.003788             NaN               NaN          NaN            NaN           NaN   0.009316 -11.325180  0.825768         NaN         NaN  0.347303  3.341802e-10  0.670772  1.993786e-09 -0.289419  9.064199e-11 -0.978531  2.720945e-28 -0.193831  0.096250
GENE_16     4         0      0    0.007576      0.034091             NaN               NaN          NaN            NaN           NaN   0.013386  -2.475089  0.014941         NaN         NaN  0.337994  7.855867e-10  0.698458  4.950072e-10 -0.282116  2.037413e-10 -0.972062  6.228652e-28 -0.154929  0.185172
GENE_59     3         0      0    0.000000      0.011364             NaN               NaN          NaN            NaN           NaN   0.014363 -10.955116  0.897996         NaN         NaN  0.333042  1.325734e-09  0.675311  1.539985e-09 -0.284952  1.417334e-10 -0.961612  3.308430e-28 -0.176737  0.128946
GENE_15     2         0      0    0.000000      0.007576             NaN               NaN          NaN            NaN           NaN   0.027556 -10.728126  0.894869         NaN         NaN  0.329435  1.331849e-09  0.697045  6.262195e-10 -0.281625  2.635041e-10 -0.956502  6.064632e-28 -0.185669  0.109201

If the user ran the tests without covariates, then the returned results :py:mod:`DataFrame` would not include the *llr_cov_p* and *PC* covariate columns, and the columns *rare_c* ... *mod_rare_p* would correspond to the coefficient/p-value in the *h0* logit model.

relatedness
===========

The potential for genetic discovery in DNA sequencing data is reduced when samples are mislabelled.
Hence, necessary quality control steps include identifying duplicates, and in the case of familial samples, verifying the ascertained familial relationships described in the pedigrees.
The :py:mod:`relatedness` module facilitates these quality control steps and is used in conjunction with `KING <http://people.virginia.edu/~wc9c/KING/>`_ software :cite:`Manichaikul2010`.
Given genotypes for relatively common variants, *KING* can efficiently calculate a kinship coefficient for each sample pair.
The :py:mod:`relatedness` module can then map each kinship coefficient to a degree of relationship and check it corresponds with the pedigree.
*KING* is often already part of NGS analysis pipelines, so incorporating :py:mod:`relatedness` is straightforward.
`Peddy <https://github.com/brentp/peddy>`_ :cite:`Pedersen2017` is an alternative which does not require *KING*.

As input, the :py:mod:`relatedness` module requires a pedigree information file in `fam format <https://www.cog-genomics.org/plink2/formats#fam>`_ and a kinship coefficients file from *KING*, either containing within or between-family sample pairs.
The example script :file:`4_example_relatedness.py` uses :file:`data/cohort.fam` for the former and :file:`data/relatedness/king.kinship.ibs` (within-family sample pairs) for the latter.
It creates a :py:obj:`relatedness.Relatedness` object with these paths, and then calls the object's :py:func:`relatedness.Relatedness.find_duplicates` and :py:func:`relatedness.Relatedness.get_exp_obs_df` methods.

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

The :py:meth:`relatedness.Relatedness.find_duplicates` method returns a list of any duplicate samples, and the :py:meth:`relatedness.Relatedness.get_exp_obs_df` method returns a Pandas :py:mod:`DataFrame` containing the expected and observed degree of relationship for each within-family sample pair.

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

The user can modify the mapping from kinship coefficient to relationship degree, but by default it is as specified in *KING* documentation: > 0.354 for duplicate samples/monozygotic twins, 0.177–0.354 for 1st degree relatives, 0.0884–0.177 for 2nd degree relatives, 0.0442–0.0884 for 3rd degree relatives, and < 0.0442 for unrelated.
The final line of the script prints the sample pairs which have a different expected and observed degree of relationship.

.. code-block:: python
   
   print(exp_obs_df.loc[(exp_obs_df["EXP_REL"]!=exp_obs_df["OBS_REL"]) & (pd.notnull(exp_obs_df["Kinship"])),:])

sge
===

The :py:mod:`sge` module has general utility in analysing NGS data, and indeed any big data on computer clusters.
Many NGS data analyses can be cast as "embarassingly parallel problems" and hence executed more efficiently on a computer cluster via a *MapReduce pattern*: the overall task is decomposed into independent sub-tasks (*map tasks*) which run in parallel, and on their completion, a *reduce action* merges/filters/summarises the results.
For example, gene burden testing across the whole exome can be decomposed into independent sub-tasks by splitting the exome into sub-units e.g. chromosomes.
Sun Grid Engine (SGE) is a widely used batch-queueing system, and analyses can be performed in a *MapReduce pattern* on SGE via *array jobs*.
Given a list of *map tasks* and the *reduce task(s)*, the :py:mod:`sge` module can create the scripts for submitting and running an *array job*.

The script :file:`5_example_sge.py` provides an example. 
It first makes lists of *map tasks* (*map_task_l*) and *map tasks* to execute (*map_task_exec_l*) via the custom *get_map_task_l* function (see the script), and then a *reduce task* string (*reduce_tasks*).
While *map_task_l* contains all *map tasks*, *map_task_exec_l* contains the subset which have not yet completed successfully and hence need to run.

.. code-block:: python

   from seqfam.sge import SGE   
   ...
   print("Making map and reduce tasks...")
   chr_l = [str(chrom) for chrom in range(1,23)] + ["X","Y"]
   [map_task_l, map_task_exec_l] = get_map_task_l(chr_l)
   reduce_tasks = "\n".join(["python 2_merge_results.py","python 3_summarise_results.py"])

Next, the script creates an :py:obj:`sge.SGE` object which stores the directory where job scripts will be written (the variable *script_dir* which here has the value :file:`data/sge`).
Finally it calls the object's :py:meth:`sge.SGE.make_map_reduce_jobs` method with the following arguments: a job script name prefix (here *test*), *map_task_l*, *map_task_exec_l* and *reduce_tasks*.

.. code-block:: python
   
   sge = SGE(script_dir)
   sge.make_map_reduce_jobs("test", map_task_l, reduce_tasks, map_task_exec_l)

This writes the job scripts, and were they for a real array job (they are not), the user could submit it to the job scheduler by running the master executable submit script :file:`submit_map_reduce.sh`.
The generated file :file:`test.map_task_exec.txt` specifies which map tasks to run (*map_tasks_exec_l*).

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
 