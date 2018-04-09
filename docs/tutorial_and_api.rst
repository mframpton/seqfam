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

Having cloned the repository, run the command python setup.py install.

Tutorial
########

gene_drop
=========

For a rare variant to be considered potentially causal of a particular trait/disease based on in silico analysis, it must satisfy various criteria, such as being biologically plausible and predicted to be pathogenic.
The user can run analyses with the gene_drop.py and pof.py modules to acquire additional evidence.
Given the structure of the families, Monte Carlo gene dropping can assess whether a variant is enriched in the cohort relative to the general population, and assuming the trait/disease is more prevalent in the cohort, such enrichment supports causality.
The user can use the pof.py module to identify variants which are carried by most or all affected members of a family, or even which segregate between affected and unaffected members.
The authors are unaware of existing software packages for performing these analyses in familial cohorts, so gene_drop.py and pof.py fulfil this need.
The gene_drop.py module can be considered complementary to the RVsharing R package (Bureau et al., 2014), which calculates the probability of multiple affected relatives sharing a rare variant under the assumption of no disease association or linkage.

By default, for each variant of interest, the gene_drop.py module performs 10,000 iterations of gene dropping in the familial cohort.
In each iteration it gene drops in each family once, seeding the founder genotypes based on the population allele frequency.
It then calculates the resulting simulated cohort allele frequency from samples specified by the user i.e. those which were sequenced.
After completing all iterations of gene dropping, gene_drop.py outputs the proportion of iterations in which the true cohort allele frequency is less than or equal to the simulated cohort allele frequency: a low proportion, e.g. < 5%, is evidence of enrichment.

The module gene drops in a family in the following way.
First, it assigns a genotype (number of copies of the mutant allele) to each founder using a Binomial distribution where the number of trials is 2 and the probability of success in each trial is the population allele frequency.
Hence the founders are assumed to be unrelated.
It then performs a depth-first traversal starting from each founder (1 per spousal pair), and for heterozygotes, uses a random number generator to determine which parental allele to pass onto the child.
Thus, every individual in the family is assigned a genotype.

gene_burden
===========

This module implements the Combined Multivariate and Collapsing (CMC) burden test (Li & Leal, 2008) for detecting rare causal variants, where the multivariate test is a log-likelihood ratio test.
The user can supply covariates to control for potential confounders such as divergent ancestry.
This burden test should be applied to unrelated samples, and hence is of no use for cohorts containing few families.
However, for cohorts containing a relatively large number of families, a sufficient number of unrelated cases can be extracted and combined with a separate set of unrelated controls.
Burden tests aggregate rare variants in a gene or functional unit into a single score (Li & Leal, 2008; Madsen & Browning, 2009; Morris & Zeggini, 2010; Price et al., 2010), and are one broad class of statistical methods which combine the effects of rare variants in order to increase power over single marker approaches.
Sequence kernel association testing (SKAT) (Wu et al., 2011) is another widely-used sub-category of such methods.
In general, burden testing is more powerful than SKAT when a large proportion of variants are causal and are all deleterious/protective.

To use the gene_burden.py module, the user must first read the various required data into Pandas data frames.
These data include variant annotations by which to group (e.g. gene/functional unit) and aggregate (e.g. population allele frequency), and the genotypes, affection status and covariates for the unrelated samples.
The user can specify multiple categories in which to aggregate variants (e.g. into population allele frequency ranges of 0–1% and 1–5%), and variants outside these categories (e.g. more common variants) remain unaggregated.
An aggregated variant category takes the value 0 or 1.
For each variant group, having aggregated the variants, gene_burden.py will perform a multivariate test, which is a log-likelihood ratio test based on Wilk’s theorem:

                     X2=2(llh0 − llh1); df = dfh1 − dfh0

where ll is log-likelihood, h1 is the alternative hypothesis, h0 is the null hypothesis and df is degrees of freedom.
Specifically, it is a log-likelihood ratio test on null and alternative hypothesis logit models where the dependent variable is derived from affection status, the variant variables (aggregated and/or unaggregated) are independent variables in the alternative model and the covariates are independent variables in both.
The logit models are fitted using the Broyden–Fletcher–Goldfarb–Shanno (BFGS) algorithm.

Since aggregation by population allele frequency is common, the module provides a function to assign variants to allele frequency ranges based on multiple columns (derived from different population databases such as gnomAD (Lek et al., 2016)).
The user specifies a preference order in case the variant is absent from the most preferred database(s).

pof
===

For a rare variant to be considered potentially causal of a particular trait/disease based on in silico analysis, it must satisfy various criteria, such as being biologically plausible and predicted to be pathogenic.
The user can run analyses with the gene_drop.py and pof.py modules to acquire additional evidence.
Given the structure of the families, Monte Carlo gene dropping can assess whether a variant is enriched in the cohort relative to the general population, and assuming the trait/disease is more prevalent in the cohort, such enrichment supports causality.
The user can use the pof.py module to identify variants which are carried by most or all affected members of a family, or even which segregate between affected and unaffected members.
The authors are unaware of existing software packages for performing these analyses in familial cohorts, so gene_drop.py and pof.py fulfil this need.
The gene_drop.py module can be considered complementary to the RVsharing R package (Bureau et al., 2014), which calculates the probability of multiple affected relatives sharing a rare variant under the assumption of no disease association or linkage.

For each family, the user can use the pof.py module to define a variant pattern of occurrence rule and check whether any supplied variants pass.
The rule can specify a minimum value for the proportion of affected members (As) who are carriers (A.carrier.p), and/or a minimum difference between the proportion of As and unaffected members (Ns) who are carriers (AN.carrier.diff).
Constraints for the number of genotyped As and Ns can also be added.

As an illustrative example, consider a cohort in which families can be categorised as follows based on their number of As and Ns:

    1. “A4N1”: ≥ 4 As and ≤ 1 N

    2. “A3N2”: ≥ 3 As and ≥ 2 Ns

For the A4N1 families, the user may be interested in variants carried by all As and so require A.carrier.p = 1, while for the A3N2 families, they may be interested in variants which are more prevalent in As than Ns and so require AN.carrier.diff ≥ 0.5.

relatedness
===========

The potential for genetic discovery in DNA sequencing data is reduced when samples are mislabelled.
Hence, necessary quality control steps include identifying duplicates, and in the case of familial samples, verifying the ascertained familial relationships described in the pedigrees.
The relatedness.py module facilitates these quality control steps and is used in conjunction with KING software (Manichaikul et al., 2010).
Given genotypes for relatively common variants, KING can efficiently calculate a kinship coefficient for each sample pair.
The relatedness.py module can then map each kinship coefficient to a degree of relationship and check it corresponds with the pedigree. KING is often already part of NGS analysis pipelines, so incorporating relatedness.py is straightforward. Peddy (Pedersen & Quinlan, 2017) is an alternative which does not require KING.

As input, the relatedness.py module requires pedigree information and a file containing kinship coefficients outputted by KING.
For each sample pair, the module will map the pedigree information and kinship coefficient to an expected and observed degree of relationship respectively.
The mapping from kinship coefficient to relationship is as specified in KING documentation: > 0.354 for duplicate samples/monozygotic twins, 0.177–0.354 for 1st degree relatives, 0.0884–0.177 for 2nd degree relatives, 0.0442–0.0884 for 3rd degree relatives, and < 0.0442 for unrelated.
The user can change this mapping if they wish.

sge
===

The final module, sge.py, has general utility in running analyses of NGS data (and indeed any “big data”) on computer clusters.
Many NGS data analyses can be cast as “embarrassingly parallel problems” and hence executed more efficiently on a computer cluster using a “MapReduce pattern”: the overall task is decomposed into independent sub-tasks (“map” tasks), then the map tasks run in parallel and after their completion, a “reduce” action merges/filters/summarises the results.
For example, gene burden testing across the whole exome can be decomposed into independent sub-tasks by splitting the exome into sub-units e.g. chromosomes.
Sun Grid Engine (SGE) is a widely used batch-queueing system, and analyses can be performed in a MapReduce pattern on SGE via so-called array jobs.
The sge.py module can be used to automatically create the scripts required for submitting and running an array job.

To use the sge.py module, the user must first create lists of map tasks, map tasks requiring execution and reduce tasks.
The map tasks requiring execution are map tasks which have not previously completed successfully and hence need to run.
Given these lists, the sge.py module can create all necessary scripts/files for submitting and running an array job.
This includes scripts for all map tasks, a text file specifying that only the map tasks requiring execution should run, and a master executable submit script for submitting the array job to the job scheduler.

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
 