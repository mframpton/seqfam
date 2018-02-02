# CohPy

Requirements
============

*CohPy* is compatible with Windows, Mac OS X and Linux operating systems. It is coded using Python 3.6, but can also be run by Python 2.7. It requires various data analysis libraries/packages such as *Pandas*. The user can acquire almost all of these by downloading and installing the *Anaconda* python distribution. The *StatsModels* module, which is not in the *Anaconda* distribution, is also required. Having cloned the repository, the user should add the repository *src* directory to their python path.

Modules
=======

*CohPy* has 5 modules:
1.	*gene_drop.py*: Monte Carlo gene dropping;
2.	*pof.py*: variant pattern of occurrence in families;
3.	*gene_burden.py*: regression-based gene burden testing;
4.	*relatedness.py*: identification of duplicates and verification of ascertained pedigree information via kinship coefficients;
5.	*sge.py*: Sun Grid Engine (SGE) array job creation.

Tests
=====

The repository contains additional scripts in *src/examples* which demonstrate the functionality of the modules on example data, including files in the *data* directory. They are *1_test_gene_drop.py*, *2_test_pof.py*, *3_test_gene_burden.py*, *4_test_relatedness.py*, and *5_test_sge.py*. The user should first change directory into *src/examples* before running the scripts because some use relative paths to files in the *data* directory.

If the package is used to generate data for a publication, you must cite the following article:

Matthew Frampton, Elena Schiff, Nikolas Pontikos, Adam P. Levine, Anthony Segal (2018) CohPy: a python package for analysis of next generation sequencing DNA data in families. *(Under review)*
