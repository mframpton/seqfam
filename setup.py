from setuptools import setup

setup(
    name='Seqfam',
    version='0.1',
    packages=['seqfam',],
    license='Creative Commons Attribution-Noncommercial-Share Alike license',
    long_description=open('README.md').read(),
    author="Matthew Frampton",
    author_email="mjeframpton@gmail.com",
    description="The seqfam package is primarily designed for analysing next generation sequencing DNA data from families with pedigree information in order to identify rare variants that are potentially causal of a disease/trait.",
    url='https://github.com/mframpton/seqfam',
    install_requires=['pandas==0.20.3','scipy==0.19.1','natsort==5.1.1','numpy==1.13.3','setuptools==38.4.0','statsmodels==0.8.0'],
)
