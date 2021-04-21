# Deng_Benchmarking_Validation

Code used to verify results in Deng et al. https://doi.org/10.1007/s10008-019-04402-6 

Code imports data taken from the Deng paper for lithium concentration and stress in a cylindrical anode

We compute the analytical solutions to these results as given in Deng and independently verified by myself

We also compute the result for the stress using the finite element solver COMSOL

COMSOL and Deng data are imported from CSV files

Analytic solutions to concentration and stress are computed using functions due to the combersome nature and the occurence of Bessel functions in the solution

Main file is DengAnalytic3.m
