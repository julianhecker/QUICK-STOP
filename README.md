sequential.cpp:
Implementation of a nearly optimal Stopping rule for randomized tests


Usage of the QUICK-STOP tool:

 ./quick_stop p p_1 p_2 alpha seed number_of_replications
 
for example
 
 ./quick_stop 0.1 0.01 0.04 0.0000000001 1 10000
 
should give

'QUICK-STOP:
p: 0.1
p_1: 0.01
p_2: 0.04
alpha: 1e-10
average number of permutations/simulations: 181.796
lower bound: 159.371
ratio QUICK-STOP/lower bound: 1.1407'
