#!/bin/bash

R CMD BATCH --vanilla '--args sparsity.level=15 sample.size=200' Simulation_main.R R1.out&
R CMD BATCH --vanilla '--args sparsity.level=25 sample.size=200' Simulation_main.R R2.out&
R CMD BATCH --vanilla '--args sparsity.level=50 sample.size=200' Simulation_main.R R3.out&
R CMD BATCH --vanilla '--args sparsity.level=100 sample.size=200' Simulation_main.R R4.out&
R CMD BATCH --vanilla '--args sparsity.level=15 sample.size=500' Simulation_main.R R5.out&
R CMD BATCH --vanilla '--args sparsity.level=25 sample.size=500' Simulation_main.R R6.out&
R CMD BATCH --vanilla '--args sparsity.level=50 sample.size=500' Simulation_main.R R7.out&
R CMD BATCH --vanilla '--args sparsity.level=100 sample.size=500' Simulation_main.R R8.out&
R CMD BATCH --vanilla '--args sparsity.level=15 sample.size=1000' Simulation_main.R R9.out&
R CMD BATCH --vanilla '--args sparsity.level=25 sample.size=1000' Simulation_main.R R10.out&
R CMD BATCH --vanilla '--args sparsity.level=50 sample.size=1000' Simulation_main.R R11.out&
R CMD BATCH --vanilla '--args sparsity.level=100 sample.size=1000' Simulation_main.R R12.out&
