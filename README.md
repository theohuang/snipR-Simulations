# snipR-Simulations
 Simulations for the SNIP algorithm

This repository contains simulation results for the paper "SNIP: An Adaptation of Sorted Neighborhood Methods for Deduplicating Pedigree Data". There are 2 types of simulations: (1) simulations to motivate the importance of deduplication (Section 2 in the paper), and (2) simulations to assess the performance of the SNIP algorithm (Section 3 in the paper).

1. To replicate the simulation results for the first type, run "Sim_Effect_a01.R" (allele frequency of 0.1) or "Sim_Effect_a001.R" (allele frequency of 0.01). The results are in "Res_Sim_Effect_a01.RData" and "Res_Sim_Effect_a001.RData", respectively. In addition, to add the step of running SNIP on these pedigrees, run "Sim_Effect_af01_dedup.R" or "Sim_Effect_af001_dedup.R", with results in "Res_Sim_Effect_af01_dedup.RData" and "Res_Sim_Effect_af001_dedup.RData", respectively.
2. To replicate the simulation results for the second type, run "Sim_Dup.R" (for simulated data with duplicates) or "Sim_Dup_nodup.R" (for simulated data without duplicates). These files were created to run on the Harvard FAS Research Computing cluster, and may have to be altered to replicate the results on your local machine. The results are in the "Sim_Param_Results" subfolder. "Sim_Param_Analysis.R" obtains performance metrics based on the results. There are also simulations on smaller pedigrees ("Sim_Dup_small.R"), with results in "res_sim_small.RData".

All code is written by either Theodore Huang or Matthew Ploenzke.
