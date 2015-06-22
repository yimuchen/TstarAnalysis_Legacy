# TstarAnalysis

This code includes the code for the tstar mass reconsturcting tstarmass. 

## Prerequisites 

This code must be ran on a ntugrid machine with a root environment setup:


## Standard procedure:
   1. From this directory, run the command 
   ```
   ./prep.csh 
   ```
   This generates the file required for parallel computation.
   2. Run the command 
   ```
   ./run.csh
   ```
   log files are stored in `*/log*` files. Processes are all placed in the background.
   When it finishes running, the final lines should read how much time the process took.
   3. Run the command
   ```
   root -l -b -q merge.C++
   ```
   to merge the results from the various files into single `SumNTuple.root` file.
   4. From the directory `./TgammaAnalysis_TemplateFit/` run the command
   ```
   root -l -q -b TemplateFit_v21.C++
   ```
   for the fitting
