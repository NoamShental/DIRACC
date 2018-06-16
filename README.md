**DIRACC feasibility study results**
====================================

This package contains the basic Matlab code used to test the DIRACC approach for efficient detection of carriers of rare copy number variations (CNV). The approach combines the mathematical field of compressed sensing and digital PCR (dPCR) for detecting carriers of germline deletions (heterozygous or homozygous) or of germline duplications of rare CNVs. 

DIRACC stands for DIgital PCR for rare CNV carrier detection via Compressed Sensing. 

The code allows to design DIRACC measuring matrices and simulate their performance in different settings of either deletions or duplications. Simulations model the dPCR noise in different settings and thus allow to check DIRACC in specific set points. 


Detecting carriers of rare germline deletions
------------------------------------------------ 

The basic function in the case of deletions is DIRACC_deletions.m whose input parameters are described below, followed by an example function call. 
Basically, for given a specific set point, random input vectors containing a certain number of homozygous or heterozygous carriers are selected. dPCR droplets are then simulated according to the error model (see manuscript) and DIRACC is applied to detect the carriers. The function returns the success rate of each set point.


‫**‬Input parameters‫**‬

*M*: the Bernoulli sensing matrix, whose columns correspond to the tested individuals, and the rows are the pools. An entry in M is one when a certain individual takes part of a certain pool, and zero otherwise 
                                  
*NumHeterozygousCarriers*: the number of carriers of a heterozygous deletion out of all samples

*NumHomozygousCarriers*: the number of carriers of a homozygous deletion out of all samples

*D*: the total number of dPCR droplets used

*NumOfTests*: the number of simulations in each setup point 

*LimitNumOfPools*: A limit on the number of pools used (i.e. only a certain subset of pools is applied) 

*AddDNAPreperationNoise*: The standard deviation of DNA concentration, referred to as DNA preparation error in the manuscript. For a noiseless case, set this parameter to zero

*AddSampelingNoise*: true/false, whether to consider sampling noise 

*MinFr*: minimal value of f to consider (f is the average fraction of occupied droplets, as calibrated by diluting the DNA according to the available number of droplets D)

*StepFr*: the delta between consecutive values of f

*MaxFr*: the maximum value of f to consider

*MinSuccess*: the minimal required success rate. Testing a certain set point would terminate if the rate is lower than this value (e.g. 97%)

*Verbose*: enable progress prints (for debug)

*‫*‬Output parameters*‫*‬
*P* - the vector containing the % success of exact detections for each tested fraction. For example if MinFr = 0.3 , MaxFr = 0.51 , StepFr = 0.1 and the results were
      98% success for f=0.3 , 100% success for f=0.4 and 97% success for f=0.5 the output vector would be [ 0.97 , 1.00 , 0.98 ]  

*Pc* - Percentage of success in the same format as P. However, in this case successful detection is declared if the correct carrier is found regardless of its deletion state (homozygous or heterozygous)

**Example function call**: 

[ P , Pc ] = DIRACC_deletions( B1024_36 , 2 , 0 , 10^7 , 500 , 36 , 0.1 , true , 0.3 , 0.1 , 0.71 , 0.97 , true)

B1024_36: Is a Bernoulli sensing matrix of 1024 individuals with 36 pools, available by loading Bernulli1024_36_0125.mat (a Bernoulli matrix with probability p=0.125)

NumHeterozygousCarriers = 2, i.e, two heterozygous carriers

NumHomozygousCarriers = 0, i.e., zero homozygous carriers

D = 10000000, i.e., 10,000,000 droplets per pool.

NumOfTests = 500, i.e, run 500 simulations for each value of f

LimitNumOfPools = 36, all provided 36 pools are applied

AddDNAPreperationNoise = 0.1, i.e., DNA preparation error ~N(0,0.1^2)

AddSampelingNoise = true, sampling noise is taken into consideration

MinFr = 0.3, start running the simulations from fraction of 30%

StepFr = 0.1, continue to the next fraction in steps of 10%

MaxFr = 0.71, stop the simulation when the fraction exceeds 71% 

MinSuccess = 0.97, Simulations of a specific set point are terminated if the number of errors exceeds 3% of the total number of simulations fo that 
set point 

Verbose = true, i.e., enables debug prints while running the simulation

Detecting carriers of rare germline duplications
--------------------------------------------------- 

The basic function in the case of deletions is DIRACC_duplications.m, whose input parameters are described below, followed by an example function call. The function is analogous to the case of deletions, while applying a duplication-relevant error model. Basically, for given a specific set point, random input vectors containing a certain number of carriers of specific duplications are selected. dPCR droplets are then simulated according to the error model (see manuscript) and DIRACC is applied to detect the carriers. The function returns the success rate of each set point.

*‫*‬Input parameters‫*‬*

*M*: the Bernoulli sensing matrix, whose columns correspond to the tested individuals, and the rows are the pools. An entry in M is one when a certain individual takes part of a certain pool, and zero otherwise 

*CarrierVec*: A vector containing ordered pairs of the number of carriers and their corresponding copy number. For example, [1 3 10 2], correspond to a single carrier with three additional copies (i.e., five copies in total), and ten carriers each having two additional gene copies (i.e, each of these carriers holds five gene copies)

*D*: the total number of dPCR droplets used

*NumOfTests*: the number of simulations in each setup point 

*LimitNumOfPools*: A limit on the number of pools used (i.e. only a certain subset of pools is applied) 

*AddDNAPreperationNoise*: The standard deviation of DNA concentration, referred to as DNA preparation error in the manuscript. For a noiseless case, set this parameter to zero

*AddSampelingNoise*: true/false, whether to consider sampling noise 

*AllowedDelta*: Allowed delta from the correct copy number value. For example in case AllowedDelta = -1 then finding the right carrier but assigning a copy number smaller by one is considered correct

*MinFr*: minimal value of f to consider (f is the average fraction of occupied droplets, as calibrated by diluting the DNA according to the available number of droplets D)

*StepFr*: the delta between consecutive values of f

*MaxFr*: the maximum value of f to consider

*MinSuccess*: the minimal required success rate. Testing a certain set point would terminate if the rate is lower than this value (e.g. 97%)

*Verbose*:enable progress prints (for debug)

**Output parameters**
*P* - the vector containing the % success of CNV carrier detection for different f values. Detection is considered correct if the correct sample is detected while the copy number may differ according to "AllowedDelta". 

*Pc* - Percentage of success in the same format as P. However, in this case successful detection is declared if the correct carrier is found while ignoring the exact copy number (any value Y>0 is considered correct).

**Example function call**: 

CarrierVec=[1,4];
[ P , Pc ] = DIRACC_duplications( B1024_36 , CarrierVec , 10000000 , 500 , 36 , 0.1 , true , -2 , 0.3 , 0.1 , 0.91 , 0.97 , true  )

B1024_36: Is a Bernoulli sensing matrix of 1024 individuals with 36 pools, available by loading Bernulli1024_36_0125.mat

CarrierVec = [1,4] i.e., one carrier having six gene copies (Y=4) 

D = 10000000, i.e., 10,000,000 droplets per pool.

NumOfTests = 500, i.e, run 500 simulations for each value of f

LimitNumOfPools = 36, all provided 36 pools are applied

AddDNAPreperationNoise = 0.1, i.e., DNA preparation error ~N(0,0.1^2)

AddSampelingNoise = true, sampling noise is taken into consideration

AllowedDelta = -2, finding the right carrier but assigning a copy number smaller by at most two is considered correct 

MinFr=0.3 , start running the simulations from fraction of 30%


StepFr=0.1 , continue to the next fraction in steps of 10%

MaxFr=0.71, stop the simulation when we exceed the fraction of 71% 

MinSuccess = 0.97, Simulations of a specific set point are terminated if the number of errors exceeds 3% of the total number of simulations fo that set point 

Verbose = true, i.e., enables debug prints while running the simulation


References
------------
Please send any comments/bug reports to: Noam Shental, shental@openu.ac.il


