% 
% The basic function in the case of deletions is DIRACC_deletions.m whose input parameters are described below, followed by an example function call. 
% Basically, for given a specific set point, random input vectors containing a certain number of homozygous or heterozygous carriers are selected. dPCR droplets are then simulated according to the error model (see manuscript) and DIRACC is applied to detect the carriers. The function returns the success rate of each set point.
% 
% 
% ** input parameters **
% 
% M: the Bernoulli sensing matrix, whose columns correspond to the tested individuals, and the rows are the pools. An entry in M is one when a certain individual takes part of a certain pool, and zero otherwise 
%                                   
% NumHeterozygousCarriers*: the number of carriers of a heterozygous deletion out of all samples
% 
% NumHomozygousCarriers*: the number of carriers of a homozygous deletion out of all samples
% 
% D: the total number of dPCR droplets used
% 
% NumOfTests: the number of simulations in each setup point 
% 
% LimitNumOfPools: A limit on the number of pools used (i.e. only a certain subset of pools is applied) 
% 
% AddDNAPreperationNoise: The standard deviation of DNA concentration, referred to as DNA preparation error in the manuscript. For a noiseless case, set this parameter to zero
% 
% AddSampelingNoise: true/false, whether to consider sampling noise 
% 
% MinFr: minimal value of f to consider (f is the average fraction of occupied droplets, as calibrated by diluting the DNA according to the available number of droplets D)
% 
% StepFr: the delta between consecutive values of f
% 
% MaxFr: the maximum value of f to consider
% 
% MinSuccess: the minimal required success rate. Testing a certain set point would terminate if the rate is lower than this value (e.g. 97%)
% 
% Verbose: enable progress prints (for debug)
% 
% ** Output parameters **
% P - the vector containing the % success of exact detections for each tested fraction. For example if MinFr = 0.3 , MaxFr = 0.51 , StepFr = 0.1 and the results were
%       98% success for f=0.3 , 100% success for f=0.4 and 97% success for f=0.5 the output vector would be [ 0.97 , 1.00 , 0.98 ]  
% 
% Pc - Percentage of success in the same format as P. However, in this case successful detection is declared if the correct carrier is found regardless of its deletion state (homozygous or heterozygous)

function [ P , Pc ] = DIRACC_deletions( M , NumHeterozygousCarriers , NumHomozygousCarriers , D , NumOfTests , LimitNumOfPools , AddDNAPreparationNoise , AddSamplingNoise , MinFr , StepFr , MaxFr , MinSuccess , Verbouse)

NumOfFractions = ceil((MaxFr - MinFr)/StepFr);
P = zeros (NumOfFractions,1);
Pc = zeros(NumOfFractions,1);
index = 1;
[~,NumOfPersons] = size(M);

% adding the cvx path
addpath('./cvx');
cvx_setup;
    
%the loop for running the code test and keeping the resuilts for each
%fraction
for Fr = MinFr:StepFr:MaxFr

    Count = 0;
    CarrierCount = 0;
    Disp = zeros(1,4);
    for i=1:NumOfTests
        
        % generate an the input vector 
        I = GenerateInputVector( NumOfPersons , NumHeterozygousCarriers , NumHomozygousCarriers );
        
        %set the sorted carrier vector 
        [Isort , Iindex] = sort(I);
        Iindex = Iindex(Isort>0);
        Iindex = sort(Iindex);
        [Hr, Ht] = SimPoolsHits( M , 0 , I , D , Fr , LimitNumOfPools , AddDNAPreparationNoise , AddSamplingNoise );
        [CNV] = ReconstructCnvVec( M , 0 , Hr , Ht , 0.2 );
        
        %if the result vector is equal to the
        %input vector then increment the test success counter        
        if CNV == I
            Count = Count + 1;                        
            CarrierCount = CarrierCount + 1;
            
        %if the result is not equal check if the carrier is detected     
        else            
            [SortCNV,CarrierIndex]=sort(CNV);
            CarrierIndex = CarrierIndex(SortCNV > 0);
            CarrierIndex = sort(CarrierIndex);
            if size(CarrierIndex) == size(Iindex) 
                if CarrierIndex == Iindex
                    CarrierCount = CarrierCount + 1;
                end
            end
        end
        %display current progress and check if we succeeded limitation 
        if  mod(i,5) == 0 
            Disp(1) = i;
            Disp(2) = Count/i;
            Disp(3) = CarrierCount/i;
            
            %calculating max posible success assuming all next test will succeed
            Disp(4) = (NumOfTests - (i - CarrierCount))/NumOfTests;
            
            %if we reached the minimal success ratio terminate
            if Disp(4) < MinSuccess
                if Verbouse == true
                    disp(' reached minimal success ratio : test terminated')
                end
                break;
            end
            
            if Verbouse == true 
                disp(Disp);
            end
        end
        
    end % for i=1:NumOfTests
    
    % set the success value only if all tests were executed , i.e if the 
    % ratio passes the minimum  
    if i == NumOfTests
        P(index) = Count/NumOfTests;
        Pc(index) = CarrierCount/NumOfTests;
    end
    index = index + 1;
    
    if Verbouse == true
        disp( Fr );
    end
end
end

