% The basic function in the case of deletions is DIRACC_duplications.m, whose input parameters are described below, followed by an example function call. The function is analogous to the case of deletions, while applying a duplication-relevant error model. Basically, for given a specific set point, random input vectors containing a certain number of carriers of specific duplications are selected. dPCR droplets are then simulated according to the error model (see manuscript) and DIRACC is applied to detect the carriers. The function returns the success rate of each set point.
% 
% ** input parameters **
% 
% M: the Bernoulli sensing matrix, whose columns correspond to the tested individuals, and the rows are the pools. An entry in M is one when a certain individual takes part of a certain pool, and zero otherwise 
% 
% CarrierVec: A vector containing ordered pairs of the number of carriers and their corresponding copy number. For example, [1 3 10 2], correspond to a single carrier with three additional copies (i.e., five copies in total), and ten carriers each having two additional gene copies (i.e, each of these carriers holds five gene copies)
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
% AllowedDelta: Allowed delta from the correct copy number value. For example in case AllowedDelta = -1 then finding the right carrier but assigning a copy number smaller by one is considered correct
% 
% MinFr: minimal value of f to consider (f is the average fraction of occupied droplets, as calibrated by diluting the DNA according to the available number of droplets D)
% 
% StepFr: the delta between consecutive values of f
% 
% MaxFr: the maximum value of f to consider
% 
% MinSuccess: the minimal required success rate. Testing a certain set point would terminate if the rate is lower than this value (e.g. 97%)
% 
% Verbose:enable progress prints (for debug)
% 
% ** Output parameters **
% P - the vector containing the % success of CNV carrier detection for different f values. Detection is considered correct if the correct sample is detected while the copy number may differ according to "AllowedDelta". 
% 
% Pc - Percentage of success in the same format as P. However, in this case successful detection is declared if the correct carrier is found while ignoring the exact copy number (any value Y>0 is considered correct).

function [ P , Pc ] = DIRACC_duplications( M , CarrierVec , D , NumOfTests , LimitNumOfPools , AddDNAPreparationNoise , AddSamplingNoise , AllowedDelta , MinFr , StepFr , MaxFr , MinSuccess , Verbouse  )

NumOfFractions = ceil((MaxFr - MinFr)/StepFr);
P = zeros (NumOfFractions,1);
Pc = zeros(NumOfFractions,1);
index = 1;
[~,NumOfPersons] = size(M);
Tresh = 0.2; % set the treshold for the GPSR algorithm for 0.2 

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
        I = GenerateDuplicateInputVector( NumOfPersons , CarrierVec );
        
        %set the sorted carrier vector 
        [Isort , Iindex] = sort(I);
        Iindex = Iindex(Isort~=0);
        Iindex = sort(Iindex);
        [Hr, Ht] = SimPoolsHits( M , 1 , I , D , Fr , LimitNumOfPools , AddDNAPreparationNoise , AddSamplingNoise );
        [CNV] = ReconstructCnvVec( M , 1 , Hr , Ht , Tresh );
        
        %if the result vector is equal to the
        %input vector then increment the test success counter        
        if CNV == I
            Count = Count + 1;                        
            CarrierCount = CarrierCount + 1;
            
        %check if all of the carriers are identified and if there are in the delta range      
        else      
            
            [SortCNV,CarrierIndex]=sort(CNV);
            CarrierIndex = CarrierIndex(SortCNV > 0);
            CarrierIndex = sort(CarrierIndex);
            if size(CarrierIndex) == size(Iindex) 
                if CarrierIndex == Iindex
                    CarrierCount = CarrierCount + 1;
                    
                    %check if all carriers are in the delta range
                    accuracy = 1;
                    for j = 1:size(CarrierIndex)
                       if AllowedDelta < 0 
                           if CNV(CarrierIndex(j)) < I(CarrierIndex(j)) + AllowedDelta || CNV(CarrierIndex(j)) > I(CarrierIndex(j))
                               accuracy = 0;
                               break;
                           end
                       else % end of AllowedDelta < 0 
                           
                           if CNV(CarrierIndex(j)) > I(CarrierIndex(j)) + AllowedDelta || CNV(CarrierIndex(j)) < I(CarrierIndex(j))
                               accuracy = 0;
                               break;     
                           end                          
                       end % end of AllowedDelta >=0                        
                    end %end of for j = 1:size(CarrierIndex)
                    
                    % if we are in the boundries of delta then increment the
                    % accuracy count
                    if accuracy == 1
                        Count = Count + 1;  
                    end
                    
                end %end of if CarrierIndex == Iindex
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

