% this function generates an input vector according to the number of
% mutliplied carriers

%input:
%NumOfPersons - the number of persons in the test
%CarrierVec - a vector containing: [ NumOfCarrier1 , CNV1 , NumOfCarrier2 , CNV2 , .... ]
%             NumOfCarrierk - the number of carriers of the spesifed CNV(k)
%             CNVk - the CNV of the carrier ( 0 - normal , 1 aditional
%             allele ,2 aditional 2 alleles , k - aditional k alleles )

function [ I ] = GenerateDuplicateInputVector( NumOfPersons , CarrierVec )

% ...... prpeare the input vector 
while 1

    %init the input vector 
    I = zeros(NumOfPersons,1);
    
    %main loop for randomizing the positions of all carriers according to
    %there types
    [~,NumOfCarrierTypes] = size(CarrierVec);
    for j = 1:2:NumOfCarrierTypes
   
        %randomize the current CNV carriers and set them in a vector
        CNV = randi(NumOfPersons,1,CarrierVec(j));

        %load the CNV carriers to the vector
        AllocationSuccess = 1;
        for i = 1:CarrierVec(j)
        
            %if the carrier is captured then we failed and we need to try to
            %randomize again
            if I(CNV(i)) ~= 0
                AllocationSuccess = 0;
                break;
            
            % set the Heterozygous carrier to the vector     
            else
                I(CNV(i)) = CarrierVec(j+1);
            end
        
         end %end of for i = 1:CarrierVec(j)   
    
         % if we have failed then break
         if AllocationSuccess == 0
            break;
         end
         
    end % end of j = 1:2:NumOfCarrierTypes 
    
    % if we have failed then continue
    if AllocationSuccess == 0
        continue;
    end
    
    % we succeeded break the loop
    break;
    
end % end of while 1 

end % end of function [ I ] = GenerateMultiInputVector( NumOfPersons , CarrierVec )

