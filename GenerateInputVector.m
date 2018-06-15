% this function generates an input vector according to the number of
% heterozygous and homozygous carriers

%input:
%NumOfPersons - the number of persons in the test
%NumHeterozygousCarriers - the number of Heterozygous in the Subject groups
%NumHomozygousCarriers - the number of Homozygous in the Subject groups

function [ I ] = GenerateInputVector( NumOfPersons , NumHeterozygousCarriers , NumHomozygousCarriers )

% ...... prpeare the input vector 
while 1

    %init the input vector 
    I = zeros(NumOfPersons,1);
    
    %randomize the Heterozygous Carriers carriers and set them in a vector
    HeterozygousCarriers = randi(NumOfPersons,1,NumHeterozygousCarriers);
    
    %load the Heterozygous carriers to the vector
    AllocationSuccess = 1;
    for i = 1:NumHeterozygousCarriers
        
        %if the carrier is captured then we failed and we need to try to
        %randomize again
        if I(HeterozygousCarriers(i)) ~= 0
            AllocationSuccess = 0;
            break;
            
        % set the Heterozygous carrier to the vector     
        else
            I(HeterozygousCarriers(i)) = 1;
        end
        
    end %end of for i = 1:NumHeterozygousCarriers   
    
    % if we have failed then continue
    if AllocationSuccess == 0
        continue;
    end

    %randomize the Homozygous Carriers carriers and set them in a vector
    HomozygousCarriers = randi(NumOfPersons,1,NumHomozygousCarriers);
    
    %load the Homozygous carriers to the vector
    AllocationSuccess = 1;
    for i = 1:NumHomozygousCarriers
        
        %if the carrier is captured then we failed and we need to try to
        %randomize again
        if I(HomozygousCarriers(i)) ~= 0
            AllocationSuccess = 0;
            break;
            
        % set the Homozygous carrier to the vector     
        else
            I(HomozygousCarriers(i)) = 2;
        end
        
    end %end of for i = 1:NumHeterozygousCarriers   
    
    % if we have failed then continue
    if AllocationSuccess == 0
        continue;
    end
    
    % we succeeded break the loop
    break;
end % end of while Success == 0 

end

