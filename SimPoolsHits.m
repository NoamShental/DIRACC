% This function simulates pools of a sensing matrix taking into account both of types of noise
% 1. measurment noise 
% 2. DNA preparation errors
%
% input:
% 
% M - the sensing matrix
% Mode - deletion or duplication CNV : 0 - deletion , 1 - duplication  
% InputVec - The examined persons input vector : In deletion mode : 0 - not carrier , 1 - heterozygous carrier , 2 homozygous carrier
%            In multiplication mode : 0 - no carrier , 1 - additional allele ( 3CNV ) , n ( n+2 ) CNV                 
% D - total number of droplets
% Fr - input sample fraction
% LimitNumOfPools - limit the number of pools to use
% AddDNAPreparationNoise - 0 - No DNA Preparation  noise added , any other value the DNA Preparation noise std
%                         
% AddSamplingNoise - false - no sampling noise , true - add sampling noise 
%
% Output:
%
% Hr - the refernce number of hits for each pool vector
% Ht - the target number of hits for each pool vector

function [ Hr,Ht ] = SimPoolsHits( M , Mode , InputVec , D , Fr , LimitNumOfPools , AddDNAPreparationNoise , AddSamplingNoise )

%set the number of persons and pools 
[MaxNumOfPools , NumOfPersons] = size(M);

%set the num of pools according to the limit
NumOfPools = min( MaxNumOfPools , LimitNumOfPools );

%calculate the number of personess in each pool line 
NumOfPersonsInPool = zeros( NumOfPools , 1 );
for i = 1 : NumOfPools
    NumOfPersonsInPool(i) = sum(M(i,1:NumOfPersons));
end

%calculate the number of droplets occupided by each refernce and target person in average
%in each pool (the ones with no carrier )
Savg = zeros( NumOfPools , 1);
for i = 1 : NumOfPools
    Savg(i) = (Fr * D) / NumOfPersonsInPool(i);
end

%if the AddDNAPreparationNoise is not zero
%calculate an error factor caused due to creation error, it is normaly
%distrobuted with std of 0.1
TargetCreationErrorFactor = ones( NumOfPersons , 1 );
ReferenceCreationErrorFactor = ones( NumOfPersons , 1 );

Sr = zeros(NumOfPools,1);
St = zeros(NumOfPools,1);

% ......... calculate the expected hit rates in the pools
%           on the tagert test we shall add 1/2 of expected amount of droplets if the person is a carrier
        
for i = 1: NumOfPools
   
    % recalcualting the creation noise for the next pool 
    if AddDNAPreparationNoise ~= 0
        for k=1:NumOfPersons
            TargetCreationErrorFactor(k) = normrnd(1,AddDNAPreparationNoise);
        end
        ReferenceCreationErrorFactor = TargetCreationErrorFactor;    
    end
    
    %if the mode is deletion then:
    %for each person in the pool add to 
    %the refernce : the average multiplied by the creation error factor
    %the target : if normal person : the average multiplied by the creation error factor 
    %             if a Homozygous carrier then half the amount 
    %             if a Heterozygous nothing is added
    if Mode == 0
        for j = 1 : NumOfPersons
            if M(i,j) == 1
           
                Sr(i) = Sr(i) + (Savg(i)*ReferenceCreationErrorFactor(j));
           
                %if a Homozygous carrier person add only half of the average Hits
                %amount
                if InputVec(j) == 1
                    St(i) = St(i) + (Savg(i)*TargetCreationErrorFactor(j))/2;
                end
           
                %if a non carirer person add tha average hits
                if InputVec(j) == 0
                    St(i) = St(i) + (Savg(i)*TargetCreationErrorFactor(j));
                end
           
                %if a Heterozygous carrier we shall not add enything .....
           
            end % M(i,j) == 1
        end % for j = 1 : NumOfPerson
    
    %if the mode is multiplication then:
    %for each person in the pool add to 
    %the refernce : the average multiplied by the creation error factor
    %the target : if normal person InputVec(j) = 0 : the average multiplied by the creation error factor 
    %             if the input is not normal InputVec(j) > 0 then 
    %             add an additional half of the average multiplied by the creation error factor multiplied by InputVec(j)     
    else
        for j = 1 : NumOfPersons
            if M(i,j) == 1
           
                Sr(i) = Sr(i) + (Savg(i)*ReferenceCreationErrorFactor(j));
                St(i) = St(i) + (InputVec(j)+2)*((Savg(i)*TargetCreationErrorFactor(j))/2);
                     
            end % M(i,j) == 1
        end % for j = 1 : NumOfPerson        
    end % Mode
    
    % verification we dont exceed D
    if Sr(i) > D
        Sr(i) = D;
    end
    if St(i) > D
        St(i) = D;
    end
end %for i = 1: NumOfPools

if AddSamplingNoise == true
    
    % ........ handle the fraction (reading) noise error on deletion mode ( mode 0 )
    % step 1 run a multinomial distribution for 4 posibilties on each droplet
    % p1 - at least one molecule containing target and refernce alleles and one molecule containing only refernce alleles in the droplet 
    %      p1 =1-(1-1/D)^St(i)-(1-1/D)^(Sr(i)-St(i))+(1-1/D)^Sr(i)  
    % p2 - only molecules containg target and refernce alleles in the droplet 
    %      p2 =(1-1/D)^(Sr(i)-St(i))-(1-1/D)^Sr(i)         
    % p3 - only molecules containg only refernce alleles in the droplet 
    %      p3 =(1-1/D)^St(i)-(1-1/D)^Sr(i)         
    % p4 - no molecule in the droplet 
    %      p4 =(1-1/D)^Sr(i)             
    % [h1,h2,h3,h4 ] ~ MultiNomial(D,p1,p2,p3,p4 )        
    %
    % step 2 finalize the refernce and target hits
    % Hr(i) = h1 + h2 + h3  - the sum of all droplets with any molecule
    % Ht(i) = h1 + h2  - the sum pf droplets only with molecules that contain
    % also the target
    if Mode == 0
        Ht = zeros(NumOfPools,1);
        Hr = zeros(NumOfPools,1);
        p = zeros (1,4);
        for i = 1 : NumOfPools
            p(1) =1-(1-1/D)^St(i)-(1-1/D)^(Sr(i)-St(i))+(1-1/D)^Sr(i);              
            p(2) =(1-1/D)^(Sr(i)-St(i))-(1-1/D)^Sr(i);
            p(3) =(1-1/D)^St(i)-(1-1/D)^Sr(i);
            p(4) =(1-1/D)^Sr(i);        
            %make sure p is not below 0 - due calculation saturation error 
            p(1) = max(0,p(1));
            p(2) = max(0,p(2));
            p(3) = max(0,p(3));
            p(4) = max(0,p(4));
        
            h = mnrnd(D,p);    
            Hr(i) = h(1)+h(2)+h(3);
            Ht(i) = h(1)+h(2);
        end
        
    % ........ handle the fraction (reading) noise error on multiplication mode ( mode 1 )
    % step 1 run a multinomial distribution for 4 posibilties on each droplet
    % p1 - at least one molecule containing target allele and one molecule containing only target alleles in the droplet 
    %      p1 =1-(1-1/D)^Sr(i)-(1-1/D)^St(i)+(1-1/D)^(Sr(i)+St(i))  
    % p2 - only molecules containg target and refernce alleles in the droplet 
    %      p2 =(1-1/D)^St(i)-(1-1/D)^(Sr(i)+St(i))         
    % p3 - only molecules containg only target alleles in the droplet 
    %      p3 =(1-1/D)^Sr(i)-(1-1/D)^(Sr(i)+St(i))         
    % p4 - no molecule in the droplet 
    %      p4 =(1-1/D)^(Sr(i)+St(i))             
    % [h1,h2,h3,h4 ] ~ MultiNominal(D,p1,p2,p3,p4 )        
    %
    % step 2 finalize the refernce and target hits
    % Hr(i) = h1 + h2  - the sum of all droplets with that contains the reference
    % Ht(i) = h1 + h3  - the sum pf droplets only with molecules that contais the target
    else %Mode == 1
        Ht = zeros(NumOfPools,1);
        Hr = zeros(NumOfPools,1);
        p = zeros (1,4);
        for i = 1 : NumOfPools
            p(1) =1-(1-1/D)^Sr(i)-(1-1/D)^St(i)+(1-1/D)^(Sr(i)+St(i));              
            p(2) =(1-1/D)^St(i)-(1-1/D)^(Sr(i)+St(i));
            p(3) =(1-1/D)^Sr(i)-(1-1/D)^(Sr(i)+St(i));
            p(4) =(1-1/D)^(Sr(i)+St(i));        
            %make sure p is not below 0 - due calculation saturation error 
            p(1) = max(0,p(1));
            p(2) = max(0,p(2));
            p(3) = max(0,p(3));
            p(4) = max(0,p(4));
        
            h = mnrnd(D,p);    
            Hr(i) = h(1)+h(2);
            Ht(i) = h(1)+h(3);
        end
    end %if Mode == 0
    
% ....... no noise case     
else
    Hr = Sr;
    Ht = St;
end

end % function [ Hr,Ht ] = SimPoolsHits( M , InputVec , D , Fr , LimitNumOfPools , AddDNAPreparationNoise , CreationNoiseMode , AddMeasurementNoise )

