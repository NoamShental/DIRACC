%this function reconstructs the CNV vector using the refernce and target
%hits measurments in the pools.

% input:
% 
% M - the Bernoulli sensing matrix
% Mode - deletion or duplication CNV : 0 - deletion , 1 - duplication 
% Hr - the refernce number of hits for each pool vector
% Ht - the target number of hits for each pool vector
% Tresh - the treshold of filtering noise from the GPSR output 
%
% The function supports up to 20 carriers
%
% Output:
%
% CNV - The reconstructed CNV for each person vector
function [ CNV ] = ReconstructCnvVec( M , Mode , Hr , Ht , Tresh )

%set the number of persons and pools 
[~ , NumOfPersons] = size(M);

%set the num of pools according to the limit
[NumOfPools,~] = size(Hr);

% ........ calculate the ratio for all measurments 
Ratio = zeros(NumOfPools,1);
for i = 1 : NumOfPools
    Ratio(i) = Ht(i)/Hr(i);        
end

% ......... perform reconstruction
Mr = M(1:NumOfPools,1:NumOfPersons);

%Normilize the vector to the number of missing alleles (CNV's)
%and multiple it by the number of alleles in a pool (2 per person)
if Mode == 1
    PoolCNV = Ratio - 1;
else
    PoolCNV = abs(1 - Ratio);
end

for i=1:NumOfPools
    [~,NumPfPersonInPool]=size(find(M(i,:)));
    PoolCNV(i) = 2* PoolCNV(i) * NumPfPersonInPool;
end

%clip the PoolCNV to reduce measurment noise
for i=1:NumOfPools
    if abs(PoolCNV(i)) < Tresh
       PoolCNV(i) = 0;
    end
end
        
% fractional solution
tau = 0.005*max(abs(Mr'*PoolCNV));
fractionalOutput = applyGPSR(PoolCNV,Mr,tau);

% ........ analyze the fractional output and find the CNV

% we shall eliminate all values close to zero from the fractional output 
% we shall solve the minimal norm ||Mx-PoolCNV||_2 directly using the cvx library 

%1.sort the fractional solution
[OutFracSort,Index] = sort(fractionalOutput,'descend');

%2.clear noisy values from the vector - assume each absolute value below treshold is zero
Index = Index(abs(OutFracSort) > Tresh);
[NumOfIndexes,~] = size(Index);

%if the vector is empty then exit with CNV all zeros
if NumOfIndexes == 0
    CNV = zeros(NumOfPersons,1);
    return;
end
  
% now we shall use only the colomns of Mr that are choosen as
% candidates
Mc = Mr(:,Index);
    
% find the minimal norm2 for the error cost : ||Mc*X-PoolCNV||_2 
cvx_begin quiet
    variable X(NumOfIndexes);
    minimize (norm(Mc*X-PoolCNV,2));
cvx_end
    
%set the CNV vector
CNV = zeros(NumOfPersons,1);
for i=1:NumOfIndexes
    CNV(Index(i)) = floor (X(i)+0.5);
end

end % function [ CNV , Success ] = ReconstructCnvVec( M , Hr , Ht )

