clear all;
clc;
disIntMap.is = 'Yes';           
disIntMap.range = [0 20];      

%% Segment integral maps with ROIs

Nexp = length(Masks);
nROI=2;
for j=1:nROI;
    for i=1:Nexp
    masked_SFmap{i,j} = SFmap{1,i}.Map.*Masks{1,i}.RoiMask(:,:,j);
    masked_SFmap{i,j}(masked_SFmap{i,j}==0)=nan;
    end
end


%% Calculate mean value
for j=1:nROI;
    for i=1:Nexp
    Mean_masked_SFmap{i,j} = nanmean(nanmean(masked_SFmap{i,j}));
    end
end


%% Calculate maximum value    
for j=1:nROI;
    for i=1:Nexp
    Max_masked_SFmap{i,j} = nanmax(nanmax(masked_SFmap{i,j}));
    end
end
