%%
% Organise_data script, Maria Yanez Lopez, 
% This version reads all the .mat files coming from
% LCmodel2_matlab_version2, which have been renamed 1.mat, 2.mat, etc
% Input: .mat files (1:15), from LCmodel2_matlab_version2
% It reads FWHM, FWHM_Hz and SN (all 15x2 double)
% It organises the data from the 11 reliable metabolites (GABA, Gln, Glu,Ins,Tau,NAA,tNAA,tCr,tCho,Glx, ML9) in three ways:
%       abs_metabolite (15x2 double) = absolute value
%       ratio_metabolite (15x2 double) = ratio to tCr
%       ratio_sum_metabolite (15x2 double) = ratio to sum of 11 reliable
%       metabolites

% OUTPUT: 
        % Multiple graphs plus Ins_data.mat containing the Ins data
        
% NOTE: This script has now been extended to do the same for all 11 selected metabolites.
% New outputs are: 
% Glu_data.mat 
% Glx_data.mat 
% NAA_data.mat 
% tNAA_data.mat 
% tCho_data.mat 
% tCr_data.mat 
% Tau_data.mat
% Gln_data.mat 
% GABA_data.mat 
% ML9_data.mat 

clear all;
clc;
% set(0,'DefaultFigureVisible','off');  % all subsequent figures "off"

%Load the whole dataset
for i = 1:15
    iname = num2str(i);
    filename = strcat(iname,'.mat');
    Mouse(i) = load(filename);
    
end
timepoints = 2;
metabolites = 11;
clear('i','iname','filename');
%Predefine variables
% Scan_ID = char(zeros(15,1));
FWHM = zeros((size(Mouse,2)),timepoints);
FWHM_Hz = zeros((size(Mouse,2)),timepoints);
SN = zeros((size(Mouse,2)),timepoints);
PBS_ratio = zeros((size(Mouse,2)),metabolites);
LPS_ratio = zeros((size(Mouse,2)),metabolites);

PBS = zeros((size(Mouse,2)),metabolites);
LPS = zeros((size(Mouse,2)),metabolites);


% Metabolites = char('Glu','Ins','NAA','Tau','tCho','tNAA','tCr','Glx');
for i = [1:15]
    for j=1:timepoints
        FWHM(i,j) = Mouse(1,i).Whole_dataset{1,j}.FWHM;
        SN(i,j) = Mouse(1,i).Whole_dataset{1,j}.SN;
        Scan_ID{i,j} = Mouse(1,i).Whole_dataset{1,j}.Filename;
    end
end
for i = 1:size(Mouse,2)
    
    for k=1:metabolites
        PBS_ratio(i,k) =  Mouse(1,i).Reliable_dataset{1,1}.Concentrations(k,2);
        LPS_ratio(i,k) =  Mouse(1,i).Reliable_dataset{1,2}.Concentrations(k,2);
        PBS(i,k) =  Mouse(1,i).Reliable_dataset{1,1}.Concentrations(k,1);
        LPS(i,k) =  Mouse(1,i).Reliable_dataset{1,2}.Concentrations(k,1);
       
    end
end


FWHM_Hz = 400.*FWHM;

abs = cat(3,PBS,LPS);
ratios = cat(3,PBS_ratio,LPS_ratio);
clear('PBS','LPS','PBS_ratio','LPS_ratio','i','j','k');

% Metabolites = char('Glu','Ins','NAA','Tau','tCho','tNAA','tCr','Glx');
%Predefine variables 
ratios_Glu = zeros((size(Mouse,2)),timepoints);
ratios_Ins = zeros((size(Mouse,2)),timepoints);
ratios_NAA = zeros((size(Mouse,2)),timepoints);
ratios_Tau = zeros((size(Mouse,2)),timepoints);
ratios_tCho = zeros((size(Mouse,2)),timepoints);
ratios_tNAA = zeros((size(Mouse,2)),timepoints);
ratios_tCr = zeros((size(Mouse,2)),timepoints);
ratios_Glx = zeros((size(Mouse,2)),timepoints);
ratios_GABA = zeros((size(Mouse,2)),timepoints);
ratios_tNAA = zeros((size(Mouse,2)),timepoints);
ratios_ML9 = zeros((size(Mouse,2)),timepoints);

ratios_sum_Glu = zeros((size(Mouse,2)),timepoints);
ratios_sum_Ins = zeros((size(Mouse,2)),timepoints);
ratios_sum_NAA = zeros((size(Mouse,2)),timepoints);
ratios_sum_Tau = zeros((size(Mouse,2)),timepoints);
ratios_sum_tCho = zeros((size(Mouse,2)),timepoints);
ratios_sum_tNAA = zeros((size(Mouse,2)),timepoints);
ratios_sum_tCr = zeros((size(Mouse,2)),timepoints);
ratios_sum_Glx = zeros((size(Mouse,2)),timepoints);
ratios_sum_GABA = zeros((size(Mouse,2)),timepoints);
ratios_sum_tNAA = zeros((size(Mouse,2)),timepoints);
ratios_sum_ML9 = zeros((size(Mouse,2)),timepoints);

abs_Glu = zeros((size(Mouse,2)),timepoints);
abs_Ins = zeros((size(Mouse,2)),timepoints);
abs_NAA = zeros((size(Mouse,2)),timepoints);
abs_Tau = zeros((size(Mouse,2)),timepoints);
abs_tCho = zeros((size(Mouse,2)),timepoints);
abs_tNAA = zeros((size(Mouse,2)),timepoints);
abs_tCr = zeros((size(Mouse,2)),timepoints);
abs_Glx = zeros((size(Mouse,2)),timepoints);
abs_GABA = zeros((size(Mouse,2)),timepoints);
abs_tNAA = zeros((size(Mouse,2)),timepoints);
abs_ML9 = zeros((size(Mouse,2)),timepoints);

for k=1:(size(Mouse,2))
    for h=1:timepoints
        ratios_Glu(k,h) = ratios(k,3,h);
        ratios_Ins(k,h) = ratios(k,4,h);
        ratios_NAA(k,h) = ratios(k,5,h);
        ratios_Tau(k,h) = ratios(k,6,h);
        ratios_tCho(k,h) = ratios(k,7,h);
        ratios_tNAA(k,h) = ratios(k,8,h);
        ratios_tCr(k,h) = ratios(k,9,h);
        ratios_Glx(k,h) = ratios(k,10,h);
        ratios_GABA(k,h) = ratios(k,1,h);
        ratios_Gln(k,h) = ratios(k,2,h);
        ratios_ML9(k,h) = ratios(k,11,h);
        
        abs_Glu(k,h) = abs(k,3,h);
        abs_Ins(k,h) = abs(k,4,h);
        abs_NAA(k,h) = abs(k,5,h);
        abs_Tau(k,h) = abs(k,6,h);
        abs_tCho(k,h) = abs(k,7,h);
        abs_tNAA(k,h) = abs(k,8,h);
        abs_tCr(k,h) = abs(k,9,h);
        abs_Glx(k,h) = abs(k,10,h);
        abs_GABA(k,h) = abs(k,1,h);
        abs_Gln(k,h) = abs(k,2,h);
        abs_ML9(k,h) = abs(k,11,h);
    end  
end
clear('h','k','ratios','abs');

for k=1:(size(Mouse,2))
    for h=1:timepoints
        ratios_sum_Glu(k,h) = abs_Glu(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_Ins(k,h) = abs_Ins(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_NAA(k,h) = abs_NAA(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_Tau(k,h) = abs_Tau(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_tCho(k,h) = abs_tCho(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_tNAA(k,h) = abs_tNAA(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_tCr(k,h) = abs_tCr(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_Glx(k,h) = abs_Glx(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_tNAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_GABA(k,h) = abs_GABA(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_NAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_Gln(k,h) = abs_Gln(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_NAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));
        ratios_sum_ML9(k,h) = abs_ML9(k,h)./(abs_Glx(k,h) + abs_tCr(k,h) + abs_NAA(k,h) + abs_tCho(k,h) + abs_Tau(k,h) + abs_Ins(k,h));

    end  
end

clear('h','k', 'Mouse');


%% Ins plots

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_Ins(k,h) == 0;
            ratios_Ins(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_Ins(k,h) == 0;
            ratios_sum_Ins(k,h) = nan;
        end
    end
end


%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};

ratios_Ins_AD = nan(15,2);
ratios_sum_Ins_AD = nan(15,2);
ratios_Ins_WT = nan(15,2);
ratios_sum_Ins_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_Ins_AD(k,:) = ratios_Ins(k,:);
        ratios_sum_Ins_AD(k,:) = ratios_sum_Ins(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_Ins_WT(k,:) = ratios_Ins(k,:);
        ratios_sum_Ins_WT(k,:) = ratios_sum_Ins(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end



%% Save Ins data
save('Ins_data.mat','ratios_Ins*','ratios_sum_Ins*','abs_Ins*');

%% Glu

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_Glu(k,h) == 0;
            ratios_Glu(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_Glu(k,h) == 0;
            ratios_sum_Glu(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};


ratios_Glu_AD = nan(15,2);
ratios_sum_Glu_AD = nan(15,2);
ratios_Glu_WT = nan(15,2);
ratios_sum_Glu_WT = nan(15,2);



for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_Glu_AD(k,:) = ratios_Glu(k,:);
        ratios_sum_Glu_AD(k,:) = ratios_sum_Glu(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_Glu_WT(k,:) = ratios_Glu(k,:);
        ratios_sum_Glu_WT(k,:) = ratios_sum_Glu(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

%% Save Glu data
save('Glu_data.mat','ratios_Glu*','ratios_sum_Glu*','abs_Glu*');

%% Glx

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_Glx(k,h) == 0;
            ratios_Glx(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_Glx(k,h) == 0;
            ratios_sum_Glx(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};


ratios_Glx_AD = nan(15,2);
ratios_sum_Glx_AD = nan(15,2);
ratios_Glx_WT = nan(15,2);
ratios_sum_Glx_WT = nan(15,2);


for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_Glx_AD(k,:) = ratios_Glx(k,:);
        ratios_sum_Glx_AD(k,:) = ratios_sum_Glx(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_Glx_WT(k,:) = ratios_Glx(k,:);
        ratios_sum_Glx_WT(k,:) = ratios_sum_Glx(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end


%% Save Glx data
save('Glx_data.mat','ratios_Glx*','ratios_sum_Glx*','abs_Glx*');

%% Tau

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_Tau(k,h) == 0;
            ratios_Tau(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_Tau(k,h) == 0;
            ratios_sum_Tau(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};

ratios_Tau_AD = nan(15,2);
ratios_sum_Tau_AD = nan(15,2);
ratios_Tau_WT = nan(15,2);
ratios_sum_Tau_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_Tau_AD(k,:) = ratios_Tau(k,:);
        ratios_sum_Tau_AD(k,:) = ratios_sum_Tau(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_Tau_WT(k,:) = ratios_Tau(k,:);
        ratios_sum_Tau_WT(k,:) = ratios_sum_Tau(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

%% Save Tau data
save('Tau_data.mat','ratios_Tau*','ratios_sum_Tau*','abs_Tau*');

%% NAA

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_NAA(k,h) == 0;
            ratios_NAA(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_NAA(k,h) == 0;
            ratios_sum_NAA(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};

ratios_NAA_AD = nan(15,2);
ratios_sum_NAA_AD = nan(15,2);
ratios_NAA_WT = nan(15,2);
ratios_sum_NAA_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_NAA_AD(k,:) = ratios_NAA(k,:);
        ratios_sum_NAA_AD(k,:) = ratios_sum_NAA(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_NAA_WT(k,:) = ratios_NAA(k,:);
        ratios_sum_NAA_WT(k,:) = ratios_sum_NAA(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

%% Save NAA data
save('NAA_data.mat','ratios_NAA*','ratios_sum_NAA*','abs_NAA*');

%% tNAA

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_tNAA(k,h) == 0;
            ratios_tNAA(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_tNAA(k,h) == 0;
            ratios_sum_tNAA(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};

ratios_tNAA_AD = nan(15,2);
ratios_sum_tNAA_AD = nan(15,2);
ratios_tNAA_WT = nan(15,2);
ratios_sum_tNAA_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_tNAA_AD(k,:) = ratios_tNAA(k,:);
        ratios_sum_tNAA_AD(k,:) = ratios_sum_tNAA(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_tNAA_WT(k,:) = ratios_tNAA(k,:);
        ratios_sum_tNAA_WT(k,:) = ratios_sum_tNAA(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

%% Save tNAA data
save('tNAA_data.mat','ratios_tNAA*','ratios_sum_tNAA*','abs_tNAA*');

%% tCr

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_tCr(k,h) == 0;
            ratios_tCr(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_tCr(k,h) == 0;
            ratios_sum_tCr(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};

ratios_tCr_AD = nan(15,2);
ratios_sum_tCr_AD = nan(15,2);
ratios_tCr_WT = nan(15,2);
ratios_sum_tCr_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_tCr_AD(k,:) = ratios_tCr(k,:);
        ratios_sum_tCr_AD(k,:) = ratios_sum_tCr(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_tCr_WT(k,:) = ratios_tCr(k,:);
        ratios_sum_tCr_WT(k,:) = ratios_sum_tCr(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

%% Save tCr data
save('tCr_data.mat','ratios_tCr*','ratios_sum_tCr*','abs_tCr*');

%% tCho

%Ratios to creatine
for k=1:15
    for h=1:2 
        if ratios_tCho(k,h) == 0;
            ratios_tCho(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_tCho(k,h) == 0;
            ratios_sum_tCho(k,h) = nan;
        end
    end
end

%Differentiating type
Type = {'WT','AD','WT','AD','WT','WT','WT','WT', 'AD', 'AD', 'WT', 'WT', 'AD', 'AD', 'WT'};

ratios_tCho_AD = nan(15,2);
ratios_sum_tCho_AD = nan(15,2);
ratios_tCho_WT = nan(15,2);
ratios_sum_tCho_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_tCho_AD(k,:) = ratios_tCho(k,:);
        ratios_sum_tCho_AD(k,:) = ratios_sum_tCho(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_tCho_WT(k,:) = ratios_tCho(k,:);
        ratios_sum_tCho_WT(k,:) = ratios_sum_tCho(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

%% Save tCho data
save('tCho_data.mat','ratios_tCho*','ratios_sum_tCho*','abs_tCho*');
%% Gln

%Ratios to creatine
for k=1:15
    for h=1:2
        if ratios_Gln(k,h) == 0;
            ratios_Gln(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_Gln(k,h) == 0;
            ratios_sum_Gln(k,h) = nan;
        end
    end
end


ratios_Gln_LPS = nan(15,2);
ratios_sum_Gln_LPS = nan(15,2);
ratios_Gln_PBS = nan(15,2);
ratios_sum_Gln_PBS = nan(15,2);

ratios_Gln_AD = nan(15,2);
ratios_sum_Gln_AD = nan(15,2);
ratios_Gln_WT = nan(15,2);
ratios_sum_Gln_WT = nan(15,2);

ratios_Gln_LPS_AD = nan(15,2);
ratios_Gln_PBS_AD = nan(15,2);
ratios_Gln_LPS_WT = nan(15,2);
ratios_Gln_PBS_WT = nan(15,2);

ratios_sum_Gln_LPS_AD = nan(15,2);
ratios_sum_Gln_PBS_AD = nan(15,2);
ratios_sum_Gln_LPS_WT = nan(15,2);
ratios_sum_Gln_PBS_WT = nan(15,2);


for k=1:15
    if char(Type(1,k)) == 'AD'
        ratios_Gln_AD(k,:) = ratios_Gln(k,:);
        ratios_sum_Gln_AD(k,:) = ratios_sum_Gln(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_Gln_WT(k,:) = ratios_Gln(k,:);
        ratios_sum_Gln_WT(k,:) = ratios_sum_Gln(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

% Save Gln data
save('Gln_data.mat','ratios_Gln*','ratios_sum_Gln*','abs_Gln*');
%% GABA

%Ratios to creatine
for k=1:15
    for h=1:2
        if ratios_GABA(k,h) == 0;
            ratios_GABA(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2
        if ratios_sum_GABA(k,h) == 0;
            ratios_sum_GABA(k,h) = nan;
        end
    end
end


ratios_GABA_LPS = nan(15,2);
ratios_sum_GABA_LPS = nan(15,2);
ratios_GABA_PBS = nan(15,2);
ratios_sum_GABA_PBS = nan(15,2);

ratios_GABA_AD = nan(15,2);
ratios_sum_GABA_AD = nan(15,2);
ratios_GABA_WT = nan(15,2);
ratios_sum_GABA_WT = nan(15,2);

ratios_GABA_LPS_AD = nan(15,2);
ratios_GABA_PBS_AD = nan(15,2);
ratios_GABA_LPS_WT = nan(15,2);
ratios_GABA_PBS_WT = nan(15,2);

ratios_sum_GABA_LPS_AD = nan(15,2);
ratios_sum_GABA_PBS_AD = nan(15,2);
ratios_sum_GABA_LPS_WT = nan(15,2);
ratios_sum_GABA_PBS_WT = nan(15,2);


for k=1:15
    if char(Type(1,k)) == 'AD'
        ratios_GABA_AD(k,:) = ratios_GABA(k,:);
        ratios_sum_GABA_AD(k,:) = ratios_sum_GABA(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_GABA_WT(k,:) = ratios_GABA(k,:);
        ratios_sum_GABA_WT(k,:) = ratios_sum_GABA(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end


% Save GABA data
save('GABA_data.mat','ratios_GABA*','ratios_sum_GABA*','abs_GABA*');
%% ML9

%Ratios to creatine
for k=1:15
    for h=1:2
        if ratios_ML9(k,h) == 0;
            ratios_ML9(k,h) = nan;
        end
    end
end

%Ratios to sum of metabolites
for k=1:15
    for h=1:2 
        if ratios_sum_ML9(k,h) == 0;
            ratios_sum_ML9(k,h) = nan;
        end
    end
end


ratios_ML9_LPS = nan(15,2);
ratios_sum_ML9_LPS = nan(15,2);
ratios_ML9_PBS = nan(15,2);
ratios_sum_ML9_PBS = nan(15,2);

ratios_ML9_AD = nan(15,2);
ratios_sum_ML9_AD = nan(15,2);
ratios_ML9_WT = nan(15,2);
ratios_sum_ML9_WT = nan(15,2);

ratios_ML9_LPS_AD = nan(15,2);
ratios_ML9_PBS_AD = nan(15,2);
ratios_ML9_LPS_WT = nan(15,2);
ratios_ML9_PBS_WT = nan(15,2);

ratios_sum_ML9_LPS_AD = nan(15,2);
ratios_sum_ML9_PBS_AD = nan(15,2);
ratios_sum_ML9_LPS_WT = nan(15,2);
ratios_sum_ML9_PBS_WT = nan(15,2);

for k=1:15 
    if char(Type(1,k)) == 'AD'
        ratios_ML9_AD(k,:) = ratios_ML9(k,:);
        ratios_sum_ML9_AD(k,:) = ratios_sum_ML9(k,:);
    elseif char(Type(1,k)) == 'WT'  
        ratios_ML9_WT(k,:) = ratios_ML9(k,:);
        ratios_sum_ML9_WT(k,:) = ratios_sum_ML9(k,:);
    else
        printf('Warning, the type factor has not been understood');
    end
end

% Save ML9 data
save('ML9_data.mat','ratios_ML9*','ratios_sum_ML9*','abs_ML9*');