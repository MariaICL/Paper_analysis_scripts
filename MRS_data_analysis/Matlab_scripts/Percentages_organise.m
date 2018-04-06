%%
% Percentages_organise script, Maria Yanez Lopez, 
% This version reads %_data.mat coming from Percentages
% Input: %_data.mat coming from Percentages
% It organises the data in groups:
%     Perc_ratios_(sum_)metabolite
%     Perc_ratios_(sum_)metabolite_AD
%     Perc_ratios_(sum_)metabolite_WT
%     Perc_ratios_(sum_)metabolite_WT
%     Perc_ratios_(sum_)metabolite_LPS
%     Perc_ratios_(sum_)metabolite_AD
%     Perc_ratios_(sum_)metabolite_WT
%     Perc_ratios_(sum_)metabolite_WT
%     Perc_ratios_(sum_)metabolite_PBS
%     Perc_ratios_(sum_)metabolite_AD
%     Perc_ratios_(sum_)metabolite_WT
%     Perc_ratios_(sum_)metabolite_WT
% OUTPUT: 
        % %_data.mat containing the metabolite data


clear all;
clc;

%Load data
load('%_data.mat');
load('Type.mat');


%% Ins

%Defining variables
Perc_ratios_Ins_AD = nan(15,1);
Perc_ratios_Ins_WT = nan(15,1);

Perc_ratios_sum_Ins_AD = nan(15,1);
Perc_ratios_sum_Ins_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_Ins_AD(k,:) = Perc_ratios_Ins(k,:);
        Perc_ratios_sum_Ins_AD(k,:) = Perc_ratios_sum_Ins(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Ins_WT(k,:) = Perc_ratios_Ins(k,:);
        Perc_ratios_sum_Ins_WT(k,:) = Perc_ratios_sum_Ins(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_Ins_AD(k,:) = Perc_ratios_Ins(k,:);
        Perc_ratios_sum_Ins_AD(k,:) = Perc_ratios_sum_Ins(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Ins_WT(k,:) = Perc_ratios_Ins(k,:);
        Perc_ratios_sum_Ins_WT(k,:) = Perc_ratios_sum_Ins(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Ins_WT(k,:) = Perc_ratios_Ins(k,:);
        Perc_ratios_sum_Ins_WT(k,:) = Perc_ratios_sum_Ins(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_Ins_WT(k,:) = Perc_ratios_Ins(k,:);
        Perc_ratios_sum_Ins_WT(k,:) = Perc_ratios_sum_Ins(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end

%% Glu
%Defining variables
Perc_ratios_Glu_AD = nan(15,1);
Perc_ratios_Glu_WT = nan(15,1);

Perc_ratios_sum_Glu_AD = nan(15,1);
Perc_ratios_sum_Glu_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_Glu_AD(k,:) = Perc_ratios_Glu(k,:);
        Perc_ratios_sum_Glu_AD(k,:) = Perc_ratios_sum_Glu(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Glu_WT(k,:) = Perc_ratios_Glu(k,:);
        Perc_ratios_sum_Glu_WT(k,:) = Perc_ratios_sum_Glu(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_Glu_AD(k,:) = Perc_ratios_Glu(k,:);
        Perc_ratios_sum_Glu_AD(k,:) = Perc_ratios_sum_Glu(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Glu_WT(k,:) = Perc_ratios_Glu(k,:);
        Perc_ratios_sum_Glu_WT(k,:) = Perc_ratios_sum_Glu(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Glu_WT(k,:) = Perc_ratios_Glu(k,:);
        Perc_ratios_sum_Glu_WT(k,:) = Perc_ratios_sum_Glu(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_Glu_WT(k,:) = Perc_ratios_Glu(k,:);
        Perc_ratios_sum_Glu_WT(k,:) = Perc_ratios_sum_Glu(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% Glx
%Defining variables
Perc_ratios_Glx_AD = nan(15,1);
Perc_ratios_Glx_WT = nan(15,1);

Perc_ratios_sum_Glx_AD = nan(15,1);
Perc_ratios_sum_Glx_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_Glx_AD(k,:) = Perc_ratios_Glx(k,:);
        Perc_ratios_sum_Glx_AD(k,:) = Perc_ratios_sum_Glx(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Glx_WT(k,:) = Perc_ratios_Glx(k,:);
        Perc_ratios_sum_Glx_WT(k,:) = Perc_ratios_sum_Glx(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_Glx_AD(k,:) = Perc_ratios_Glx(k,:);
        Perc_ratios_sum_Glx_AD(k,:) = Perc_ratios_sum_Glx(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Glx_WT(k,:) = Perc_ratios_Glx(k,:);
        Perc_ratios_sum_Glx_WT(k,:) = Perc_ratios_sum_Glx(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Glx_WT(k,:) = Perc_ratios_Glx(k,:);
        Perc_ratios_sum_Glx_WT(k,:) = Perc_ratios_sum_Glx(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_Glx_WT(k,:) = Perc_ratios_Glx(k,:);
        Perc_ratios_sum_Glx_WT(k,:) = Perc_ratios_sum_Glx(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% Tau
%Defining variables
Perc_ratios_Tau_AD = nan(15,1);
Perc_ratios_Tau_WT = nan(15,1);

Perc_ratios_sum_Tau_AD = nan(15,1);
Perc_ratios_sum_Tau_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_Tau_AD(k,:) = Perc_ratios_Tau(k,:);
        Perc_ratios_sum_Tau_AD(k,:) = Perc_ratios_sum_Tau(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Tau_WT(k,:) = Perc_ratios_Tau(k,:);
        Perc_ratios_sum_Tau_WT(k,:) = Perc_ratios_sum_Tau(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_Tau_AD(k,:) = Perc_ratios_Tau(k,:);
        Perc_ratios_sum_Tau_AD(k,:) = Perc_ratios_sum_Tau(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Tau_WT(k,:) = Perc_ratios_Tau(k,:);
        Perc_ratios_sum_Tau_WT(k,:) = Perc_ratios_sum_Tau(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Tau_WT(k,:) = Perc_ratios_Tau(k,:);
        Perc_ratios_sum_Tau_WT(k,:) = Perc_ratios_sum_Tau(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_Tau_WT(k,:) = Perc_ratios_Tau(k,:);
        Perc_ratios_sum_Tau_WT(k,:) = Perc_ratios_sum_Tau(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% NAA
%Defining variables
Perc_ratios_NAA_AD = nan(15,1);
Perc_ratios_NAA_WT = nan(15,1);

Perc_ratios_sum_NAA_AD = nan(15,1);
Perc_ratios_sum_NAA_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_NAA_AD(k,:) = Perc_ratios_NAA(k,:);
        Perc_ratios_sum_NAA_AD(k,:) = Perc_ratios_sum_NAA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_NAA_WT(k,:) = Perc_ratios_NAA(k,:);
        Perc_ratios_sum_NAA_WT(k,:) = Perc_ratios_sum_NAA(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_NAA_AD(k,:) = Perc_ratios_NAA(k,:);
        Perc_ratios_sum_NAA_AD(k,:) = Perc_ratios_sum_NAA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_NAA_WT(k,:) = Perc_ratios_NAA(k,:);
        Perc_ratios_sum_NAA_WT(k,:) = Perc_ratios_sum_NAA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_NAA_WT(k,:) = Perc_ratios_NAA(k,:);
        Perc_ratios_sum_NAA_WT(k,:) = Perc_ratios_sum_NAA(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_NAA_WT(k,:) = Perc_ratios_NAA(k,:);
        Perc_ratios_sum_NAA_WT(k,:) = Perc_ratios_sum_NAA(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% tNAA
%Defining variables
Perc_ratios_tNAA_AD = nan(15,1);
Perc_ratios_tNAA_WT = nan(15,1);

Perc_ratios_sum_tNAA_AD = nan(15,1);
Perc_ratios_sum_tNAA_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_tNAA_AD(k,:) = Perc_ratios_tNAA(k,:);
        Perc_ratios_sum_tNAA_AD(k,:) = Perc_ratios_sum_tNAA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tNAA_WT(k,:) = Perc_ratios_tNAA(k,:);
        Perc_ratios_sum_tNAA_WT(k,:) = Perc_ratios_sum_tNAA(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_tNAA_AD(k,:) = Perc_ratios_tNAA(k,:);
        Perc_ratios_sum_tNAA_AD(k,:) = Perc_ratios_sum_tNAA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tNAA_WT(k,:) = Perc_ratios_tNAA(k,:);
        Perc_ratios_sum_tNAA_WT(k,:) = Perc_ratios_sum_tNAA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tNAA_WT(k,:) = Perc_ratios_tNAA(k,:);
        Perc_ratios_sum_tNAA_WT(k,:) = Perc_ratios_sum_tNAA(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_tNAA_WT(k,:) = Perc_ratios_tNAA(k,:);
        Perc_ratios_sum_tNAA_WT(k,:) = Perc_ratios_sum_tNAA(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% tCr
%Defining variables
Perc_ratios_tCr_AD = nan(15,1);
Perc_ratios_tCr_WT = nan(15,1);

Perc_ratios_sum_tCr_AD = nan(15,1);
Perc_ratios_sum_tCr_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_tCr_AD(k,:) = Perc_ratios_tCr(k,:);
        Perc_ratios_sum_tCr_AD(k,:) = Perc_ratios_sum_tCr(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tCr_WT(k,:) = Perc_ratios_tCr(k,:);
        Perc_ratios_sum_tCr_WT(k,:) = Perc_ratios_sum_tCr(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_tCr_AD(k,:) = Perc_ratios_tCr(k,:);
        Perc_ratios_sum_tCr_AD(k,:) = Perc_ratios_sum_tCr(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tCr_WT(k,:) = Perc_ratios_tCr(k,:);
        Perc_ratios_sum_tCr_WT(k,:) = Perc_ratios_sum_tCr(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tCr_WT(k,:) = Perc_ratios_tCr(k,:);
        Perc_ratios_sum_tCr_WT(k,:) = Perc_ratios_sum_tCr(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_tCr_WT(k,:) = Perc_ratios_tCr(k,:);
        Perc_ratios_sum_tCr_WT(k,:) = Perc_ratios_sum_tCr(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% tCho
%Defining variables
Perc_ratios_tCho_AD = nan(15,1);
Perc_ratios_tCho_WT = nan(15,1);

Perc_ratios_sum_tCho_AD = nan(15,1);
Perc_ratios_sum_tCho_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_tCho_AD(k,:) = Perc_ratios_tCho(k,:);
        Perc_ratios_sum_tCho_AD(k,:) = Perc_ratios_sum_tCho(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tCho_WT(k,:) = Perc_ratios_tCho(k,:);
        Perc_ratios_sum_tCho_WT(k,:) = Perc_ratios_sum_tCho(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_tCho_AD(k,:) = Perc_ratios_tCho(k,:);
        Perc_ratios_sum_tCho_AD(k,:) = Perc_ratios_sum_tCho(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tCho_WT(k,:) = Perc_ratios_tCho(k,:);
        Perc_ratios_sum_tCho_WT(k,:) = Perc_ratios_sum_tCho(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_tCho_WT(k,:) = Perc_ratios_tCho(k,:);
        Perc_ratios_sum_tCho_WT(k,:) = Perc_ratios_sum_tCho(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_tCho_WT(k,:) = Perc_ratios_tCho(k,:);
        Perc_ratios_sum_tCho_WT(k,:) = Perc_ratios_sum_tCho(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% GABA
%Defining variables
Perc_ratios_GABA_AD = nan(15,1);
Perc_ratios_GABA_WT = nan(15,1);

Perc_ratios_sum_GABA_AD = nan(15,1);
Perc_ratios_sum_GABA_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_GABA_AD(k,:) = Perc_ratios_GABA(k,:);
        Perc_ratios_sum_GABA_AD(k,:) = Perc_ratios_sum_GABA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_GABA_WT(k,:) = Perc_ratios_GABA(k,:);
        Perc_ratios_sum_GABA_WT(k,:) = Perc_ratios_sum_GABA(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_GABA_AD(k,:) = Perc_ratios_GABA(k,:);
        Perc_ratios_sum_GABA_AD(k,:) = Perc_ratios_sum_GABA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_GABA_WT(k,:) = Perc_ratios_GABA(k,:);
        Perc_ratios_sum_GABA_WT(k,:) = Perc_ratios_sum_GABA(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_GABA_WT(k,:) = Perc_ratios_GABA(k,:);
        Perc_ratios_sum_GABA_WT(k,:) = Perc_ratios_sum_GABA(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_GABA_WT(k,:) = Perc_ratios_GABA(k,:);
        Perc_ratios_sum_GABA_WT(k,:) = Perc_ratios_sum_GABA(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% Gln
%Defining variables
Perc_ratios_Gln_AD = nan(15,1);
Perc_ratios_Gln_WT = nan(15,1);

Perc_ratios_sum_Gln_AD = nan(15,1);
Perc_ratios_sum_Gln_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_Gln_AD(k,:) = Perc_ratios_Gln(k,:);
        Perc_ratios_sum_Gln_AD(k,:) = Perc_ratios_sum_Gln(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Gln_WT(k,:) = Perc_ratios_Gln(k,:);
        Perc_ratios_sum_Gln_WT(k,:) = Perc_ratios_sum_Gln(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_Gln_AD(k,:) = Perc_ratios_Gln(k,:);
        Perc_ratios_sum_Gln_AD(k,:) = Perc_ratios_sum_Gln(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Gln_WT(k,:) = Perc_ratios_Gln(k,:);
        Perc_ratios_sum_Gln_WT(k,:) = Perc_ratios_sum_Gln(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_Gln_WT(k,:) = Perc_ratios_Gln(k,:);
        Perc_ratios_sum_Gln_WT(k,:) = Perc_ratios_sum_Gln(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_Gln_WT(k,:) = Perc_ratios_Gln(k,:);
        Perc_ratios_sum_Gln_WT(k,:) = Perc_ratios_sum_Gln(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% ML9
%Defining variables
Perc_ratios_ML9_AD = nan(15,1);
Perc_ratios_ML9_WT = nan(15,1);

Perc_ratios_sum_ML9_AD = nan(15,1);
Perc_ratios_sum_ML9_WT = nan(15,1);

for k=1:15
    if  strcmp(Type(1,k),'AD')
        Perc_ratios_ML9_AD(k,:) = Perc_ratios_ML9(k,:);
        Perc_ratios_sum_ML9_AD(k,:) = Perc_ratios_sum_ML9(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_ML9_WT(k,:) = Perc_ratios_ML9(k,:);
        Perc_ratios_sum_ML9_WT(k,:) = Perc_ratios_sum_ML9(k,:);
    elseif strcmp(Type(1,k),'AD')
        Perc_ratios_ML9_AD(k,:) = Perc_ratios_ML9(k,:);
        Perc_ratios_sum_ML9_AD(k,:) = Perc_ratios_sum_ML9(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_ML9_WT(k,:) = Perc_ratios_ML9(k,:);
        Perc_ratios_sum_ML9_WT(k,:) = Perc_ratios_sum_ML9(k,:);
    elseif strcmp(Type(1,k),'WT') 
        Perc_ratios_ML9_WT(k,:) = Perc_ratios_ML9(k,:);
        Perc_ratios_sum_ML9_WT(k,:) = Perc_ratios_sum_ML9(k,:);
    elseif strcmp(Type(1,k),'WT')
        Perc_ratios_ML9_WT(k,:) = Perc_ratios_ML9(k,:);
        Perc_ratios_sum_ML9_WT(k,:) = Perc_ratios_sum_ML9(k,:);
    else
        printf('Warning, the type/treatment factor have not been understood');
    end
end
%% Save percentages data
save('%_data.mat','Perc*');

