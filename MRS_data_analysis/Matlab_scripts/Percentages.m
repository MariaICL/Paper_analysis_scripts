%%
% Percentages script, Maria Yanez Lopez, 
% This version reads all the files coming from
% Organise_data, such as Ins_data, Tau_data etc
% Input: metabolite_data files (8), from Organise_data
% It calculates the % increase ratio of LPS injection voxel compared with PBS.
%OUTPUT: 
        % %_data.mat containing the metabolite data

clear all;
clc;

%Load data
load('Ins_data');
load('Glu_data');
load('Glx_data');
load('Tau_data');
load('NAA_data');
load('tNAA_data');
load('tCho_data');
load('tCr_data');
load('Gln_data');
load('GABA_data');
load('ML9_data');

%Create variables        
Perc_ratios_Glu = zeros(15,1);
Perc_ratios_Ins = zeros(15,1);
Perc_ratios_NAA = zeros(15,1);
Perc_ratios_Tau = zeros(15,1);
Perc_ratios_tCho = zeros(15,1);
Perc_ratios_tNAA = zeros(15,1);
Perc_ratios_tCr = zeros(15,1);
Perc_ratios_Glx = zeros(15,1);
Perc_ratios_GABA = zeros(15,1);
Perc_ratios_Gln = zeros(15,1);
Perc_ratios_ML9 = zeros(15,1);

Perc_ratios_sum_Glu = zeros(15,1);
Perc_ratios_sum_Ins = zeros(15,1);
Perc_ratios_sum_NAA = zeros(15,1);
Perc_ratios_sum_Tau = zeros(15,1);
Perc_ratios_sum_tCho = zeros(15,1);
Perc_ratios_sum_tNAA = zeros(15,1);
Perc_ratios_sum_tCr = zeros(15,1);
Perc_ratios_sum_Glx = zeros(15,1);
Perc_ratios_sum_GABA = zeros(15,1);
Perc_ratios_sum_Gln = zeros(15,1);
Perc_ratios_sum_ML9 = zeros(15,1);

Perc_abs_Glu = zeros(15,1);
Perc_abs_Ins = zeros(15,1);
Perc_abs_NAA = zeros(15,1);
Perc_abs_Tau = zeros(15,1);
Perc_abs_tCho = zeros(15,1);
Perc_abs_tNAA = zeros(15,1);
Perc_abs_tCr = zeros(15,1);
Perc_abs_Glx = zeros(15,1);
Perc_abs_GABA = zeros(15,1);
Perc_abs_Gln = zeros(15,1);
Perc_abs_ML9 = zeros(15,1);

%Calculate percentage change
for i = 1:15
    
        Perc_ratios_Glu(i) = 100.*((ratios_Glu(i,2) - ratios_Glu(i,1))./ratios_Glu(i,1));
        Perc_ratios_Ins(i) = 100.*((ratios_Ins(i,2) - ratios_Ins(i,1))./ratios_Ins(i,1));
        Perc_ratios_NAA(i) = 100.*((ratios_NAA(i,2) - ratios_NAA(i,1))./ratios_NAA(i,1));
        Perc_ratios_Tau(i) = 100.*((ratios_Tau(i,2) - ratios_Tau(i,1))./ratios_Tau(i,1));
        Perc_ratios_tCho(i) = 100.*((ratios_tCho(i,2) - ratios_tCho(i,1))./ratios_tCho(i,1));
        Perc_ratios_tNAA(i) = 100.*((ratios_tNAA(i,2) - ratios_tNAA(i,1))./ratios_tNAA(i,1));
        Perc_ratios_tCr(i) = 100.*((ratios_tCr(i,2) - ratios_tCr(i,1))./ratios_tCr(i,1));
        Perc_ratios_Glx(i) = 100.*((ratios_Glx(i,2) - ratios_Glx(i,1))./ratios_Glx(i,1));
        Perc_ratios_GABA(i) = 100.*((ratios_GABA(i,2) - ratios_GABA(i,1))./ratios_GABA(i,1));
        Perc_ratios_ML9(i) = 100.*((ratios_ML9(i,2) - ratios_ML9(i,1))./ratios_ML9(i,1));
        Perc_ratios_Gln(i) = 100.*((ratios_Gln(i,2) - ratios_Gln(i,1))./ratios_Gln(i,1));
        
        Perc_ratios_sum_Glu(i) = 100.*((ratios_sum_Glu(i,2) - ratios_sum_Glu(i,1))./ratios_sum_Glu(i,1));
        Perc_ratios_sum_Ins(i) = 100.*((ratios_sum_Ins(i,2) - ratios_sum_Ins(i,1))./ratios_sum_Ins(i,1));
        Perc_ratios_sum_NAA(i) = 100.*((ratios_sum_NAA(i,2) - ratios_sum_NAA(i,1))./ratios_sum_NAA(i,1));
        Perc_ratios_sum_Tau(i) = 100.*((ratios_sum_Tau(i,2) - ratios_sum_Tau(i,1))./ratios_sum_Tau(i,1));
        Perc_ratios_sum_tCho(i) = 100.*((ratios_sum_tCho(i,2) - ratios_sum_tCho(i,1))./ratios_sum_tCho(i,1));
        Perc_ratios_sum_tNAA(i) = 100.*((ratios_sum_tNAA(i,2) - ratios_sum_tNAA(i,1))./ratios_sum_tNAA(i,1));
        Perc_ratios_sum_tCr(i) = 100.*((ratios_sum_tCr(i,2) - ratios_sum_tCr(i,1))./ratios_sum_tCr(i,1));
        Perc_ratios_sum_Glx(i) = 100.*((ratios_sum_Glx(i,2) - ratios_sum_Glx(i,1))./ratios_sum_Glx(i,1));
        Perc_ratios_sum_GABA(i) = 100.*((ratios_sum_GABA(i,2) - ratios_sum_GABA(i,1))./ratios_sum_GABA(i,1));
        Perc_ratios_sum_ML9(i) = 100.*((ratios_sum_ML9(i,2) - ratios_sum_ML9(i,1))./ratios_sum_ML9(i,1));
        Perc_ratios_sum_Gln(i) = 100.*((ratios_sum_Gln(i,2) - ratios_sum_Gln(i,1))./ratios_sum_Gln(i,1));
        
        Perc_abs_Glu(i) = 100.*((abs_Glu(i,2) - abs_Glu(i,1))./abs_Glu(i,1));
        Perc_abs_Ins(i) = 100.*((abs_Ins(i,2) - abs_Ins(i,1))./abs_Ins(i,1));
        Perc_abs_NAA(i) = 100.*((abs_NAA(i,2) - abs_NAA(i,1))./abs_NAA(i,1));
        Perc_abs_Tau(i) = 100.*((abs_Tau(i,2) - abs_Tau(i,1))./abs_Tau(i,1));
        Perc_abs_tCho(i) = 100.*((abs_tCho(i,2) - abs_tCho(i,1))./abs_tCho(i,1));
        Perc_abs_tNAA(i) = 100.*((abs_tNAA(i,2) - abs_tNAA(i,1))./abs_tNAA(i,1));
        Perc_abs_tCr(i) = 100.*((abs_tCr(i,2) - abs_tCr(i,1))./abs_tCr(i,1));
        Perc_abs_Glx(i) = 100.*((abs_Glx(i,2) - abs_Glx(i,1))./abs_Glx(i,1));
        Perc_abs_GABA(i) = 100.*((abs_GABA(i,2) - abs_GABA(i,1))./abs_GABA(i,1));
        Perc_abs_ML9(i) = 100.*((abs_ML9(i,2) - abs_ML9(i,1))./abs_ML9(i,1));
        Perc_abs_Gln(i) = 100.*((abs_Gln(i,2) - abs_Gln(i,1))./abs_Gln(i,1));
end

%% Save percentages data
save('%_data.mat','Perc*');