clear all;
clc;
%Load data
load('SN.mat');
load('FWHM_Hz.mat');
load('Type.mat');
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
load('%_data');
ID = [1:15]';

%Create excell sheet
filename = 'MRS_full_dataset.xlsx';

% Write Excluded mice, Type, Subtype, Treatment to first sheet
xlswrite(filename,Type',1,'A2');
xlswrite(filename,FWHM_Hz,1,'B2');
xlswrite(filename,SN,1,'D2');
% xlswrite(filename,Excluded_mice',1,'H2');
xlswrite(filename,{'Type','FWHM_PBS','FWHM_LPS','SN_PBS','SN_LPS','ID','Mouse ID','Wax ID','Date','Mice excluded','Selected metabolites'},1,'A1');
xlswrite(filename,{'Glu','Ins','Tau','NAA','tNAA','tCr','tCho','Glx','Gln','GABA','ML9'}',1,'K2');
xlswrite(filename,ID,1,'F2');

%%
%Writing bulk data for 11 selected metabolites in 11 separate sheets
for i=2:12
    xlswrite(filename,ID,i,'A3');    
    xlswrite(filename,{'ID','PBS','LPS'},i,'A2');
    xlswrite(filename,{'PBS','LPS'},i,'L2');
    xlswrite(filename,{'PBS','LPS'},i,'E2');    
    xlswrite(filename,{'% LPS-PBS'},i,'H2');
    xlswrite(filename,{'% LPS-PBS'},i,'J2'); 
    xlswrite(filename,{'% LPS-PBS'},i,'O2');
end

%Ins
xlswrite(filename,ratios_Ins,2,'B3');
xlswrite(filename,{'Ins ratios to tCr'},2,'A1');

xlswrite(filename,ratios_sum_Ins,2,'E3');
xlswrite(filename,{'Ins ratios to sum selected metabolites'},2,'E1');

xlswrite(filename,abs_Ins,2,'L3');
xlswrite(filename,{'Ins absolute values'},2,'L1');

%Glu
xlswrite(filename,ratios_Glu,3,'B3');
xlswrite(filename,{'Glu ratios to tCr'},3,'A1');

xlswrite(filename,ratios_sum_Glu,3,'E3');
xlswrite(filename,{'Glu ratios to sum selected metabolites'},3,'E1');

xlswrite(filename,abs_Glu,3,'L3');
xlswrite(filename,{'Glu absolute values'},3,'L1');

%Tau
xlswrite(filename,ratios_Tau,4,'B3');
xlswrite(filename,{'Tau ratios to tCr'},4,'A1');

xlswrite(filename,ratios_sum_Tau,4,'E3');
xlswrite(filename,{'Tau ratios to sum selected metabolites'},4,'E1');

xlswrite(filename,abs_Tau,4,'L3');
xlswrite(filename,{'Tau absolute values'},4,'L1');

%NAA
xlswrite(filename,ratios_NAA,5,'B3');
xlswrite(filename,{'NAA ratios to tCr'},5,'A1');

xlswrite(filename,ratios_sum_NAA,5,'E3');
xlswrite(filename,{'NAA ratios to sum selected metabolites'},5,'E1');

xlswrite(filename,abs_NAA,5,'L3');
xlswrite(filename,{'NAA absolute values'},5,'L1');

%tNAA
xlswrite(filename,ratios_tNAA,6,'B3');
xlswrite(filename,{'tNAA ratios to tCr'},6,'A1');

xlswrite(filename,ratios_sum_tNAA,6,'E3');
xlswrite(filename,{'tNAA ratios to sum selected metabolites'},6,'E1');

xlswrite(filename,abs_tNAA,6,'L3');
xlswrite(filename,{'tNAA absolute values'},6,'L1');

%tCho
xlswrite(filename,ratios_tCho,7,'B3');
xlswrite(filename,{'tCho ratios to tCr'},7,'A1');

xlswrite(filename,ratios_sum_tCho,7,'E3');
xlswrite(filename,{'tCho ratios to sum selected metabolites'},7,'E1');

xlswrite(filename,abs_tCho,7,'L3');
xlswrite(filename,{'tCho absolute values'},7,'L1');

%tCr
xlswrite(filename,ratios_tCr,8,'B3');
xlswrite(filename,{'tCr ratios to tCr'},8,'A1');

xlswrite(filename,ratios_sum_tCr,8,'E3');
xlswrite(filename,{'tCr ratios to sum selected metabolites'},8,'E1');

xlswrite(filename,abs_tCr,8,'L3');
xlswrite(filename,{'tCr absolute values'},8,'L1');

%Glx
xlswrite(filename,ratios_Glx,9,'B3');
xlswrite(filename,{'Glx ratios to tCr'},9,'A1');

xlswrite(filename,ratios_sum_Glx,9,'E3');
xlswrite(filename,{'Glx ratios to sum selected metabolites'},9,'E1');

xlswrite(filename,abs_Glx,9,'L3');
xlswrite(filename,{'Glx absolute values'},9,'L1');

%Gln
xlswrite(filename,ratios_Gln,10,'B3');
xlswrite(filename,{'Gln ratios to tCr'},10,'A1');

xlswrite(filename,ratios_sum_Gln,10,'E3');
xlswrite(filename,{'Gln ratios to sum selected metabolites'},10,'E1');

xlswrite(filename,abs_Gln,10,'L3');
xlswrite(filename,{'Gln absolute values'},10,'L1');

%GABA
xlswrite(filename,ratios_GABA,11,'B3');
xlswrite(filename,{'GABA ratios to tCr'},11,'A1');

xlswrite(filename,ratios_sum_GABA,11,'E3');
xlswrite(filename,{'GABA ratios to sum selected metabolites'},11,'E1');

xlswrite(filename,abs_GABA,11,'L3');
xlswrite(filename,{'GABA absolute values'},11,'L1');

%ML9
xlswrite(filename,ratios_ML9,12,'B3');
xlswrite(filename,{'ML9 ratios to tCr'},12,'A1');

xlswrite(filename,ratios_sum_ML9,12,'E3');
xlswrite(filename,{'ML9 ratios to sum selected metabolites'},12,'E1');

xlswrite(filename,abs_ML9,12,'L3');
xlswrite(filename,{'ML9 absolute values'},12,'L1');

%% Writing % bulk data for 11 selected metabolites in 11 separate sheets
%Ins
xlswrite(filename,Perc_ratios_Ins,2,'H3');
xlswrite(filename,{'% change from baseline Ins ratios to tCr'},2,'H1');

xlswrite(filename,Perc_ratios_sum_Ins,2,'J3');
xlswrite(filename,{'% change from baseline Ins ratios to sum selected metabolites'},2,'J1');

xlswrite(filename,Perc_abs_Ins,2,'O3');
xlswrite(filename,{'% change from baseline Ins absolute values'},2,'O1');

%Glu
xlswrite(filename,Perc_ratios_Glu,3,'H3');
xlswrite(filename,{'% change from baseline Glu ratios to tCr'},3,'H1');

xlswrite(filename,Perc_ratios_sum_Glu,3,'J3');
xlswrite(filename,{'% change from baseline Glu ratios to sum selected metabolites'},3,'J1');

xlswrite(filename,Perc_abs_Glu,3,'O3');
xlswrite(filename,{'% change from baseline Glu absolute values'},3,'O1');

%Tau
xlswrite(filename,Perc_ratios_Tau,4,'H3');
xlswrite(filename,{'% change from baseline Tau ratios to tCr'},4,'H1');

xlswrite(filename,Perc_ratios_sum_Tau,4,'J3');
xlswrite(filename,{'% change from baseline Tau ratios to sum selected metabolites'},4,'J1');

xlswrite(filename,Perc_abs_Tau,4,'O3');
xlswrite(filename,{'% change from baseline Tau absolute values'},4,'O1');

%NAA
xlswrite(filename,Perc_ratios_NAA,5,'H3');
xlswrite(filename,{'% change from baseline NAA ratios to tCr'},5,'H1');

xlswrite(filename,Perc_ratios_sum_NAA,5,'J3');
xlswrite(filename,{'% change from baseline NAA ratios to sum selected metabolites'},5,'J1');

xlswrite(filename,Perc_abs_NAA,5,'O3');
xlswrite(filename,{'% change from baseline NAA absolute values'},5,'O1');

%tNAA
xlswrite(filename,Perc_ratios_tNAA,6,'H3');
xlswrite(filename,{'% change from baseline tNAA ratios to tCr'},6,'H1');

xlswrite(filename,Perc_ratios_sum_tNAA,6,'J3');
xlswrite(filename,{'% change from baseline tNAA ratios to sum selected metabolites'},6,'J1');

xlswrite(filename,Perc_abs_tNAA,6,'O3');
xlswrite(filename,{'% change from baseline tNAA absolute values'},6,'O1');

%tCho
xlswrite(filename,Perc_ratios_tCho,7,'H3');
xlswrite(filename,{'% change from baseline tCho ratios to tCr'},7,'H1');

xlswrite(filename,Perc_ratios_sum_tCho,7,'J3');
xlswrite(filename,{'% change from baseline tCho ratios to sum selected metabolites'},7,'J1');

xlswrite(filename,Perc_abs_tCho,7,'O3');
xlswrite(filename,{'% change from baseline tCho absolute values'},7,'O1');

%tCr
xlswrite(filename,Perc_ratios_tCr,8,'H3');
xlswrite(filename,{'% change from baseline tCr ratios to tCr'},8,'H1');

xlswrite(filename,Perc_ratios_sum_tCr,8,'J3');
xlswrite(filename,{'% change from baseline tCr ratios to sum selected metabolites'},8,'J1');

xlswrite(filename,Perc_abs_tCr,8,'O3');
xlswrite(filename,{'% change from baseline tCr absolute values'},8,'O1');

%Glx
xlswrite(filename,Perc_ratios_Glx,9,'H3');
xlswrite(filename,{'% change from baseline Glx ratios to tCr'},9,'H1');

xlswrite(filename,Perc_ratios_sum_Glx,9,'J3');
xlswrite(filename,{'% change from baseline Glx ratios to sum selected metabolites'},9,'J1');

xlswrite(filename,Perc_abs_Glx,9,'O3');
xlswrite(filename,{'% change from baseline Glx absolute values'},9,'O1');

%Gln
xlswrite(filename,Perc_ratios_Gln,10,'H3');
xlswrite(filename,{'% change from baseline Gln ratios to tCr'},10,'H1');

xlswrite(filename,Perc_ratios_sum_Gln,10,'J3');
xlswrite(filename,{'% change from baseline Gln ratios to sum selected metabolites'},10,'J1');

xlswrite(filename,Perc_abs_Gln,10,'O3');
xlswrite(filename,{'% change from baseline Gln absolute values'},10,'O1');

%GABA
xlswrite(filename,Perc_ratios_GABA,11,'H3');
xlswrite(filename,{'% change from baseline GABA ratios to tCr'},11,'H1');

xlswrite(filename,Perc_ratios_sum_GABA,11,'J3');
xlswrite(filename,{'% change from baseline GABA ratios to sum selected metabolites'},11,'J1');

xlswrite(filename,Perc_abs_GABA,11,'O3');
xlswrite(filename,{'% change from baseline GABA absolute values'},11,'O1');

%ML9
xlswrite(filename,Perc_ratios_ML9,12,'H3');
xlswrite(filename,{'% change from baseline ML9 ratios to tCr'},12,'H1');

xlswrite(filename,Perc_ratios_sum_ML9,12,'J3');
xlswrite(filename,{'% change from baseline ML9 ratios to sum selected metabolites'},12,'J1');

xlswrite(filename,Perc_abs_ML9,12,'O3');
xlswrite(filename,{'% change from baseline ML9 absolute values'},12,'O1');