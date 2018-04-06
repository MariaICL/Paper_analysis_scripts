% Studying linewidths and reorganising data, Maria Yanez Lopez, 5/03/14
% This version reads.mat files from LCmodel2matlab_version2.m, renamed to
% the order of scanning (n=38 mice).
% Input: .mat files from LCmodel2matlab_version2.m. Ex: 1.mat
% NOTES: This script requires distinguishable_colors.m, aboxplot.m,
% colorgrad.m
% OUTPUT: 
        % FWHM (38,6) double (ppm)
%         FWHM_Hz (38,6) double (Hz)
%         SN (38,6) double (signal to noise ratio)
%         abs_metabolite (38,6) double (absolute values for 8 selected metabolites)
%         ratios_metabolite (38,6) double (ratios to tCr values for 8 selected metabolites) 
% Several graphs of FWHM and SN.


clear all;
clc;
%Load the whole dataset
for i = [1:38] %Discarding 6,15:19 mice
    iname = num2str(i);
    filename = strcat(iname,'.mat');
    Mouse(i) = load(filename);
    
end

clear('i','iname','filename');
%Predefine variables 
FWHM = zeros(38,6);
SN = zeros(38,6);

for i = [1:5 7:16 18:38]
    for j=1:6
        FWHM(i,j) = Mouse(1,i).Whole_dataset{1,j}.FWHM;
        SN(i,j) = Mouse(1,i).Whole_dataset{1,j}.SN;
    end
end

for i = [1:5 7:14 20:38]
    for j=1:6
        FWHM(i,j) = Mouse(1,i).Whole_dataset{1,j}.FWHM;
        SN(i,j) = Mouse(1,i).Whole_dataset{1,j}.SN;
    end
end

%Remove metabolites:
ind = logical([0 0 0 0 0 1 0 0 0 0 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]);
FWHM = removerows(FWHM,ind);
SN = removerows(SN,ind);

FWHM2 = reshape(FWHM,32*6,1);
SN2 = reshape(SN,32*6,1);


% figure(1)
% h = {FWHM};
% aboxplot(h,'labels',[0,0,1,2,3,4],'colorgrad','orange_down','OutlierMarker','h','OutlierMarkerSize',10);
% title(sprintf('FWHM'),'FontSize',16,'FontWeight','bold','interpreter','none')
% xlabel('hours','FontSize',16,'FontWeight','bold'); % Set the X-axis label
% ylabel('ppm','FontSize',16,'FontWeight','bold'); % Set the Y-axis label
% ylim([min(min(FWHM)) max(max(FWHM))]);
% % legend('%f Responders','%f Non responders'); % Add a legend
% legend('FWHM','Location','best'); 

% figure(2)
% h = {SN};
% aboxplot(h,'labels',[0,0,1,2,3,4],'colorgrad','green_down','OutlierMarker','h','OutlierMarkerSize',10);
% title(sprintf('SN'),'FontSize',16,'FontWeight','bold','interpreter','none')
% xlabel('hours','FontSize',16,'FontWeight','bold'); % Set the X-axis label
% ylabel('ratio','FontSize',16,'FontWeight','bold'); % Set the Y-axis label
% ylim([min(min(SN)) max(max(SN))]);
% % legend('%f Responders','%f Non responders'); % Add a legend
% legend('SN','Location','best'); 

FWHM_Hz = zeros(38,6);
FWHM_Hz = FWHM*300;

% figure(3)
% h = {FWHM_Hz;SN};
% aboxplot(h,'labels',[0,0,1,2,3,4],'colorgrad','orange_down','OutlierMarker','h','OutlierMarkerSize',10);
% title(sprintf('FWHM (Hz) and SNR'),'FontSize',16,'FontWeight','bold','interpreter','none')
% xlabel('timepoints','FontSize',16,'FontWeight','bold'); % Set the X-axis label
% ylabel('','FontSize',16,'FontWeight','bold'); % Set the Y-axis label
% ylim([min(min(SN)) max(max(SN))]);
% % legend('%f Responders','%f Non responders'); % Add a legend
% legend('FWHM','SN','Location','best'); 

figure(4)
h = {FWHM_Hz};
aboxplot(h,'labels',[0,0,1,2,3,4],'colorgrad','orange_down','OutlierMarker','h','OutlierMarkerSize',10);
title(sprintf('FWHM full dataset, 38 mice'),'FontSize',16,'FontWeight','bold','interpreter','none')
xlabel('hours','FontSize',16,'FontWeight','bold'); % Set the X-axis label
ylabel('Hz','FontSize',16,'FontWeight','bold'); % Set the Y-axis label
ylim([min(min(FWHM_Hz)) max(max(FWHM_Hz))]);
% legend('%f Responders','%f Non responders'); % Add a legend
legend('FWHM','Location','best'); 

