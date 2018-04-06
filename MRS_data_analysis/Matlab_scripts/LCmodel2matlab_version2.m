% LC model to matlab script, Maria Yanez Lopez, 27/02/14
% This version reads 2 spectra: right and left hand side (named
% 1-2). Selected metabolites are
% Glu,Ins,Tau,NAA,tNAA,tCr,tCho,Glx,GABA,Gln,ML9
% Input: Varian scan file, with LCmodel table file inside. Ex: HF_2014060901
% NOTES: This script requires distinguishable_colors.m
% OUTPUT: 
        % .mat with the name of the scan: 
%             --Whole_dataset {1,2} cell: including all metabolites
%                 -Metabolites: {35,1} cell (name of the metabolites)
%                 -Concentrations: {35,2} double (first column: abs values,
%                 -second column: ratios to tCr).
%                 -SD: {35,1} double (Cramer Rao)
%                 -FWHM: {1,1} double (linewidth estimation from LC model, in ppm)
%                 -SN: {1,1} double (signal to noise ratio)
%                 -Filename: {1,1} char (name of the experiment folder)
%             --Reliable dataset {1,2} cell: same but only including
%             Glu,Ins,Tau,NAA,tNAA,tCr,tCho,Glx,Gln, GABA, ML9
        % Graph with 11 metabolites in 2 points

clear all;
clc;
places = 2; %Right and left hand side

%Predefining cells for speed
SN = cell(1,places);
FWHM = cell(1,places);
metabolites = cell(1,places);
concentrations = cell(1,places);
Sdev = cell(1,places);


for scans=1:2 %Adjusments need to be made here when not a full dataset
    % The LC model results might be directly inside the HF_ folder or in the subfolder LC_model inside of it.      
    myfile = ['\Varian_LPS_experiment_JanFeb_2014\HF_2014013101\' strcat('0',num2str(scans),'.fid') '\LCmodel\table'];
    fileID = fopen(myfile);

    %Read introduction lines, to later extract the name of the scan:
    InputText = textscan(fileID,'%s',6,'delimiter','\n');  % Read strings delimited by a carriage return
    Intro = InputText{1};


    % Read header line
    InputText = textscan(fileID,'%s',1,'delimiter','\n');
    HeaderLines = InputText{1};


    %Read data
    InputText = textscan(fileID,'%f %s %f %s');

    %Read FWHM and S/N
    InputTextextra = textscan(fileID,'%s',1,'HeaderLines',1,'delimiter','\n');
    InputTextextra = InputTextextra{1};
    InputTextextra = char(InputTextextra);
    fclose(fileID);

    %Remove the % sign from the percentages column
    for i=1:length(InputText{1,2})
        temp1 = InputText{1,2}(i);
        temp1 = cell2mat(temp1);
        temp2 = temp1(1:end-1);
        temp2 = {temp2};
        InputText{1,2}(i) = temp2;
    end

    %Name and organise the data
    Concentrations = zeros(length(InputText{1,2}),2);
    Concentrations (:,1) = InputText{1,1};
    Concentrations (:,2) = InputText{1,3};

    Metabolites = InputText{1,4};
    SD_temp = InputText{1,2};
    for i=1:length(SD_temp)
        SD_temp{i,:} = str2double(SD_temp{i,:});
    end

    SD = cell2mat(SD_temp);
    Filename = Intro{2,1}(1:13);
    Reliable_Concentrations = Concentrations;
    Reliable_Metabolites = Metabolites;
    FWHM_temp = str2double(InputTextextra(8:12));
    SNratio = str2double(InputTextextra(28:29));
    
    %in cells
    SN{1,scans} = SNratio;
    FWHM{1,scans} = FWHM_temp;
    Fname{1,scans} = Intro{2,1};
    metabolites{1,scans} = Metabolites;
    concentrations{1,scans} = Concentrations;
    Sdev{1,scans} = SD;

    clear('fileID','temp1','temp2','InputText','InputTextextra','Intro','SD_temp', 'myfile','SNratio', 'FWHM_temp','HeaderLines');

    %Reject all metabolites but GABA, Gln, Glu,Ins,Tau,NAA,tNAA,tCr,tCho,Glx
    ind = logical([1 1 1 1 0 1 0 0 1 1 1 0 1 0 1 1 0 1 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 0 1]);
        for i=1:length(ind)
            if  ind(i)==1;
                Reliable_Concentrations(i,:) = 0;
                Reliable_Metabolites{i,1} = 'unreliable';
            end
        end
   
    %Name and organise the reduced dataset
    Reliable_Concentrations = removerows(Reliable_Concentrations,ind);
    Reliable_SD = removerows(SD,ind);
    ind = strmatch('unreliable',Reliable_Metabolites(:,1));
    Reliable_Metabolites = char(Reliable_Metabolites);
    Reliable_Metabolites(ind,:)=[];
    
    %in cells
    Reliable_metabolites{1,scans} = Reliable_Metabolites;
    Reliable_concentrations{1,scans} = Reliable_Concentrations;
    Reliable_Sdev{1,scans} = Reliable_SD;

    clear('ind','i','col','Metabolites','Concentrations','SD','Reliable_Metabolites','Reliable_Concentrations','Reliable_SD');

end

clear('scans','ans');

figure('units','normalized','outerposition',[0 0 1 1]);
Timepoints = char('PBS','LPS');
Metabolite_time_course = zeros(11,places);
colours = distinguishable_colors(11);


for h=1:1:11
    linecolour=colours(h,:);
    for k=1:1:length(Reliable_concentrations)
        Metabolite_time_course(h,k) = Reliable_concentrations{1,k}(h,2);     
    end 
    h1a = plot(1:1:places,Metabolite_time_course(h,:),'-mo');
    set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour, 'LineWidth',3);
    hold on;
end

title(sprintf('Mouse experiment %s metabolites PBS and LPS sites, ratio to tCr',Filename),'FontSize',16,'FontWeight','bold','interpreter','none')
set(gca, 'XTick',1:length(Reliable_concentrations), 'XTickLabel',Timepoints,'FontSize',16);
legend1=legend(gca,Reliable_metabolites{1,1});
set(legend1,'Location','EastOutside','FontSize',16,'Box','off');
ylim([0 max(max(Metabolite_time_course))]);

clear('Timepoints','colours','h','h1a','k','legend1','linecolour');

%Saving the data
    %Scans 
for i = 1:places    
   Whole_dataset{i}.Metabolites = metabolites{1,1};
   Whole_dataset{i}.Concentrations = concentrations{1,i};
   Whole_dataset{i}.SD = Sdev{1,i};
   Reliable_dataset{i}.Metabolites = Reliable_metabolites{1,1};
   Reliable_dataset{i}.Concentrations = Reliable_concentrations{1,i};
   Reliable_dataset{i}.SD = Reliable_Sdev{1,i};
   Whole_dataset{i}.FWHM = FWHM{1,i};
   Whole_dataset{i}.SN = SN{1,i};
   Whole_dataset{i}.Filename = Fname{1,i};
end

filename = strcat(Filename, '.mat');

save(filename, 'Whole_dataset', 'Reliable_dataset');

clear('places','i');