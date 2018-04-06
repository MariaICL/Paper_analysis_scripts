function [B1lowI, B1highI,FlipAngle,scandirectory]=loadVarianB1data(LFscan,HFscan,scandirectory)

%% ****  B1 data load Function for Varian **** %%

%INPUT: 
%       LFscan = Low Flip Angle scan number.
%
%       HFscan = High Flip Angle scan number 
%
%       scandirectory = path to the scan folder
%                       eg: 'C:\Users\msxmy\Desktop\CEST\Raw_data\HF_121101_2.gJ1'
%                       You are not obliged to give this input: if not, you
%                       will be asked to choose the scan directory manually.


%OUTPUT:
%       B1lowI = images obtained with the Low Flip Angle 
%                (d1 x d2 x d3 array) with: [d1,d2] the dimensions of each image
%                                            d3 the number of slices
%
%       B1highI = images obtained with the High Flip Angle 
%                (d1 x d2 x d3 array) with: [d1,d2] the dimensions of each image
%                                           d3 the number of slices
%
%      FlipAngle = vector containing the values (in �) of the Low Flip angle 
%                  and High Flip Angle (2-elements vector)
%
%      scandirectory = same as the input argument.



%% Selecting the scan main folder

%     if nargin<4
%         
%        scandirectory = uigetdir('','Select scan directory');
%        
%        if scandirectory == 0
%           error('User cancelled main folder selection')
%        end
%           
%     end    
    
%% Check the Low flip angle scan folder
 
    if LFscan < 10
       LFpath = fullfile(scandirectory,sprintf('0%d.img',LFscan));
    else 
       LFpath = fullfile(scandirectory,sprintf('%d.img',LFscan));
    end
    
    % Retrieve the names of the images in the LFA scan folder
    dirlist = dir(LFpath);
    LFfiles = {};
    m = 1;
    
    for i = 3:length(dirlist) 
        
         [~,~,ext] = fileparts(dirlist(i).name);
         
         if strcmp(ext,'.fdf') || strcmp(ext,'.nii')|| strcmp(ext,'.gz')|| strcmp(ext,'.hdr')
            LFfiles{m} = dirlist(i).name;
            m = m+1;
         end
         
    end

    if isempty(LFfiles) 
       error('The LFA directory is empty')
    end

    % Retrieve the number of slices/images/echoes
    LFext = cell(1,length(LFfiles));
    imageLF = zeros(1,length(LFfiles));
    echoLF = zeros(1,length(LFfiles));

    for i=1:length(LFfiles)
        
        [~,name,LFext{i}] = fileparts(LFfiles{i});
        idx1 = find(name == 'i',1,'last');
        idx2 = find(name=='e',1,'last');
        imageLF(i) = str2num(name(idx1+5:idx2-1));
        echoLF(i) = str2num(name(idx2+4:end));
        
    end

    NimLF = max(imageLF);
    NechoLF = max(echoLF);

    % Check the number of images and echos
    if NimLF ~= 1
        error('The LFA folder must only contain one image per slice')
    end
    
    if NechoLF ~= 1
        error('LoadVarianB1Data does not work for several echo images')
    end

                    
%% Check the High flip angle scan folder
 
    if HFscan < 10
       HFpath = fullfile(scandirectory,sprintf('0%d.img',HFscan));
    else 
       HFpath = fullfile(scandirectory,sprintf('%d.img',HFscan));
    end
        
    % Retrieve the names of the images in the LFA scan folder
    dirlist = dir(HFpath);
    HFfiles = {};
    m = 1;
    
    for i=3:length(dirlist)   
        
         [~,~,ext] = fileparts(dirlist(i).name);
         
         if strcmp(ext,'.fdf') || strcmp(ext,'.nii')|| strcmp(ext,'.gz')|| strcmp(ext,'.hdr')
            HFfiles{m} = dirlist(i).name;
            m = m+1;
         end
         
    end

    if isempty(HFfiles) 
       error('The HFA directory is empty')
    end

    % Retrieve the number of slices/images/echoes
    HFext = cell(1,length(HFfiles));
    imageHF = zeros(1,length(HFfiles));
    echoHF = zeros(1,length(HFfiles));

    for i = 1:length(HFfiles)
        
        [~,name,HFext{i}] = fileparts(HFfiles{i});
        idx1 = find(name == 'i',1,'last');
        idx2 = find(name=='e',1,'last');
        imageHF(i) = str2num(name(idx1+5:idx2-1));
        echoHF(i) = str2num(name(idx2+4:end));
        
    end

    NimHF = max(imageHF);
    NechoHF = max(echoHF);

    % Check the number of images and echos
    if NimHF ~= 1
        error('The HFA folder must only contain one image per slice')
    end
    
    if NechoHF ~= 1
        error('LoadVarianB1Data does not work for several echo images')
    end

    
%% Read procpar files to determine flip angles
    
    B1procpar.file{1} = fullfile(LFpath,'procpar');
    B1procpar.file{2} = fullfile(HFpath,'procpar');
    
    FlipAngle=zeros(1,2);
    
    for b=1:2
        
        if b == 1
            fprintf('     Reading procpar file for low FA data')
        elseif b == 2
            fprintf('     Reading procpar file for high FA data')
        end
        
        B1procpar.fileID(b) = fopen(B1procpar.file{b});
        
        if B1procpar.fileID(b) ~= -1
            
            while(feof(B1procpar.fileID(b)) == 0)
                 line = fgetl(B1procpar.fileID(b));
                
                 % Excitation flip angle
                 if strncmp(line,'flip1 ',6) == 1
                    line = fgetl(B1procpar.fileID(b));
                    flip1_tmp = str2num(line);
                    FlipAngle(b) = flip1_tmp(1,2);
                end    
                                
            end
            
        elseif b == 1
            error(' Procpar file not found in LFA folder')
        elseif b == 2
            error('Procpar file not found in HFA folder')
        end
        
    end
    
    % Checking errors
    if FlipAngle(1) ~= FlipAngle(2)/2
        error('The Low Flip Angle ( %d ) must be equal to half of the High Flip Angle ( %d )',...
               FlipAngle(1),FlipAngle(2))
    end

    
%% Load Low FA scans

    disp('     Loading Low Flip Angle images')
    
    for i = 1:length(LFfiles)
        
        if strcmpi(LFext{i},'.nii') || strcmpi(LFext{i},'.gz') || strcmpi(LFext{i},'.hdr')
            
           LFfile = fullfile(LFpath,LFfiles{i});
           LF_tmp  = aedes_read_nifti(LFfile);
           B1lowI(:,:) = LF_tmp.FTDATA;
            
        elseif strcmpi(LFext{i},'.fdf')
            
            LFfile = fullfile(LFpath,LFfiles{i});
            LF_tmp  = aedes_readfdf(LFfile);
            B1lowI(:,:) = LF_tmp.FTDATA;
            
        end
        
    end
    
%% Load High FA scans

    disp('     Loading High Flip Angle images')
    
    for i = 1:length(HFfiles)
        
        if strcmpi(HFext{i},'.nii') || strcmpi(HFext{i},'.gz') || strcmpi(HFext{i},'.hdr')

            HFfile = fullfile(HFpath,HFfiles{i});
            HF_tmp  = aedes_read_nifti(HFfile);
            B1highI(:,:) = HF_tmp.FTDATA;

        elseif strcmpi(HFext{i},'.fdf')

            HFfile = fullfile(HFpath,HFfiles{i});
            HF_tmp  = aedes_readfdf(HFfile);
            B1highI(:,:) = HF_tmp.FTDATA;

        end
        
    end

%% Checking the dimensions

    if size(B1lowI,1) ~= size(B1highI,1) || size(B1lowI,2) ~= size(B1highI,2)
       error('The LFA array and HFA array do not have the same dimensions')
    end