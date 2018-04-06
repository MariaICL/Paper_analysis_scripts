
function [I_ref, I_norm, ppm_ord, scandirectory] = loadVarianData(refscan,...
                                                                scan,type, scandirectory)

%% ****  CEST and WASSR data load Function for Varian  **** %%

%INPUT: 
%       refscans = reference scan number
%
%       scans = CEST or WASSR scan number 
%
%       type = type of data: 'CEST' for CEST data
%                            'WASSR' for WASSR data
%
%       scandirectory = path to the scan main folder
%                       eg: 'C:\Users\msxmy\Desktop\CEST\Raw_data\HF_121101_2.gJ1'
%                       You are not obliged to give this input: if not, you
%                       will be asked to choose the scan directory manually


%OUTPUT:
%       I_ref = reference images 
%              (d1 x d2 x d3 array) with: [d1,d2] the dimensions of each image
%                                          d3 the number of slices
%
%       I_norm = CEST or WASSR normalized images for all slices and all offset frequencies.
%               (d1 x d2 x d3 x nfreq array) with [d1,d2] the dimensions of each image
%                                                  d3 the number of slices
%                                                  nfreq the number of offset frequencies
%
%       ppm_ord = offset frequencies (in ppm)by increasing order 
%                (1xnfreq vector)
%
%       scandirectory = same as the input argument.

    
%% Selecting the scan main folder

    if nargin<4
        
       scandirectory = uigetdir('/media/OS/CEST_Varian','Select scan main folder');
       
       if scandirectory == 0
          error('User cancelled main folder selection')
       end
          
    end    


%% Check the reference scan folder    
    
    if strcmp(type,'CEST')
        
        if refscan < 10
           refpath = fullfile(scandirectory,sprintf('0%d.img',refscan));
        else 
           refpath = fullfile(scandirectory,sprintf('%d.img',refscan));
        end
               
        % Retrieve the names of the images in the reference scan folder
        dirlist = dir(refpath);
        reffiles = {};
        m = 1;
        
        for i = 3:length(dirlist)  
            
             [unused,unused,ext] = fileparts(dirlist(i).name);
             
             if strcmp(ext,'.fdf') || strcmp(ext,'.nii')|| strcmp(ext,'.gz')|| strcmp(ext,'.hdr')
                reffiles{m} = dirlist(i).name;
                m = m+1;
             end
             
        end
    
        if isempty(reffiles) 
            error('The reference scan folder is empty')
        end
        
        % Retrieve the number of slices/images/echoes
        refext = cell(1,length(reffiles));
        imageR = zeros(1,length(reffiles));
        echoR = zeros(1,length(reffiles));
        
        for i=1:length(reffiles)
            
            [unused,name,refext{i}] = fileparts(reffiles{i});
            idx1 = find(name == 'i',1,'last');
            idx2 = find(name == 'e',1,'last');
            imageR(i) = str2num(name(idx1+5:idx2-1));
            echoR(i) = str2num(name(idx2+4:end));
            
        end
            

        NimR = max(imageR);
        NechoR = max(echoR);
        
        % Check the number of images and echos
        if NimR ~= 1
            error('The reference scan folder must only contain one image per slice')
        end
        
        if NechoR ~= 1
            error('LoadVarianData does not work for several echo images')
        end
              
    elseif strcmp(type,'WASSR')
           
        
    else
        error('Type not understood')
    end
    
    
%% Check the CEST or WASSR folder

    if scan<10
       scanpath = fullfile(scandirectory,sprintf('0%d.img',scan));
    else 
       scanpath = fullfile(scandirectory,sprintf('%d.img',scan));
    end
    
    % Retrieve the names of the images in the WASSR or CEST folder
    dirlist = dir(scanpath);
    scanfiles = {};
    m = 1;
     
    for i = 3:length(dirlist)   
         
        [unused,unused,ext] = fileparts(dirlist(i).name);
         
         if strcmp(ext,'.fdf') || strcmp(ext,'.nii')|| strcmp(ext,'.gz')|| strcmp(ext,'.hdr')
            scanfiles{m} = dirlist(i).name;
            m = m+1;
         end
         
         
    end
     
    if isempty(scanfiles)
       if strcmp(type,'CEST')
          error('The CEST  folder is empty')
       elseif strcmp(type,'WASSR')
          error('The WASSR folder is empty')
       end
    end
    
    % Retrieve the number of slices/images/echoes
    scanext = cell(1,length(scanfiles));
    imageS = zeros(1,length(scanfiles));
    echoS = zeros(1,length(scanfiles));
      
    for i=1:length(scanfiles)
    
        [unused,name,scanext{i}]= fileparts(scanfiles{i});
        idx1=find(name=='i',1,'last');
        idx2=find(name=='e',1,'last');
        imageS(i)=str2num(name(idx1+5:idx2-1));
        echoS(i)=str2num(name(idx2+4:end));
        
     end
            

     NimS=max(imageS);
     NechoS=max(echoS);
       
     % Check the number of echos/slices
     if NechoS~=1
        error('LoadVarianData does not work for several echo images')
     end
       
     
%% Read procpar files to find the offset frequencies
    
    if strcmp(type,'CEST')
        disp('     Reading procpar files for CEST data')
    elseif strcmp(type, 'WASSR')
         disp('     Reading procpar files for WASSR data')
    else
        error('Type not understood')
    end
    
    H1fq = 399.453222707;
    disp('H1fq = 399.453222707, To avoid cine bug in Varian');
    
    procpar = fopen(fullfile(scanpath,'procpar'),'r','l');
    
    if procpar ~= -1
        
        while (feof(procpar) == 0)
            
            line = fgetl(procpar);
            
            if (isempty(strfind(line,'mtfrq')) == 0)
                line = fgetl(procpar);
                of_tmp = str2num(line);
                offsetfqs_tmp = of_tmp;
            end
            
            clear of_tmp
            
%             if (isempty(strfind(line,'H1reffrq')) == 0)
%                 line = fgetl(procpar);
%                 H1fq_tmp = str2num(line);
%                 H1fq = H1fq_tmp(1,2);
%             end
            
        end
        
    else
        
        if strcmp(type,'CEST')
           error('Procpar file not found in CEST folder ')
        elseif strcmp(type,'WASSR')
            error('Procpar file not found in WASSR folder ')
        end
 
    end
    
    offsetfqs = offsetfqs_tmp(:,2:end);
    
    if length(offsetfqs) ~= NimS
        if strcmp(type,'CEST')
           error('The number of offset frequencies is different from the number of CEST images')
        else
           error('The number of offset frequencies is different from the number of WASSR images')
        end
    end  

    fclose('all');

   
%% Load scans
    
    if strcmp(type, 'CEST') 
       progressbar('Loading CEST images');
       bar = 1;
    elseif strcmp(type, 'WASSR') 
       progressbar('Loading WASSR images');
       bar = 1;
    end
    

    %% Load reference scans
    
    if strcmp(type,'CEST')
    
        disp('     Loading CEST reference images')
        
        for i = 1:length(reffiles)
            
            %Determine the reading function
            if strcmpi(refext{i},'.nii') || strcmpi(refext{i},'.gz') || strcmpi(refext{i},'.hdr')                
                reffile = fullfile(refpath,reffiles{i});
                I_ref_tmp = aedes_read_nifti(reffile);
                I_ref(:,:) = I_ref_tmp.FTDATA;
                
            elseif strcmpi(refext{i},'.fdf')                
                reffile = fullfile(refpath,reffiles{i});
                I_ref_tmp = aedes_readfdf(reffile);
                I_ref(:,:) = I_ref_tmp.FTDATA;
                
            end
            
        end

    elseif strcmp(type,'WASSR')
        
    else
        error('Type not understood')
    end
    
    %% Load arrays
       
    if strcmp(type,'CEST')
        disp('     Loading CEST images')
    elseif strcmp(type, 'WASSR')
         disp('     Loading WASSR images')
    else
        error('Type not understood')
    end
 
    
     for i = 1:length(scanfiles)
            
         %Determine the reading function
          if strcmpi(scanext{i},'.nii') || strcmpi(scanext{i},'.gz') || strcmpi(scanext{i},'.hdr') 
             scanfile = fullfile(scanpath,scanfiles{i});
             I_tmp = aedes_read_nifti(scanfile);
             I(:,:,imageS(i)) = I_tmp.FTDATA;
                
          elseif strcmpi(scanext{i},'.fdf')             
              scanfile = fullfile(scanpath,scanfiles{i});
              I_tmp = aedes_readfdf(scanfile);
              I(:,:,imageS(i)) = I_tmp.FTDATA;
                
          end
            
          bar = bar+1;
          progressbar(bar/(length(scanext)));
            
     end    
           
     progressbar(1);
     

%%  Separate WASSR reference image from WASSR array

    if strcmp(type,'WASSR')
        
       % Assumes WASSR reference scan is a part of the WASSR array
        disp('     Separating WASSR reference image from WASSR array')
    
        W_ref_freq = -8000; % WASSR reference scan frequency Varian  
        W_offsetfqs_tmp = offsetfqs;
        W_I_tmp = I;

        ind = W_offsetfqs_tmp==W_ref_freq;   
        
        if max(ind) == 0           
           error('No WASSR reference scan found at offset frequency = %d Hz',...
                          W_ref_freq)
        end
        
        I_ref = W_I_tmp(:,:,ind);
        offsetfqs = W_offsetfqs_tmp(~ind);
        I = W_I_tmp(:,:,~ind);
        
    elseif strcmp(type,'CEST')
        
    else
        error('Type not understood')
    end
    
    
%% Determine ppm.
    
    disp('     Calculating ppm')

    ppm(1,:) = offsetfqs/H1fq;

    
%% Order images by ppm.

    disp('     Reordering images')

    [ppm_ord,index] = sort(ppm);
    I_ord = I(:,:,index);
    

%% Calculate normalised pixel intensities

    disp('     Normalizing images')
    
    tmp_ref = I_ref;
    tmp_2Dref = repmat(tmp_ref,[1,1,length(ppm_ord)]);
    I_norm = I_ord./tmp_2Dref*100;
    I_norm(isnan(I_norm)) = 0;
    
    


