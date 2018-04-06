%% **** Compute the different steps of the data processing using the parameters **** %%


%% Loading data       

    
    Nexp = length(refscan);
    
    % CEST data      
    Iref = cell(1,Nexp);
    CESTdata = cell(1,Nexp);
    ppm_cest = cell(1,Nexp);
    
    fprintf('LOADING CEST DATA FOR EXPERIMENT %d\n', 1)    
    [Iref{1},CESTdata{1},ppm_cest{1},scandirectory] = loadVarianData(refscan{1},CESTscan{1},'CEST');    
    
    for i=2:Nexp
        
        fprintf('LOADING CEST DATA FOR EXPERIMENT %d\n', i)
        [Iref{i},CESTdata{i},ppm_cest{i}] = loadVarianData(refscan{i},CESTscan{i},'CEST',scandirectory);
        
    end
    for i=1:Nexp
        CESTdata{i} = double(CESTdata{i});
    end
    
    % WASSR data
    
        
        WASSRdata = cell(1,Nexp);
        ppm_wassr = cell(1,Nexp);
        
        for i=1:Nexp
            fprintf('LOADING WASSR DATA FOR EXPERIMENT %d\n', i)
            [unused,WASSRdata{i}, ppm_wassr{i}] = loadVarianData(refscan{i},Wscan{i},'WASSR', scandirectory);
        end
        
        if size(WASSRdata,1) ~= size(CESTdata,1) || size(WASSRdata,2) ~= size(CESTdata,2) || size(WASSRdata,3) ~= size(CESTdata,3)
           error ('The WASSR data and CEST data do not have the same dimensions')
        end
        for i=1:Nexp
        WASSRdata{i} = double(WASSRdata{i});
        end        

%% Preprocessing
    
    % CEST data
    if strcmp(SameThExp,'No')

        for i = 1:Nexp
            fprintf('PREPROCESSING CEST DATA ( EXPERIMENT %d )\n', i) 
            [ThMask{i},CEST_new{i},ppm_cest_new{i},idx{i}] = preprocessing(Iref{i},CESTdata{i},ppm_cest{i},...
                                                             [],smooth_CEST,InterpC); 
        end

    elseif strcmp(SameThExp,'Yes')

         fprintf('PREPROCESSING CEST DATA ( EXPERIMENT %d )\n', 1)
         [ThMask{1},CEST_new{1},ppm_cest_new{1},idx{1}] = preprocessing(Iref{1},CESTdata{1},ppm_cest{1},...
                                                          [],smooth_CEST,InterpC); 
         for i = 2:Nexp
            fprintf('PREPROCESSING CEST DATA ( EXPERIMENT %d )\n', i) 
            [ThMask{i},CEST_new{i},ppm_cest_new{i},idx{i}] = preprocessing(Iref{i},CESTdata{i},ppm_cest{i},...
                                                             ThMask{1},smooth_CEST,InterpC);
         end

    else
        disp('SameThExp not understood')    
    end
    


 % WASSR data


   for i=1:Nexp   
       fprintf('PREPROCESSING WASSR DATA ( EXPERIMENT %d )\n', i)
       [unused,WASSR_new{i},ppm_wassr_new{i}] = preprocessing(Iref{i},WASSRdata{i},ppm_wassr{i},...
                                                   ThMask{i},smooth_WASSR,InterpW); 
   end



%% Drawing ROI(s)  

if strcmp(DrawRoi,'Yes')

        %Checking the names and the number of ROI
    if length(ROInames) == 1
       ROInames = repmat(ROInames,1,Nexp);
    end

    if length(nROI) == 1
       nROI=repmat(nROI,1,Nexp);
    end

        %Designing the ROIs
    if strcmp(RefRoi,'No')     
        if strcmp(SameRoiExp,'Yes')

           fprintf('DRAWING ROI ( EXPERIMENT %d )\n', 1)
           [RoiMask{1},poly{1},ROInames{1}] = Draw_ROI(Iref{1},ThMask{1},nROI{1},...
                                                       ROInames{1});
           RoiMask = repmat(RoiMask,1,Nexp);
           ROInames = repmat(ROInames,1,Nexp);
           poly = repmat(poly,1,Nexp);

        elseif strcmp(SameRoiExp,'No')

            for i = 1:Nexp 
                fprintf('DRAWING ROI ( EXPERIMENT %d )\n', i) 
                [RoiMask{i},poly{i},ROInames{i}] = Draw_ROI(Iref{i},ThMask{i},nROI{i},...
                                                            ROInames{i});
            end

        else
           error('SameRoiExp not understood ')      
        end
  
    elseif strcmp(RefRoi,'Yes')
        [Iref{1}]=loadimage(1);
        Iref{1} = cell2mat(Iref{1});
        if strcmp(SameRoiExp,'Yes')

           fprintf('DRAWING ROI ( EXPERIMENT %d )\n', 1)
           [RoiMask{1},poly{1},ROInames{1}] = Draw_ROI(Iref{1},ThMask{1},nROI{1},...
                                                       ROInames{1});
           RoiMask = repmat(RoiMask,1,Nexp);
           ROInames = repmat(ROInames,1,Nexp);
           poly = repmat(poly,1,Nexp);

        elseif strcmp(SameRoiExp,'No')

            for i = 1:Nexp 
                fprintf('DRAWING ROI ( EXPERIMENT %d )\n', i) 
                [RoiMask{i},poly{i},ROInames{i}] = Draw_ROI(Iref{i},ThMask{i},nROI{i},...
                                                            ROInames{i});
            end

        else
           error('SameRoiExp not understood ')      
        end
    end
  
elseif strcmp(DrawRoi,'No')
    
    RoiMask=ThMask;
    ROInames=repmat({'ROI1'},1,Nexp);
    poly={};
    
end

%% Calculating B0 map

if strcmp(b0cor,'No') == 0 
    
    b0map=cell(1,Nexp);
    for i = 1:Nexp
        fprintf('CALCULATING B0 MAP( EXPERIMENT %d )\n', i) 
        b0map{i} = B0Mapping(WASSR_new{i},ppm_wassr_new{i},CEST_new{i},ppm_cest_new{i},b0cor);
        poolobj = gcp('nocreate');
        delete(poolobj);
    end
     
   
else
    b0map = cell(Nexp,1);
  
end

%% Calculating B0 map for WASSR

if strcmp(b0corWASSR,'No') == 0
    
    b0mapWASSR=cell(1,Nexp);
    for i = 1:Nexp
        fprintf('CALCULATING B0 WASSR MAP( EXPERIMENT %d )\n', i) 
        b0mapWASSR{i} = B0Mapping(WASSR_new{i},ppm_wassr_new{i},CEST_new{i},ppm_cest_new{i},b0corWASSR);
    end
     
else
    b0mapWASSR = cell(Nexp,1);
  
end


%% Calculating Z-spectra and MTRasym spectra

for i=1:Nexp
    fprintf('CALCULATING SPECTRA( EXPERIMENT %d )\n', i)
    [Wspec{i},Zspec{i}, AsymSpec{i}, CESTdata_cor{i}, WASSRdata_cor{i}] = Calculating_Spectra(CEST_new{i},ppm_cest{i}, ppm_cest_new{i} ,idx{i}, ...
                                     RoiMask{i}, b0map{i} ,b0mapWASSR{i},WASSR_new{i},ppm_wassr_new{i} );
end
                                            
%% Calculating maps    

    % Single Frequency Maps
if length(SFfreq) == 1
   SFfreq = repmat(SFfreq,1,Nexp); 

elseif isempty(SFfreq)
   SFfreq = repmat({[]},1,Nexp);

elseif length(SFfreq) ~= Nexp
    error('SFfreq must contain Nexp elements')

end

    % Integral maps
if length(rangeInt) == 1
   rangeInt = repmat(rangeInt ,1,Nexp); 

elseif isempty(rangeInt)
   rangeInt = repmat({[]},1,Nexp);

elseif length(rangeInt)~=Nexp
    error('rangeInt must contain Nexp elements')

end

    % Integral WASSR maps
if length(rangeWASSRInt) == 1
   rangeWASSRInt = repmat(rangeWASSRInt ,1,Nexp); 

elseif isempty(rangeWASSRInt)
   rangeWASSRInt = repmat({[]},1,Nexp);

elseif length(rangeWASSRInt)~=Nexp
    error('rangeWASSRInt must contain Nexp elements')

end

 
    %Calculating Maps
for i = 1:Nexp
    fprintf('CALCULATING MAPS( EXPERIMENT %d )\n', i)
    [SFmap{i},IntMap{i},IntWASSRMap{i}] = Calculating_Maps(CESTdata_cor{i},WASSRdata_cor{i},...
                                               ppm_cest_new{i},ppm_wassr_new{i},SFfreq{i}, rangeInt{i},rangeWASSRInt{i});
end

%% Restructuring and saving the data

disp('SAVING DATA')

    %Scans 
for i = 1:Nexp
    
       Scan{i}.Type = 'Varian';
       Scan{i}.RefScan = refscan{i};
       Scan{i}.CESTscan = CESTscan{i};
       Scan{i}.WASSRscan = Wscan{i};
       Scan{i}.DirectoryPath = scandirectory;
    
end

    % Raw data (without any threshold mask or correction)      
for i=1:Nexp
    RawData{i}.Iref = Iref{i};
    RawData{i}.CESTdata= CESTdata{i};
    RawData{i}.ppm_cest=ppm_cest{i};
    if strcmp(b0cor,'SplineWASSR') || strcmp(b0cor,'MaxSymWASSR') || strcmp(b0cor,'Lorentz')
       RawData{i}.WASSRdata= WASSRdata{i};
       RawData{i}.ppm_wassr=ppm_wassr{i};
    end
    
end

    % Corrected data (after preprocessing)
for i=1:Nexp
    CorData{i}.CESTcor = CEST_new{i};
    CorData{i}.ppm_cest_new = ppm_cest_new{i};
    if ~isempty(WASSR_new{i})
       CorData{i}.WASSRcor= WASSR_new{i};
       CorData{i}.ppm_wassr_new=ppm_wassr_new{i};
    end
    
end

    % Masks
for i=1:Nexp
    Masks{i}.ThresholdMask = ThMask{i};
    %Masks{i}.RoiMask = RoiMask{i};
    %if ~isempty(poly)
    %   Masks{i}.RoiCoordinates=poly{i};
    %end
    %Masks{i}.NumberRoi = nROI{i};
    %Masks{i}.RoiNames = ROInames{i};
end


    % B0 correction
for i=1:Nexp
    if strcmp(b0cor,'No')==0
       B0{i}.TypeCor = b0cor;
       B0{i}.Map= b0map{i};
    else
       B0{i} = 'NoB0Correction';
    end
end

    % B0 WASSR correction
for i=1:Nexp
    if strcmp(b0corWASSR,'No')==0
       B0WASSR{i}.TypeCor = b0corWASSR;
       B0WASSR{i}.Map= b0mapWASSR{i};
    else
       B0WASSR{i} = 'NoB0WASSRCorrection';
    end
end

    % idx
for i=1:Nexp
    
    idx{i} = idx{i};
    
end

    %Reference image for future ROI display
for i=1:Nexp
   
    Iref{i} =  Iref{i};
   
end
   
    %Saving
save(filename, 'Wspec','Zspec','AsymSpec','SFmap','IntMap','IntWASSRMap','Scan',...
     'RawData','Masks','B0','B0WASSR','CorData','idx','Iref')
    

%% Displaying
 
disp('DISPLAYING DATA' ) 

Display(disWSpectra,disROI,disB0,disSFmap,disIntMap,disIntWASSRMap,disSpectra,B0,Iref,Masks,...
                  Zspec,Wspec,AsymSpec,SFmap,IntMap,IntWASSRMap) 
              
clear all

