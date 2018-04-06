clear%% **** Parameters for the CEST data processing **** %%

clear all,
close all,
clc,

addpath('B0 correction','Calculating Spectra', 'Calculating Maps',...
        'Loading','Preprocessing', 'Whole CEST process','Drawing ROI',...
        'Displaying');
    

%% Scan parameters

  
       refscan =  {5};                            
       CESTscan = {6};
       Wscan =    {7}; 
                           
                              
                             
%% Preprocessing

    % Thresholding
    
       
          SameThExp='No';  % 'Yes' if you wanna use the same threshold mask for all experiments,
                            % 'No' otherwise.                      
                             

    % CEST data preprocessing                 
    smooth_CEST = 1; % smooth_CEST = N>0 smoothes the Zspectrum in each pixel by averaging each
                     % element along with the 2*N elements at its sides.
                     % If you don't want to smooth the spectra, put N=0.
                     
                     
    InterpC.apply = 'Yes';     % 'Yes' if you wanna perform CEST data interpolation, 'No' otherwise                      
    InterpC.NpInt = 1000;      % Number of points used to interpolate the CEST data
    InterpC.type = 'Spline';   % Type of interpolation : 'Spline' or 'Cubic'
    
                  
   % WASSR data preprocessing                     
    smooth_WASSR = 1;          % same thing as smooth_CEST but for the WASSR data
                     
    InterpW.apply= 'Yes';      % 'Yes' if you wanna perform WASSR data interpolation, 'No' otherwise                     
    InterpW.NpInt= 1000;       % Number of points to use to interpolate the WASSR data
    InterpW.type= 'Spline';    % Type of interpolation : 'Spline' or 'Cubic'
    
     
%% Drawing ROI(s)   
   
   DrawRoi = 'No'; % 'Yes' if you want to perform a ROI analysis
                    % 'No' otherwise;
                    

                   
           % If  DrawRoi = 'Yes', fill the following parameters:
           SameRoiExp='Yes';   % 'Yes' if you wanna use the same ROI for all
                               % experiments,
                               % 'No' otherwise.
                               
           RefRoi = 'No';       %'Yes' load an specific reference for drawing rois
                                %'No'use Iref
           
           nROI = {2};            % Number of ROI for each experiment (nExp cell array)
                                  % If you wanna use the same number of ROI for
                                  % all experiments, just enter the number once.
                                  % NB: it is assumed that the number of ROI is the 
                                  % same for all slices. 

           ROInames = {char()};     % Cell array containing the names of the 
                                    % ROIs for each experiment.  
                                    % ROInames{i} must either contain
                                    % nROI{i} names, either be empty.
                                    % In that case, the default 
                                    % names used are'ROI1', 'ROI2' ...
                                    % If you wanna use the same names for
                                    % all experiments, just enter the names
                                    % once.
                                    % NB: it is assumed that the ROI names are the 
                                    % same for each slice. 

    
%% B0 Correction

    b0cor = 'MaxSymmetry_WASSR';  % Type of B0correction:
                               % 'No' if you don't want a correction,
                               % 'LoadMap' if you have already a b0map,
                               % otherwise choose one of the following types:
                               % 'SplineWASSR', 'Lorentz', 'MaxSymmetry_WASSR',
                               % 'MaxSymmetry_CEST' or 'SplineCEST'.
                               % (see B0Mapping.m for details) 
    
        
%% B0 Correction for WASSR

    b0corWASSR = 'SplineWASSR';  % Type of B0correction:
                               % 'No' if you don't want a correction,
                               % otherwise choose one of the following types:
                               % 'SplineWASSR', 'Lorentz', 'MaxSymmetry_WASSR',
                               % (see B0Mapping.m for details) 
  

                               
%% Calculating Maps

    %Single Frequency maps
    
    SFfreq = {0.6};   % nExp cell array.  
                   % The cell i contains the frequency used to calulate the SF map
                   % for the experiment i.
                   % If you wanna use the same frequency for
                   % each experiment, just enter the frequency once.
                   % If you don't wanna callculate a Single Frequency
                   % CEST Map, put SFfreq={}.
                   % NB: it is assumed that the frequency is the 
                   % same for all slices. 
                                                

                     
    % Integral maps

    rangeInt= {[0.5 1]};   % nExp cell array.  
                    % The cell i contains the frequency range (ie a 2-elements vector)
                    % used to calulate the Integral map for the experiment i.
                    % If you wanna use the same range for
                    % each experiment, just enter the range once.
                    % If you don't wanna calculate a Integral
                    % CEST Map, put rangeInt={}.
                    % NB: it is assumed that the range is the 
                    % same for each slice. 
                        
                    
    % WASSR Integral maps

    rangeWASSRInt= {[0.5 1]};   % nExp cell array.  
                    % The cell i contains the frequency range (ie a 2-elements vector)
                    % used to calulate the Integral map for the experiment i.
                    % If you wanna use the same range for
                    % each experiment, just enter the range once.
                    % If you don't wanna calculate a Integral
                    % WASSR Map, put rangeInt={}.
                    % NB: it is assumed that the range is the 
                    % same for each slice.                 
                        
%% Saving

    filename = '';  % All data are saved in a MATLAB 
                        % file (MAT-file) called filename.
                          
%% Displaying

    % Displaying ROI
        disROI='Yes';     % 'Yes' if you wanna display the reference images with the ROIs,
                         % 'No' otherwise. 
                      
   % Displaying B0 maps                                        
        disB0.is='Yes';      % 'Yes' if you wanna display the B0maps for each experiment,
                            % 'No' otherwise. 
        disB0.range=[-0.2 0.2]; % disB0.range=[a,b] normalizes the values in the B0maps to the range [a,b];
                            % disB1.range='auto' for a display by default
                            
 
    % Displaying SF maps                      
        disSFmap.is = 'Yes';           % 'Yes' if you wanna display the SF maps for each experiment,
                                      % 'No' otherwise.
        disSFmap.range = [0 20];      % disSFmap.range=[a,b] normalizes 
                                      % the values in the SF CEST maps
                                      % to the range [a,b];
                                      % disSFmap.rangeInt='auto' for a display by default  
 
    % Displaying Integral maps                     
        disIntMap.is = 'Yes';           % 'Yes' if you wanna display the Integral maps for each experiment,
                                       % 'No' otherwise.
        disIntMap.range = [-5 15];      % same thing as disSFmap.range but for 
        
    % Displaying WASSR Integral maps                     
        disIntWASSRMap.is = 'Yes';           % 'Yes' if you wanna display the WASSR Integral maps for each experiment,
                                       % 'No' otherwise.
        disIntWASSRMap.range =[25 35];      % same thing as disSFmap.range but for 
                                       % the Integral maps                                   % the Integral maps
                                       
                                       
         
    % Displaying MTR and Z Spectra
        disSpectra='OnePerExp_AllRoi'; % 'OnePerExp_OnePerRoi' if you wanna 
                                       % display one figure per each experiment
                                       % and per each ROI
                                       % 'OnePerExp_AllRoi' if you wanna 
                                       % display one figure per each experiment
                                       % (all ROI on the same figure)
                                       % 'AllExp_OnePerRoi' if you wanna 
                                       % display one figure per each ROI
                                       % (all experiments on the same figure)
                                       % 'No' if you don't wanna display
                                       % the spectra.
                                       
     % Displaying WASSR Z Spectra
        disWSpectra='No';
       
%% Processing
    mainCEST
