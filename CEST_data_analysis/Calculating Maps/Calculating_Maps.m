function [SFmap,IntMap,IntWASSRMap]=Calculating_Maps(CESTdata_cor, WASSRdata_cor,freq,freq_WASSR,SFfreq, rangeInt,rangeWASSRInt )

%% **** Calculate Single Frequency maps, Integral maps, APT maps, NOE maps **** %%

%INPUT:
%      CESTdata_cor = corrected CEST data
%                     d1 x d2  x nfreqNewC Array 
%                     with [d1,d2] the dimensions of each image
%                          
%                          nfreqNewC the number of CEST offset frequencies
%
%      freq = offset frequencies after correction of the CEST data
%             (nfreqNewC vector)
%                 
%      SFfreq = frequency for which the Single Frequency map is calculated
%               If you don't want to calculate a Single Frequency map, put
%               SFfreq=[];                   
%
%      rangeInt = max and min ppm values to integrate between to calculate an
%                 Integral map.
%                 e.g. [0.5,1.5]
%                 If you don't want to calculate an Integral map, put rangeInt=[];     
%



%OUTPUT:
%      SFmap = Structure containing the following fields :
%                    - SFmap.Map = Single Frequency CEST maps
%                                  (d1 x d2  array); 
%                    - SFmap.Freq = same thing than SFfreq (input) 
%
%      IntMap = Structure containing the following fields :
%                    - IntMap.Map = Integral CEST maps
%                                   (d1 x d2  array); 
%                    - IntMap.FreqRange = same thing than rangeInt (input)  
%


%% Dimensions

[Nx,Ny,unused] = size(CESTdata_cor);


%% Calculating Single Frequency CEST Maps

if ~isempty(SFfreq)
    
    % Check that the offset frequencies are symmetric     
    m = 1;
     for k = 1:length(freq)
        if ~isempty(find(freq(k)+freq == 0))
           idx(m) = k;
           m = m+1;
        end
    end
   
    
    % Calculating  MTRasym spectra per each pixel:
    C = CESTdata_cor(:,:,idx);   
    F = freq(idx);
    disp('Calculating  MTRasym spectra per each pixel ')
    
    
        MTRZpix(:,:,:) = (flipdim(C(:,:,F<=0),3)...
                -C(:,:,F>=0))./flipdim(C(:,:,F<=0),3)*100;
        MTRZpix(isnan(MTRZpix)) = 0;
    
    
    MTRfreq = F(F>=0);

    % Calculating Maps    
    disp('     Calculating Single Frequency maps ')
   
    
        [unused,W] = min(abs(MTRfreq-SFfreq)); 
        SF_cest_map(:,:) = medfilt2(MTRZpix(:,:,W(1)));
    
   
    SFmap.Map = SF_cest_map;
    SFmap.Freq = SFfreq;
 
else
    SFmap=[];
end
 
%% Calculating Integral maps

if  ~isempty(rangeInt)
   disp('     Calculating Integral maps')
   
      
       integral(:,:) = integrals (CESTdata_cor(:,:,:) ,rangeInt ,freq);
       integral(:,:) = medfilt2(integral(:,:));
   

   IntMap.Map = integral;
   IntMap.FreqRange = rangeInt;

else
    IntMap=[];
end
%% Calculating WASSR Integral maps

if  ~isempty(rangeWASSRInt) && ~isempty(WASSRdata_cor) 
   disp('     Calculating WASSR Integral maps')
   
     
       WASSRintegral(:,:) = integralsWASSR (CESTdata_cor(:,:,:) ,WASSRdata_cor(:,:,:),rangeWASSRInt,freq,freq_WASSR);
       WASSRintegral(:,:) = medfilt2(WASSRintegral(:,:));
   

   IntWASSRMap.Map = WASSRintegral;
   IntWASSRMap.FreqRange = rangeWASSRInt;

else
    IntWASSRMap=[];
end


