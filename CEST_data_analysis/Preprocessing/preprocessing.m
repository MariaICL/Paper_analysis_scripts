function [ThresholdMask,data_int,ppm_int,idx] = preprocessing(I_ref,data,freq,...
                                              ThresholdMask,smooth,Interp) 

%% **** Preprocessing CEST or WASSR data **** %%

%INPUT: 
%       I_ref = Reference images (d1 x d2 x d3 array)
%
%       data = CEST or WASSR data you wanna process.
%              (d1 x d2  x nfreq) Array with [d1,d2] the dimensions of each image
%                                                 
%                                                 nfreq the number of offset frequencies
%
%       freq = offset frequencies (nfreq vector)
%
%       ThresholdMask = Threshold Mask (d1 x d2 x d3 array).
%                       If you wanna use a pre-existing Threshold Mask,
%                       put it here. 
%                       Otherwise, put ThresholdMask=[];
%
%      smooth = N>0 smoothes the Zspectrum in each pixel by averaging each
%               element along with the 2*N elements at its sides.
%               If you don't want to smooth the spectra, put N=0.
%
%     Interp = Structure containing the following information about interpolation:
%                    -Interp.apply = 'Yes' if you wanna interpolate your data,
%                                    'No' otherwise.
%                    -Interp.NpInt = Number of points you wanna use to 
%                                    interpolate the spectra.
%                    -Interp.type = type of interpolation
%                                   Must be'Spline' or 'Cubic'.
%

%OUTPUT:
%      ThresholdMask = Threshold mask (d1 x d2 x d3 array)
%
%      data_new = preprocessed data
%                 (d1 x d2 x d3 x nfreqNew array)
%
%      ppm_new = new offset frequencies (nfreqNew vector)
%                If you choose to smooth without interpolating, no frequency
%                is added so ppm_new = freq.
%
%      idx = ppm_new is constructed by adding points between the values
%            contained in freq. idx is a vector which contains the
%            index of the ppm_new values corresponding to the freq
%            values (nfreq vector).
%      NB: if ppm_new = freq, idx = 1:length(freq).



%% Threshold mask

[Nx,Ny,nfreq]  = size(data);  
if isempty(ThresholdMask)     
        
            % Noise calculation
            disp('     Calculating the noise');   
            fh = figure;
            imagesc(I_ref(:,:))
            colormap(gray)
            axis off
            title('Select a region where the noise is important','FontSize',16,'FontWeight','Bold')
            Noise_mask = imcrop;
            index = Noise_mask >0;
            Noise = Noise_mask(index);
            N = std(Noise) ;
            close(fh)

            % Thresholding
            fprintf('     Thresholding images\n')
            [Th,unused] = threshold_gui(I_ref(:,:),N);
            ThresholdMask_temp = repmat(Th,[1,1]);
            
            % Masking the brain
            disp('Select the area of interest');
            f_BW = figure;
            BW_display=mat2gray(I_ref(:,:,1));
            BW = roipoly((BW_display));
            close(f_BW)
            BWMask = repmat(BW,[1,1]);
            ThresholdMask = BWMask.*ThresholdMask_temp;
        
end    
    dataTh = data.*repmat(ThresholdMask,[1,1,nfreq]);
    
    
 
%% Smoothing data

if smooth ~= 0
    
   data_S = zeros(size(dataTh));  
   
       fprintf('     Smoothing data \n');
       [row,col] = find(dataTh(:,:,1)>0);   

       for i = 1:length(row)
           m = squeeze(dataTh(row(i),col(i),:));
           mS = moving_average(m,smooth);
           data_S(row(i),col(i),:) = mS;
       end      
   
   
else
    
    data_S = dataTh;
    
end

%% Interpolating data

if strcmp(Interp.apply,'No')
    
    data_int = data_S;
    ppm_int = freq;
    idx = 1:length(freq);
    
elseif strcmp(Interp.apply,'Yes')
    
    % Calulating a new vector of frequencies which contains the original
    % frequencies plus intermediate frequencies:    
    ppm_int=[];
    N = Interp.NpInt-length(freq);
    Nint = length(freq)-1;
    Np = ceil(N/Nint);

    for i = 1:length(freq)-1
        b = linspace(freq(i),freq(i+1),2+Np);
        ppm_int = [ppm_int,b];
    end
    
    ppm_int = unique(ppm_int);
    idx = 1:Np+1:length(ppm_int);
    
    % Interpolating data :
    data_int = zeros(Nx,Ny,length(ppm_int));
    
    
    
        fprintf('     Interpolating data \n')
        
        [row,col] = find(dataTh(:,:,1)>0); 
        
        for i = 1:length(row)
            
            m = squeeze(data_S(row(i),col(i),:));
            if strcmp(Interp.type,'Cubic')
               mInt = interp1(freq,m,ppm_int,'cubic');
            elseif strcmp(Interp.type,'Spline')
                mInt = interp1(freq,m,ppm_int,'spline');
            else
                error('Interp.type not understood')
            end
            data_int(row(i),col(i),:) = mInt;
            
        end
        
    
    
else
    error('Interp.apply not understood')
end

