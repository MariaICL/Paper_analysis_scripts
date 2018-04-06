function [WZspec, Zspec, AsymSpec, data_cor,wassr_cor] = Calculating_Spectra(CESTdata,ppm_cest, ppm_int ,idx, ...
                                     phantom_mask, b0map ,b0mapWASSR,WASSRdata, ppm_int_WASSR )

%% **** Calculate Z-spectra and MTRasym spectra  **** %%

%INPUT:
%      CESTdata = preprocessed CEST data
%                 d1 x d2  x nfreqNewC Array 
%                 whith [d1,d2] the dimensions of each image
%                       
%                       nfreqNewC the number of CEST offset frequencies
%
%      ppm_cest = original offset frequencies (in ppm) by increasing order (nfreqOr vector)
%
%      ppm_int = offset frequencies after preprocessing of the CEST data
%               (nfreqNewC vector)
%      NB: if you did not interpolate your CEST data before, put
%      ppm_int = ppm_cest.
%
%      idx = ppm_int is constructed by adding points between the values
%            contained in ppm_cest. idx is a vector which contains the
%            index of the ppm_int values corresponding to the ppm_cest
%            values(nfreqOr vector).
%
%      phantom_mask = ROI mask(s)( d1 x d2  x nROI array)
%
%
%      b0map = d1 x d2  array
%              Contains the center frequency of the Z-spectrum for
%              each pixel .
%              If you don't want to apply a B0correction, put b0map=[];
%
%      b1map = d1 x d2  array
%              Contains the estimated value of the RF field for each pixel
%              
%              If you don't want to apply a B1correction, put b1map=[];
%                        

%OUTPUT:
%      Zspec = structure containing 4 fields:
%                   -Zspec.IntValues = interpolated Z-spectra values for
%                                      each ROI 
%                                      nfreqNew x nROI 
%                                      nfreqNew > nfreq if the data are interpolated
%
%                   -Zspec.IntFreq = same as the input ppm_int
%
%                   -Zspec.OriginalValues = Z -spectra values for each ROI 
%                                           at the original frequencies
%                                           nfreqOr x nROI
%                                  
%                   -Zspec.OriginalFreq = same as the ppm_cest input
%
%
%      AsymSpec = same structure than Zspec but with the MTRasym
%                 spectra values

                     
%% Number of ROI

nROI = size(phantom_mask,3);

%% B0 correction 

if ~isempty(b0map)
   disp('     Applying the B0map') 
   data_cor = B0correction_Spline(CESTdata, ppm_int, b0map);
   if ~isempty(b0mapWASSR)
       
       wassr_cor = B0correction_Spline(WASSRdata, ppm_int_WASSR, b0mapWASSR);
   else
       wassr_cor=[];
   end
else
   data_cor = CESTdata;  
   wassr_cor = WASSRdata;
end


%% Calculating Z-spectra per each ROI

disp('     Calculating Z-spectra')

[unused,unused,nfreq] = size(data_cor);
Zroi = zeros(nfreq,nROI);


    
    for k = 1:nROI

        for n = 1:nfreq
            ROI = data_cor(:,:, n).* phantom_mask(:,:,k);
            Zroi(n,k) = mean(nonzeros(ROI));
        end

    end
    


    %Store the Z-spectra interpolated values  
Zspec.IntValues = Zroi;
Zspec.IntFreq = ppm_int;

   %Retrieve and store the values of the Z-spectra at the original frequencies
Zspec.OriginalValues = Zroi(idx,:,:); %???
Zspec.OriginalFreq = ppm_cest;
   

%% Calculating WASSR Z-spectra per each ROI

if ~isempty(wassr_cor)

    disp('     Calculating WASSR Z-spectra')

    [unused,unused,Wnfreq] = size(wassr_cor);
    WZroi = zeros(Wnfreq,nROI);

    

        for k = 1:nROI

            for n = 1:Wnfreq
                ROI = wassr_cor(:,:, n).* phantom_mask(:,:,k);
                WZroi(n,k) = mean(nonzeros(ROI));
            end

        end

   

        %Store the Z-spectra interpolated values  
    WZspec.IntValues = WZroi;
    WZspec.IntFreq = ppm_int_WASSR;
    
else
    WZspec=[];
   

end

%% Calculating  MTRasym spectra per each ROI:

 disp('     Calculating MTRasym spectra')
 
    % Check that the offset frequencies are symmetric
 m = 1;
 for k = 1:length(ppm_cest)
     if ~isempty(find(ppm_cest(k)+ppm_cest == 0))
        idx_ppm(m) = k;
        m = m+1;
     end
 end
 
    % Calculate the MTRasym values 
 Z = Zspec.OriginalValues(idx_ppm,:,:);
 PPM = ppm_cest(idx_ppm);


     MTRZroi(:,:) = (flipdim(Z(PPM<=0,:),1)-Z(PPM>=0,:))./flipdim(Z(PPM<=0,:),1)*100;

 MTRw = PPM(PPM>=0);
 
     % Store the MTRasym spectra values at the original frequencies
 AsymSpec.OriginalValues = MTRZroi;
 AsymSpec.OriginalFreq = MTRw;
 
    % Interpolate and store the MTRasym spectra values 
     
      AsymSpec.IntValues(:,:) = interp1(MTRw,MTRZroi(:,:),ppm_int(ppm_int>=0)','spline');
  
  AsymSpec.IntFreq = ppm_int(ppm_int>=0);

%%%% end of main function  
          
function [newMI] = B0correction_Spline(Mz, ppm_int, B0map) 


    [Nx, Ny, nfreq] = size(Mz);
    newMI = zeros(size(Mz));

    
        [row, col] = find (Mz(:,:,1)>0);

        for ind = 1: length(row)
            m = squeeze(Mz(row(ind),col(ind),:));        
            b0 = B0map(row(ind), col(ind));
%             b0 = cell2mat(b0); %B0 map

             if b0 <= 0
                mSp = interp1(ppm_int-b0,m,ppm_int,'spline',m(1)); 
             elseif b0 > 0
                 mSp = interp1(ppm_int-b0,m,ppm_int,'spline',m(end));
             end
             
            newMI(row(ind), col(ind),:) = mSp;

        end
   


