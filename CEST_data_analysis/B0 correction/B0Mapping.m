function [ b0map]= B0Mapping(WassrTh,freq_wassr,CestTh,freq_cest,type)

%% **** B0mapping based on WASSR or CEST data set **** %%

% INPUT:
%       WassrTh = preprocessed WASSR data
%                 d1 x d2 x nfreqNewW Array 
%                 with [d1,d2] the dimensions of each image
%                        
%                        nfreqNewW the number of WASSR offset frequencies
%           
%       freq_wassr = offset frequencies (in ppm) by increasing order for the WASSR data
%                   (nfreqNewW vector)
%
%       CestTh = preprocessed CEST data
%                d1 x d2  x nfreqNewC Array 
%                with [d1,d2] the dimensions of each image
%                      
%                      nfreqNewC the number of CEST offset frequencies
%           
%       freq_cest = offset frequencies (in ppm) by increasing order for the CEST data
%                   (nfreqNewC vector)
%
%       type = Type of B0 correction:
%               'SplineWASSR' : Find the minimum of the WASSR Z-Spectrum
%                               for each pixel.
%               'SplineCEST' : Find the minimum of the CEST Z-Spectrum
%                               for each pixel.
%               'Lorentz': Try to fit a Lorentzian curve with the WASSR
%                          Z-Spectrum for each pixel (in Least Squares sense)
%                          and retrieve the minimum of the Lorentzian fitted curve.
%               'MaxSymmetryWASSR' : Maximum Symmetry Algorithm applied to
%                                    the WASSR spectra
%               'MaxSymmetryCEST' : Maximum Symmetry Algorithm applied to
%                                    a little portion of the CEST spectra


%OUTPUT:
%       b0map = d1 x d2  array
%               Contains the center frequency of the Z_spectrum for
%               each pixel .



%% Switching the methods

    %Calculating the b0map for the spline method
   
if strcmp(type,'SplineWASSR')
    [unused,unused,unused] = size(WassrTh);
    
        fprintf('     Start calculating B0map( SplineWASSR method) ')
        b0map(:,:) = B0map_Spline(WassrTh(:,:,:), freq_wassr);
        fprintf('     Spline is done \n') 
        
elseif strcmp(type,'SplineCEST')
    [unused,unused,unused] = size(CestTh);
    
        fprintf('     Start calculating B0map( SplineCEST method)')
        b0map(:,:) = B0map_Spline(CestTh(:,:,:), freq_cest);
        fprintf('     Spline is done \n') 
    
    
     %Calculating the b0map for the Lorentz method
elseif strcmp(type,'Lorentz')
    [unused,unused,unused] = size(WassrTh);
    
        fprintf('     Start calculating B0map( Lorentz method) ')
        b0map(:,:) = B0map_LZ_N(WassrTh(:,:,:), freq_wassr);
        fprintf('     Lorentz is done \n') 
    
       
      %Calculating the b0map for the Maximum Symmetry method       
elseif strcmp(type,'MaxSymmetry_WASSR')
    [unused,unused,unused] = size(WassrTh);
    
        fprintf('     Start calculating B0map( MaxSymmetry WASSR method) ')
        b0map(:,:) = B0map_MaxSym(WassrTh(:,:,:), freq_wassr,type);
        fprintf('     MaxSymmetry is done \n') 
    

elseif strcmp(type,'MaxSymmetry_CEST')
    [unused,unused,unused] = size(CestTh);
    
        fprintf('     Start calculating B0map( MaxSymmetry CEST method) ') 
        b0map(:,:) = B0map_MaxSym(CestTh(:,:,:), freq_cest,type);
        fprintf('     MaxSymmetry is done \n') 
    
       
else
    error('Type of correction not understood.\n')

end


function [b0map] = B0map_Spline(data, freq) 

    [M, N, siw] = size(data);
    b0map = zeros( M,N);
    [row, col] = find (data(:,:,1)>0); 

    for ind = 1: length(row)
        
        mSpb0 = squeeze(data(row(ind),col(ind),:));
        [unused,I] = min(mSpb0);
        b0shift = freq(I(1));
        b0map(row(ind), col(ind)) = b0shift;

    end

function [b0map] = B0map_MaxSym(data, freq,type)

    [Wx,Wy,nfreq] =size(data);

  % Determine search initialisation point for each pixel 
    disp('     Determining search initialisation point')
    dataR=reshape(data,Wx*Wy, nfreq);
    [unused,IdxMin]= min(dataR,[],2); 
    search_point_idx = reshape(IdxMin,Wx,Wy);
    search_point = freq(search_point_idx);  
  
  % Search for minimum mean squared error between original and mirrored Z-spectra    
        % Opening matalbpool for parallel computing
     parpool;
    
        % Launching the algorithm
    disp('     Searching for Maximum Symmetry centre frequencies')
    
    progressbar;
    progressbar('Searching for maximum symmetry centre frequencies');
    bar = 0;
    
    options = optimset('Display','off','TolX',0.001,'TolFun',0.001,'MaxIter',100,...
        'MaxFunEvals',100);
    
    b0map = zeros(Wx, Wy);

    for b = 1:Wx

          parfor c = 1:Wy
               m = squeeze(data(b,c,:));
               
               if strcmp(type,'MaxSymmetry_CEST')
                  th = min(m)+(max(m)-min(m))*(5/100);
                  idx = (m<=th);
                  fxi = m(idx);
                  w = freq(idx);
               else
                  fxi = m;
                  w = freq;
               end
                          
               if max(fxi) > 0
                    mscf_tmp(1,c) =...
                       fminsearch(@(C)sum((fxi-interp1(-w(:),...
                       fxi,w(:)-2*C,'spline',100)).^2),...
                       search_point(b,c),options);                
               else                  
                    mscf_tmp(1,c) = 0;
               end
      
          end

          b0map(b,:) = mscf_tmp;
          clear mscf_tmp
            
          bar = bar+1;
          progressbar(bar/sum(Wx));
            
    end
    progressbar(1);

    
function [b0map] = B0map_LZ_N(data, freq) 

    [M, N, siw] = size(data);
    b0map = zeros(M,N); 
    
    progressbar;
    progressbar('Calculating B0map(Lorentz method)');
    bar = 0;
    
for b = 1:M
    
    parfor c = 1:N
        
    % Searching for the initial parameters to launch the algorithm    
        mb0 = squeeze (data(b, c, : ));
        
        % Find par1
        [unused,yy0] = min(mb0);
        par1 = freq(yy0);
        
        % Find par3
        th = min(mb0)+(max(mb0)-min(mb0))*(50/100);
        idx = find(mb0<=th);
        f = freq(idx);
        par3 =(f(end)-f(1))/2;
        
        % Find par4
        par4 = min(mb0);
        
        % Find par2
        par2 = max(mb0)-par4;
    
    % Launching the algorithm    
        par0 = [par1,  par2, par3, par4];

        lb = [par1-max(freq), 40, 1e-3, 0 ];

        ub = [par1+max(freq), 150, 5, 50];

        options = optimset('MaxFunEvals',1000000,'TolFun',1e-10,'TolX',1e-10,  'Display',  'off' );

        par(:,c) = lsqcurvefit(@lorentz_iN,par0(:) , freq(:) , mb0(:) , lb, ub, options); 

        
    end
    
        b0map(b, :) = par(1,:);

        bar=bar+1;
        progressbar(bar/sum(M));
end

progressbar(1);


function y_fit =lorentz_iN(par, delta)

    denum=1+(par(3)./(delta-par(1))).^2;
    y_fit=par(4)+par(2)./denum;

