

function [integral] = integrals (Zspec,range ,ppm_ord)
%% **** Calculating Integral maps **** %
%
% INPUT:
%      Zspec = CEST data for one slice
%              d1 x d2 x nfreq Array with [d1,d2] the dimensions of the image
%                                         nfreq the number of offset frequencies
%
%      range = max and min ppm values to integrate between
%             e.g. [0.5,1.5]
% 
%      ppm_ord = offset frequencies (in ppm) by increasing order for the CEST
%                data (nfreq vector)

%OUPUT:
%     integral = Integral map for one slice (d1 x d2 matrix)



    %% Integral CEST effect
    
    Nx = size(Zspec(:,:,1),1);
    Ny = size(Zspec(:,:,1),2);
        
    spacing = (max(ppm_ord)-min(ppm_ord))/length(ppm_ord);
    off_range = logical((ppm_ord <=-range(1)+0.01).*(ppm_ord >= -range(2)-0.01));
    I_norm_reshape = reshape(Zspec(:,:,off_range),Nx*Ny,length(find(off_range)));
    peak_int_off = trapz(I_norm_reshape')*spacing;
    
    on_range = logical((ppm_ord>=range(1)-0.01).*(ppm_ord<=range(2)+0.01));
    I_norm_reshape2 = reshape(Zspec(:,:,on_range),Nx*Ny,length(find(on_range)));
    peak_int_on  = trapz(I_norm_reshape2')*spacing;
    
    integral = (peak_int_off-peak_int_on)./peak_int_off*100;
    integral = reshape(integral,Nx,Ny);
    
        %% Integral CEST effect
%     
%     Nx = size(Zspec(:,:,1),1);
%     Ny = size(Zspec(:,:,1),2);
%         
%     spacing = (max(ppm_ord)-min(ppm_ord))/length(ppm_ord);
%     off_range = logical((ppm_ord <=-range(1)+0.01).*(ppm_ord >= -range(2)-0.01));
%     I_norm_reshape = reshape(Zspec(:,:,off_range),Nx*Ny,length(find(off_range)));
%     peak_int_off = 1 - (trapz(I_norm_reshape')*spacing);
%     
%     on_range = logical((ppm_ord>=range(1)-0.01).*(ppm_ord<=range(2)+0.01));
%     I_norm_reshape2 = reshape(Zspec(:,:,on_range),Nx*Ny,length(find(on_range)));
%     peak_int_on  = 1 - (trapz(I_norm_reshape2')*spacing);
%     
%     integral = (peak_int_off-peak_int_on)./peak_int_off*100;
%     integral = reshape(integral,Nx,Ny);
    
end
        

   