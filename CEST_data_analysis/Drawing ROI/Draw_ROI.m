function [phantom_mask,poly,ROInames] = Draw_ROI(I_ref,ThresholdMask,nROI,ROInames)

%% **** Drawing ROI using the reference images **** %%

%INPUT: 
%      I_ref = reference images (d1 x d2 array)
%
%      ThresholdMask = Threshold Mask (d1xd2 array).
%
%      nROI = number of ROI
%
%      ROInames = character array containing the names of the ROIs.
%                 eg: ROInames = char('ROI1','ROI2')
%                 ROInames must either contain nROI elements, either be
%                 empty. In that case, the default names used are 'ROI1',
%                 'ROI2'...

%OUTPUT
%      phantom_mask = ROI mask(s)( d1 x d2 x nROI array) 
%
%      poly = structure containing the coordinates of the points which define the ROI(s):
%                    poly.X = NpointsPerRoi x nROI x d3 array
%                             For each slice there is a matrix with nROI columns, 
%                             each colum contains the x-coordinates of the corresponding ROI. 
%                    poly.Y = NpointsPerRoi x nROI x d3 array
%                             For each slice there is a matrix with nROI columns, 
%                             each colum contains the y-coordinates of the corresponding ROI. 


%% Check ROInames

if size(ROInames,1)~=nROI && ~isempty(ROInames)
   error('ROInames must contains nROI elements')
end

if isempty(ROInames)
   ROInames=strcat('ROI',num2str((1:nROI)'));
end

%% Drawing ROIs

%[Nx,Ny] = size(Iref{1,1});
%phantom_mask = zeros(Nx,Ny,nROI{1,1});

[Nx,Ny] = size(I_ref);
phantom_mask = zeros(Nx,Ny,nROI);
 
   
     % Displaying the thresholded reference image
     fprintf('     User input: select the ROIs\n');
     roi = figure; 
     imagesc(I_ref(:,:).*ThresholdMask(:,:))
     %imagesc(Iref{1,1}(:,:).*Masks{1,1}.ThresholdMask(:,:))
     colormap gray
     axis off
     axis image
     hold on;

    % Creating the ROI
    xpoly = cell(1,nROI);
    ypoly = cell(1,nROI);
    Nelem = zeros(nROI); 
    %xpoly = cell(1,nROI{1,1});
    %ypoly = cell(1,nROI{1,1});
    %Nelem = zeros(nROI{1,1}); 
    for k = 1:nROI
        name = sprintf('Select ROI %d ( %s )',k, ROInames(k,:));
        title(name,'FontSize',16,'FontWeight','Bold');
        [phantom_mask(:,:,k), xpoly{k}, ypoly{k}] = roipoly;
        Nelem(k) = length(xpoly{k});
        plot(xpoly{k},ypoly{k},'LineWidth',2,'Color','r')
        hold on
    end
    
       
    % Storing the coordinates of the ROIs (to display them later)
    Xpoly(:,:) = nan(max(max(Nelem)),nROI);
    Ypoly(:,:) = nan(max(max(Nelem)),nROI);
    for k = 1:nROI
        Xpoly(1:Nelem(k),k,1) = xpoly{k};
        Ypoly(1:Nelem(k),k,1) = ypoly{k};
    end
    
    


poly.X=Xpoly;
poly.Y=Ypoly;

close(roi)

end



