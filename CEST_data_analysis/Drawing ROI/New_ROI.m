%% Drawing ROI(s) (Copied from mainCEST.m)
Nexp = length(Masks);
DrawRoi = 'Yes'; 
SameRoiExp='Yes';
nROI = {2}; 
ROInames = {char('','')};

if strcmp(DrawRoi,'Yes')

        %Checking the names and the number of ROI
    if length(ROInames) == 1
       ROInames = repmat(ROInames,1,Nexp);
    end

    if length(nROI) == 1
       nROI=repmat(nROI,1,Nexp);
    end

        %Designing the ROIs
    if strcmp(SameRoiExp,'Yes')

       fprintf('DRAWING ROI ( EXPERIMENT %d )\n', 1)
       [RoiMask{1},poly{1},ROInames{1}] = Draw_ROI(Iref{1,1},Masks{1,1}.ThresholdMask,nROI{1},...
                                                   ROInames{1});
       RoiMask = repmat(RoiMask,1,Nexp);
       ROInames = repmat(ROInames,1,Nexp);
       poly = repmat(poly,1,Nexp);

    elseif strcmp(SameRoiExp,'No')

        for i = 1:Nexp 
            fprintf('DRAWING ROI ( EXPERIMENT %d )\n', i) 
            [RoiMask{i},poly{i},ROInames{i}] = Draw_ROI(Iref{1,i},Masks{1,i}.ThresholdMask,nROI{i},...
                                                        ROInames{i});
        end

    else
       error('SameRoiExp not understood ')      
    end

elseif strcmp(DrawRoi,'No')
    
%     RoiMask=ThMask;
%     ROInames=repmat({'ROI1'},1,Nexp);
%     poly={};
    
end
%% Calculating Z-spectra and MTRasym spectra (Copied from mainCEST.m)
for i=1:Nexp
    fprintf('CALCULATING SPECTRA( EXPERIMENT %d )\n', i)
    [Wspec{i},Zspec{i}, AsymSpec{i}, CESTdata_cor{i}, WASSRdata_cor{i}] = Calculating_Spectra(CorData{1,i}.CESTcor,RawData{1,i}.ppm_cest, CorData{1,i}.ppm_cest_new,idx{i}, ...
                                     RoiMask{i}, B0{1,i}.Map,B0WASSR{1,i}.Map,CorData{1,i}.WASSRcor,CorData{1,i}.ppm_wassr_new );
end
  %% SAVING NEW ROI DATA 
     filename = '';
         % Masks
for i=1:Nexp
%     Masks{i}.ThresholdMask = ThMask{i};
    Masks{i}.RoiMask = RoiMask{i};
    if ~isempty(poly)
       Masks{i}.RoiCoordinates=poly{i};
    end
    Masks{i}.NumberRoi = nROI{i};
    Masks{i}.RoiNames = ROInames{i};
end

    % Raw data (without any threshold mask or correction)      
for i=1:Nexp
    Iref{i} = RawData{1,i}.Iref;   
end
     
     save(filename, 'Wspec','Zspec','AsymSpec','B0','Masks','Iref');
      %% DISPLAYING NEW ROI DATA
                                       
    Nexp = length(Masks);
     
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
                                       
    
if Nexp<=3
   nsub2=Nexp;
   nsub1=1;
   
elseif  Nexp>3 && Nexp<=6
    nsub2=3;
    nsub1=2;
    
elseif  Nexp>6 && Nexp<=9
    nsub2=3;
    nsub1=3;
else
    nsub1=0;
end

if strcmp(disSpectra,'OnePerExp_AllRoi')

   
       for i=1:Nexp

           nROI=Masks{i}.NumberRoi;
           leg=strcat('ROI',num2str((1:nROI)'));

           figure
           colours=colormap(lines(nROI));
           set(gcf,'Position',get(0,'ScreenSize'),'Color',[1,1,1])

           % Z spectra 
            subplot (1,2,1)
            hold on

            for b=1:nROI
                linecolour=colours(b,:);
                h1a=plot(Zspec{i}.IntFreq, Zspec{i}.IntValues(:,b));
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);          
            end  

            for b=1:nROI
                linecolour=colours(b,:);
                h1b=plot(Zspec{i}.OriginalFreq, Zspec{i}.OriginalValues(:,b));
                set(h1b,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
                        'LineWidth',2)        
            end 

            title(sprintf('Z Spectra ( Experiment %d, all ROI)',i),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,leg);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            grid on
            xlim([min(Zspec{i}.IntFreq),max(Zspec{i}.IntFreq)])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('Normalised Signal Intensity (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')


            % MTRasym spectra
            subplot (1,2,2)
            hold on

            for b=1:nROI
                linecolour=colours(b,:);  
                 h1a = plot(AsymSpec{i}.IntFreq, AsymSpec{i}.IntValues(:,b));
                 set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);
            end

            for b=1:nROI
                 linecolour=colours(b,:);        
                 h1b = plot(AsymSpec{i}.OriginalFreq, AsymSpec{i}.OriginalValues(:,b));
                 set(h1b,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
                        'LineWidth',2)
            end

            grid on       
            title(sprintf('MTR_A_S_S_Y_M Spectra( Experiment %d ), all ROI',i),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,leg);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            xlim([min(AsymSpec{i}.IntFreq),max(AsymSpec{i}.IntFreq)])
            ylim([-20,30])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('CEST Effect (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')
       end
   
   
elseif strcmp(disSpectra,'OnePerExp_OnePerRoi')
   
    
       for i=1:Nexp

           nROI=Masks{i}.NumberRoi;

           for b=1:nROI

               figure
               colours=colormap(lines(nROI));
               linecolour=colours(b,:);
               set(gcf,'Position',get(0,'ScreenSize'),'Color',[1,1,1])

                % Z spectra
                subplot(1,2,1)
                h1a=plot(Zspec{i}.IntFreq, Zspec{i}.IntValues(:,b));
                hold on
                h1b=plot(Zspec{i}.OriginalFreq, Zspec{i}.OriginalValues(:,b));
                set(h1b,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
                        'LineWidth',2)
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);          
                title(sprintf('Z Spectra ( ROI%d Experiment%d  )',b,i),'FontSize',16,'FontWeight','bold')
                grid on
                xlim([min(Zspec{i}.IntFreq),max(Zspec{i}.IntFreq)])
                xlabel('Saturation Offset (ppm)','FontSize',16)
                ylabel('Normalised Signal Intensity (%)','FontSize',16)
                set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')


                % MTRasym spectra
                subplot (1,2,2)
                h1a = plot(AsymSpec{i}.IntFreq, AsymSpec{i}.IntValues(:,b));
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);
                hold on
                h1b = plot(AsymSpec{i}.OriginalFreq, AsymSpec{i}.OriginalValues(:,b));
                set(h1b,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
                            'LineWidth',2)
                grid on       
                title(sprintf('MTR_A_S_S_Y_M Spectra ( ROI%d Experiment%d )',b, i),'FontSize',16,'FontWeight','bold')
                xlim([min(AsymSpec{i}.IntFreq),max(AsymSpec{i}.IntFreq)])
                ylim([-20,30])
                xlabel('Saturation Offset (ppm)','FontSize',16)
                ylabel('CEST Effect (%)','FontSize',16)
                set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')
           end
       end
    
    
elseif strcmp(disSpectra,'AllExp_OnePerRoi')

   for i=1:Nexp   
       nROI(i)=Masks{i}.NumberRoi;
   end
   
   if length(unique(nROI))>1
      error('To use "AllExp_OnePerRoi" you must have the same number of ROI per each experiment')   
   end
   
   leg=strcat('Experiment',num2str((1:Nexp)'));
   
   
       for b=1:nROI(1)
           figure
           colours=colormap(lines(Nexp));
           set(gcf,'Position',get(0,'ScreenSize'),'Color',[1,1,1])

           % Z spectra 
            subplot (1,2,1)
            hold on

            for i=1:Nexp
                linecolour=colours(i,:);
                h1a=plot(Zspec{i}.IntFreq, Zspec{i}.IntValues(:,b));
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);   
            end  

            for i=1:Nexp
                linecolour=colours(i,:);
                h1b=plot(Zspec{i}.OriginalFreq, Zspec{i}.OriginalValues(:,b));
                set(h1b,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
                        'LineWidth',2)        
           end 

            title(sprintf('Z Spectra for ROI %d ',b),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,leg);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            grid on
            xlim([min(Zspec{i}.IntFreq),max(Zspec{i}.IntFreq)])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('Normalised Signal Intensity (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')


            % MTRasym spectra
            subplot (1,2,2)
            hold on

            for i=1:Nexp
                linecolour=colours(i,:);  
                h1a = plot(AsymSpec{i}.IntFreq, AsymSpec{i}.IntValues(:,b));
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);
            end

            for i=1:Nexp
                 linecolour=colours(i,:);        
                 h1b = plot(AsymSpec{i}.OriginalFreq, AsymSpec{i}.OriginalValues(:,b));
                 set(h1b,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
                        'LineWidth',2)
            end

            grid on       
            title(sprintf('MTR_A_S_S_Y_M Spectra for ROI %d ',b),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,leg);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            xlim([min(AsymSpec{i}.IntFreq),max(AsymSpec{i}.IntFreq)])
            ylim([-20,30])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('CEST Effect (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')
       end
   
   
elseif strcmp(disSpectra,'No')  
    
else
    error('disSpectra not understood')

end

%% Display Subtracting spectra: Subtracting right and left ROIs for every
%% experiment
Nexp = length(Masks);
for b=1:Nexp
    
    for i=1:length(AsymSpec{1,b}.IntFreq)
        Subtract_Zspec_int(i,b) = (AsymSpec{1,b}.IntValues(i,1) - AsymSpec{1,b}.IntValues(i,2));
        
    end
    for i=1:length(AsymSpec{1,b}.OriginalFreq)
        Subtract_Zspec_original(i,b) = (AsymSpec{1,b}.OriginalValues(i,1) - AsymSpec{1,b}.OriginalValues(i,2));
        
    end    
    
end

colours=colormap(lines(Nexp));
figure;
hold on
            for i=1:Nexp
                linecolour=colours(i,:);  
                h1a_subs = plot(AsymSpec{1,1}.IntFreq(1,:),Subtract_Zspec_int(:,i));
                set(h1a_subs,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);
            end

%             for i=1:Nexp
%                  linecolour=colours(i,:);        
%                  h1b_subs = plot(AsymSpec{1,2}.OriginalFreq(1,:),Subtract_Zspec_original(:,i));
%                  set(h1b_subs,'Marker','.','LineStyle','none','MarkerSize',15,'Color',linecolour,...
%                         'LineWidth',2)
%             end
leg=strcat('Experiment',num2str((1:Nexp)'));
            grid on       
            title(sprintf('Subtracted left and right ROI spectra '),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,leg);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            xlim([min(AsymSpec{i}.IntFreq),max(AsymSpec{i}.IntFreq)])
            ylim([-20,30])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('CEST Effect (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')

%% Display the thresholded image with the ROI
   
    
    
        if nsub1~=0
           figure
           set(gcf,'Position',get(0,'ScreenSize'))
           
            for i=1:Nexp         
                subplot(nsub1,nsub2,i)
                imagesc(Iref{1,i}(:,:).*Masks{i}.ThresholdMask(:,:));
                colormap Gray
                hold on
                plot(Masks{i}.RoiCoordinates.X(:,:),Masks{i}.RoiCoordinates.Y(:,:),'LineWidth',3);
                legend(Masks{i}.RoiNames)
                axis off
                axis image
                title(sprintf('ROIs ( Experiment %d )',i),'FontWeight','bold')
            end
%             mtit(sprintf('Slice %d',z),'FontWeight','bold','fontsize',14)
        else
        
            for i=1:Nexp
                    figure

                    imagesc(Iref{1,i}(:,:).*Masks{i}.ThresholdMask(:,:));
                    colormap Gray
                    hold on
                    plot(Masks{i}.RoiCoordinates.X(:,:),Masks{i}.RoiCoordinates.Y(:,:),'LineWidth',3);
                    legend(Masks{i}.RoiNames)
                    axis off
                    axis image
                    title(sprintf('ROIs ( Experiment %d )',i),'FontWeight','bold')
            end
        end
    


   
   

  
