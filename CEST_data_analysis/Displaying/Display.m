function Display(disWSpectra,disROI,disB0,disSFmap,disIntMap,disIntWASSRMap,disSpectra,B0,Iref,Masks,...
                  Zspec,Wspec,AsymSpec,SFmap,IntMap,IntWASSRMap)

%% ***** Displaying Function ***** %%

%INPUT:
%     disROI = 'Yes' if you wanna display the reference images with the ROI,
%              'No' otherwise. 
%
%     disB0 = Structure with the following fields
%                - disB0.is = 'Yes' or 'No' if you wanna display the B0maps or not
%                - disB0.range = [a,b] if you wanne normalize the values in B0map 
%                                to the range [a,b];
%                                'auto' for a display by default.
%
%     disB1 = Structure with the following fields
%                - disB1.is = 'Yes' or 'no' if you wanna display the B1maps or not
%                - disB1.range = [a,b] if you wanne normalize the values in B1map 
%                                to the range [a,b];
%                                'auto' for a display by default.
%
%     disSFmap = Structure with the following fields
%                - disSFmap.is = 'Yes' or 'no' if you wanna display the 
%                                 SF CEST maps or not
%                - disSFmap.range = [a,b] if you wanne normalize the values in the SF CEST maps
%                                   to the range [a,b];
%                                   'auto' for a display by default
%
%     disIntMap = Structure with the following fields
%                - disIntMap.is = 'Yes' or 'no' if you wanna display the Integral 
%                                  CEST maps or not
%                - disIntMap.range = same thing as disSFmap.range but for 
%                                    the Integral maps
%
%     disSpectra = several possibilities to display the Z and MTR spectra:
%                  'OnePerExp_OnePerRoi' if you wanna display one figure per
%                   each experiment and per each ROI.
%                  'OnePerExp_AllRoi' if you wanna display one figure per 
%                   each experiment and all ROI on the same figure.
%                  'AllExp_OnePerRoi' if you wanna display one figure per 
%                   each ROI and all experiments on the same figure.
%                  'No' if you don't wanna display the spectra
%
%     B0 = nExp cell array of structures, each cell i corresponds to one experiment
%          and contains the following fields:
%             - B0{i}.Map = B0 maps for experiment i (d1xd2xd3 array)
%             - B0{i}.TypeCor = Type of B0 correction for experiment i (string)
%
%     B1 = nExp cell array of structures, each cell i corresponds to one experiment
%          and contains the following fields:
%             - B1{i}.Map = B1 maps for experiment i (d1xd2xd3 array)
%             - B1{i}.TypeCor = Type of B1 correction for experiment i (string) 
%
%    Iref = nExp cell array of d1xd2xd3 arrays, each cell i contains the reference 
%           images for the experiment i.
%
%     Masks = nExp cell array of structures, each cell i corresponds to one experiment
%             and contains the following fields:
%             - Masks{i}.ThresholdMask = ThresholdMask for experiment i
%                                       (d1xd2xd3 array)
%             - Masks{i}.RoiMask = ROI masks for experiment i  
%                                 (d1xd2xd3xnROI array)
%             - Masks{i}.RoiCoordinates = Coordinates of each ROI for experiment i
%                                        (same as the output poly of the function Draw_Roi.m)
%             - Masks{i}.NumberRoi = Number of ROI for experiment i
%             - Masks{i}.RoiNames = Names of the ROI for experiment i
%                                   (nROI character array)
%
%     Zspec = nExp cell array of structures, each cell i corresponds to one experiment
%             and contains the following fields:
%              -Zspec{i}.IntValues = Interpolated Zspectra values for each
%                                    ROI 
%                                    nfreqNew x nROI 
%                                    nfreqNew > nfreq if the data are interpolated
%              -Zspec{i}.IntFreq = offset frequencies (in ppm) by increasing order 
%                                  after interpolation of the CEST data (nfreqNew vector)
%              -Zspec{i}.OriginalValues = Zspectra values for each ROI
%                                         at the original frequencies
%                                         nfreqOr x nROI
%              -Zspec{i}.OriginalFreq = original offset frequencies (in ppm) 
%                                       by increasing order (nfreqOr vector)
%
%     AsymSpec = same thing as Zspec but with the MTRasym
%                 spectra values.
%
%
%      SFmap = Structure containing the following fields :
%                    - SFmap.Map = Single Frequency CEST maps 
%                                  (d1 x d2  array); 
%                    - SFmap.Freq = same thing than SFfreq (input) 
%
%      IntMap = Structure containing the following fields :
%                    - IntMap.Map = Integral CEST maps
%                                   (d1 x d2  array); 
%                    - IntMap.FreqRange = same thing than rangeInt (input)
        
%%  Initial Parameters
Nexp = length(Masks);

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


%% Display the thresholded image with the ROI

if strcmp(disROI,'Yes')
    
    
    
        if nsub1~=0
           figure
           set(gcf,'Position',get(0,'ScreenSize'))
           
            for i=1:Nexp         
                subplot(nsub1,nsub2,i)
                imagesc(Iref{i}(:,:).*Masks{i}.ThresholdMask(:,:));
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

                    imagesc(Iref{i}(:,:).*Masks{i}.ThresholdMask(:,:));
                    colormap Gray
                    hold on
                    plot(Masks{i}.RoiCoordinates.X(:,:),Masks{i}.RoiCoordinates.Y(:,:),'LineWidth',3);
                    legend(Masks{i}.RoiNames)
                    axis off
                    axis image
                    title(sprintf('ROIs ( Experiment %d )',i),'FontWeight','bold')
            end
        end
    
    
end

%% Displaying the b0map

if strcmp(disB0.is,'Yes')
    
    
        if nsub1~=0
           figure
           set(gcf,'Position',get(0,'ScreenSize'))
           
            for i=1:Nexp
                subplot(nsub1,nsub2,i)
                if length(disB0.range)==2
%                    B0{i}.Map=cell2mat(B0{i}.Map); %B0 map
                   imagesc(B0{i}.Map(:,:),disB0.range);
                else
%                     B0{i}.Map=cell2mat(B0{i}.Map); %B0 map
                    imagesc(B0{i}.Map(:,:));
                end
                colorbar
                axis off
                axis image
                title(sprintf('B0 map ( Experiment %d )',i),'FontWeight','bold')
                
            end
%             mtit(sprintf('Slice %d',z),'FontWeight','bold','fontsize',14)
        else

            for i=1:Nexp
                figure
                if length(disB0.range)==2
                   imagesc(B0{i}.Map(:,:),disB0.range);
                else
                    imagesc(B0{i}.Map(:,:));
                end
                colorbar
                axis off
                axis image
                title(sprintf('B0 map ( Experiment %d ) ',i),'FontWeight','bold')          
            end

        end
    
end


%% Displaying SFmap

if strcmp(disSFmap.is,'Yes')
    
   

        if nsub1~=0
           figure
           set(gcf,'Position',get(0,'ScreenSize'))
          
            for i=1:Nexp
                subplot(nsub1,nsub2,i)
                if length(disSFmap.range)==2
                   imagesc(SFmap{i}.Map(:,:),disSFmap.range);
                else
                    imagesc(SFmap{i}.Map(:,:));
                end
                colorbar
                axis off
                axis image
                cest_freq=SFmap{i}.Freq;
                title( sprintf('Single Frequency CEST map ( %.1f ppm )\nExperiment %d',cest_freq,i),...
                       'FontWeight','bold')
            end
%             mtit(sprintf('Slice %d',z),'FontWeight','bold','fontsize',14)
            
        else

            for i=1:Nexp
                figure
                if length(disSFmap.range)==2
                   imagesc(SFmap{i}.Map(:,:),disSFmap.range);
                else
                    imagesc(SFmap{i}.Map(:,:));
                end
                colorbar
                axis off
                axis image
                cest_freq=SFmap{i}.Freq;
                title(  sprintf('Single Frequency CEST map ( %.1f ppm )\nExperiment %d ',cest_freq,i),...
                       'FontWeight','bold')              
            end

        end

    

end


%% Displaying Integral Maps

if strcmp(disIntMap.is,'Yes')

        
    

        if nsub1~=0
           figure
           set(gcf,'Position',get(0,'ScreenSize'))
           
           for i=1:Nexp
                subplot(nsub1,nsub2,i)
                if length(disIntMap.range)==2
                   imagesc(IntMap{i}.Map(:,:),disIntMap.range);
                else
                    imagesc(IntMap{i}.Map(:,:));
                end
                colorbar
                axis off
                axis image
                range=IntMap{i}.FreqRange;
                title( sprintf('Integral CEST map ( range = [ %.1f ppm %.1f ppm ] )\nExperiment %d',range(1),range(2),i),...
                       'FontWeight','bold')
           end
%            mtit(sprintf('Slice %d',z),'FontWeight','bold','fontsize',14)

        else

            for i=1:Nexp
                figure
                if length(disIntMap.range)==2
                   imagesc(IntMap{i}.Map(:,:),disIntMap.range);
                else
                    imagesc(IntMap{i}.Map(:,:));
                end
                colorbar
                axis off
                axis image
                range=IntMap{i}.FreqRange;
                title( sprintf('Integral CEST map ( range = [ %.1f ppm %.1f ppm ] )\nExperiment %d ',range(1),range(2),i),...
                       'FontWeight','bold')
            end
        end
    
end

%% Displaying WASSR Integral Maps

% if strcmp(disIntWASSRMap.is,'Yes')
% 
%         
%     
% 
%         if nsub1~=0
%            figure
%            set(gcf,'Position',get(0,'ScreenSize'))
%            
%            for i=1:Nexp
%                 subplot(nsub1,nsub2,i)
%                 if length(disIntWASSRMap.range)==2
%                    imagesc(IntWASSRMap{i}.Map(:,:),disIntWASSRMap.range);
%                 else
%                     imagesc(IntWASSRMap{i}.Map(:,:));
%                 end
%                 colorbar
%                 axis off
%                 axis image
%                 range=IntWASSRMap{i}.FreqRange;
%                 title( sprintf('New metric map ( range = [ %.1f ppm %.1f ppm ] )\nExperiment %d',range(1),range(2),i),...
%                        'FontWeight','bold')
%            end
% %            mtit(sprintf('Slice %d',z),'FontWeight','bold','fontsize',14)
% 
%         else
% 
%             for i=1:Nexp
%                 figure
%                 if length(disIntWASSRMap.range)==2
%                    imagesc(IntWASSRMap{i}.Map(:,:),disIntWASSRMap.range);
%                 else
%                     imagesc(IntWASSRMap{i}.Map(:,:));
%                 end
%                 colorbar
%                 axis off
%                 axis image
%                 range=IntWASSRMap{i}.FreqRange;
%                 title( sprintf('New metric map ( range = [ %.1f ppm %.1f ppm ] )\nExperiment %d ',range(1),range(2),i),...
%                        'FontWeight','bold')
%             end
%         end
%     
% end

%% Display Spectra

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

%% Display WASSR Spectra

if strcmp(disWSpectra,'OnePerExp_AllRoi')

   
       for i=1:Nexp

           nROI=Masks{i}.NumberRoi;

           figure
           colours=colormap(lines(nROI));
           set(gcf,'Position',get(0,'ScreenSize'),'Color',[1,1,1])

           % W Z spectra 
            subplot (1,2,1)
            hold on

            for b=1:nROI
                linecolour=colours(b,:);
                h1a=plot(Wspec{i}.IntFreq, Wspec{i}.IntValues(:,b));
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);          
            end  

 
            title(sprintf('WASSR Z Spectra ( Experiment %d)',i),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,Masks{i}.RoiNames);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            grid on
            xlim([min(Wspec{i}.IntFreq),max(Wspec{i}.IntFreq)])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('Normalised Signal Intensity (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')


           
       end
   
   
elseif strcmp(disWSpectra,'OnePerExp_OnePerRoi')
   
    
       for i=1:Nexp

           nROI=Masks{i}.NumberRoi;

           for b=1:nROI

               figure
               colours=colormap(lines(nROI));
               linecolour=colours(b,:);
               set(gcf,'Position',get(0,'ScreenSize'),'Color',[1,1,1])

                % Z spectra
                subplot(1,2,1)
                h1a=plot(Wspec{i}.IntFreq, Wspec{i}.IntValues(:,b));
    
                
               
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);          
                title(sprintf('WASSR Z Spectra ( ROI%d Experiment%d )',b,i),'FontSize',16,'FontWeight','bold')
                grid on
                xlim([min(Wspec{i}.IntFreq),max(Wspec{i}.IntFreq)])
                xlabel('Saturation Offset (ppm)','FontSize',16)
                ylabel('Normalised Signal Intensity (%)','FontSize',16)
                set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')


               
           end
       end
    
    
elseif strcmp(disWSpectra,'AllExp_OnePerRoi')

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

           % W Z spectra 
            subplot (1,2,1)
            hold on

            for i=1:Nexp
                linecolour=colours(i,:);
                h1a=plot(Wspec{i}.IntFreq, Wspec{i}.IntValues(:,b));
                set(h1a,'MarkerSize',4,'MarkerFaceColor',linecolour,'Color',linecolour,...
                        'LineWidth',2);   
            end  

           
            title(sprintf('WASSR Z Spectra for ROI %d ',b),'FontSize',16,'FontWeight','bold')
            legend1=legend(gca,leg);
            set(legend1,'Location','EastOutside','FontSize',14,'Box','off');
            grid on
            xlim([min(Wspec{i}.IntFreq),max(Wspec{i}.IntFreq)])
            xlabel('Saturation Offset (ppm)','FontSize',16)
            ylabel('Normalised Signal Intensity (%)','FontSize',16)
            set(gca,'Xdir', 'reverse','FontSize',14,'LineWidth',2,'FontWeight','demi','XMinorTick','on','YMinorTick','on')


           
       end
   
   
elseif strcmp(disWSpectra,'No')  
    
else
    error('disWSpectra not understood')

end