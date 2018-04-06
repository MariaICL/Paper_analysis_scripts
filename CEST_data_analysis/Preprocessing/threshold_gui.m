function [mask,threshold]=threshold_gui(input_image,Noise)

%% **** Threshold GUI: Interactive image thresholding - Ryan Bendell 2012 **** %%
% If the image noise has already been calculated, this tool enables the user to 
% choose a level of SNR and to see the corresponding thresholded image. 

%INPUT: 
%      input_image = image to be thresholded 
%                   (d1 x d2 x d3 array) with [d1,d2] the dimensions of each image
%                                              d3 the number of slices
%
%      Noise = value of the input_image noise (calculated by preprocessing.m)

%OUTPUT:
%      mask = Threshold mask (d1 x d2 x d3 array)
%
%      threshold = Value of the threshold (= SNR*Noise)
     

%% ------------------------------------------------------------------------
% Data Initialisation
% ------------------------------------------------------------------------

DATA.input_image  = input_image;
DATA.output_image = input_image;
DATA.mask         = zeros(size(input_image));
DATA.threshold    = 0;

% ------------------------------------------------------------------------
% GUI Construction
% ------------------------------------------------------------------------

% Main figure
h_fig = figure;           
set(h_fig,...
    'Menubar','none',...
    'Name','Threshold Masking GUI',...
    'NumberTitle','off',...
    'Toolbar','none',...
    'Units','normalized',...
    'Color',get(0,'defaultUicontrolBackgroundColor'),...
    'Resize','off',...
    'Colormap',gray);

% Instruction text
h_inst = uicontrol(...    
    'Parent',h_fig,...
    'Style','text',...
    'String',['Adjust the SNR using the bottom ',...
        'text box or slider. Adjust image scaling using the box in ',...
        'the top right.Press "Done" to save.'],...
    'Units','normalized',...
    'Position',[.01,.87,.77,.12]);

% Image axes
h_axes=axes(...        
    'Parent',h_fig',...
    'Position',[.1,.15,.8,.7]);

% Text edit for setting maximum of image scaling
h_scal = uicontrol(...    
    'Parent',h_fig,...
    'Style','edit',...
    'String',10000,...
    'TooltipString','Image Scale Maximum',...
    'BackgroundColor',[1 1 1],...
    'Units','normalized',...
    'Position',[0.79,0.92,0.13,0.05],...
    'Callback',@h_scal_callback);

% Text edit for setting slider
h_sled = uicontrol(...    
    'Parent',h_fig,...
    'Style','edit',...
    'String',0,...
    'TooltipString','Masking Threshold (Percentile)',...
    'BackgroundColor',[1 1 1],...
    'Units','normalized',...
    'Position',[0.07,0.05,0.12,0.05],...
    'Callback',@h_sled_callback);

% Threshold slider
h_slid = uicontrol(...    
    'Parent',h_fig,...
    'Style','slider',...
    'Max',100,...
    'Min',0,...
    'Value',DATA.threshold,...
    'SliderStep',[0.01,0.1],...
    'Units','normalized',...
    'Position',[.2,.05,.58,.05],...
    'Callback',@h_slid_callback);

% Done button
h_done = uicontrol(...     
    'Parent',h_fig,...
    'Style','pushbutton',...
    'String','Done',...
    'Units','normalized',...
    'Position',[.8,.05,.1,.05],...
    'Callback',@h_done_callback);


% ------------------------------------------------------------------------
% GUI Initialisation
% ------------------------------------------------------------------------

displayimage(DATA.input_image)          % Display input image in axes
set(h_fig,'Visible','on')               % Make figure visible
uiwait(h_fig)

% ------------------------------------------------------------------------
% Callback functions
% ------------------------------------------------------------------------

    function h_scal_callback(hObject,eventdata)
        displayimage(DATA.output_image)
    end

    function h_sled_callback(hObject,eventdata)
        th = str2num(get(hObject,'String'));
        set(h_slid,'Value',th)
        DATA.threshold = Noise*th;
        DATA.mask = DATA.input_image>=DATA.threshold;
        DATA.output_image = DATA.input_image.*DATA.mask;
        displayimage(DATA.output_image)
    end

    function h_slid_callback(hObject,eventdata)
        th = get(hObject,'Value');
        set(h_sled,'String',th)
        DATA.threshold = Noise*th;
        DATA.mask = DATA.input_image>=DATA.threshold;
        DATA.output_image = DATA.input_image.*DATA.mask;
        displayimage(DATA.output_image)
    end

    function h_done_callback(hObject,eventdata)
        mask = DATA.mask;
        threshold = DATA.threshold;
        close(h_fig)
    end

% ------------------------------------------------------------------------
% Utility functions
% ------------------------------------------------------------------------

      function displayimage(image)
        scalemax = str2num(get(h_scal,'String'));         
        try
           h_input = imagesc((image),[0,scalemax]);
        end
        
        try
           h_input = imagesc(image,[0,scalemax]);
        end
        
        axis image off
        colorbar
    end

end