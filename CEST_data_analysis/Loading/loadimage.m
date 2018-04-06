function [I_ref]=loadimage(images)


% User select image
disp('User input required: Select image files.')
start_path='/media/OS/CEST_Varian';

for a=1:images
    [file,path] = uigetfile({'*.fdf;*.FDF;*.nii;*.NII;*.nii.gz;*.hdr;*.HDR',...
        'Varian FDF-files, NIfTI and Analyze 7.5 files (*.fdf,*.nii,*.nii.gz,*.hdr)';...
        '*.*', 'All Files (*.*)'},...
        ['Select NIfTI, Analyze or Varian FDF file for reference scan ',num2str(a)],start_path);
    if file == 0
        error('loadCESTvarian:LoadRefScans','User cancelled file selection.')
    end
    
    disp(['Reference scan ',num2str(a),': ',path,file,'...'])
    
    reffiles{a,1}=fullfile(path,file);
    
    if ispc
        ind = find(reffiles{a,1} == '\', 1, 'last');
    else
        ind = find(reffiles{a,1} == '/', 1, 'last');
    end
    start_path = reffiles{a,1}(1:ind);
end
% Load images scans

% -- Find file extension --
for a=1:images
    disp(['Loading image ',num2str(a),'...'])
    ind = find(reffiles{a,1} == '.', 1, 'last');
    if isempty(ind)
        error('loadimagevarian','Cannot find file extension!');
    else
        ext = reffiles{a,1}(ind:end);
    end
    
    % -- Determine read function --
    if strcmpi(ext,'.nii') || strcmpi(ext,'.gz') || strcmpi(ext,'.hdr')
        I_ref_tmp  = aedes_read_nifti(reffiles{a,1});
        I_ref{a,1} = I_ref_tmp.FTDATA;
    elseif strcmpi(ext,'.fdf')
        I_ref_tmp  = aedes_readfdf(reffiles{a,1});
        I_ref{a,1} = I_ref_tmp.FTDATA;
    end
end
end


