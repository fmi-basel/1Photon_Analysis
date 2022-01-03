%% Analyze 1P Imaging Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you want to restart the analysis, please load the Restart_(n).mat file and
% The Current_Status variable indicates which processing step the code is
% finished with. If in doubt, execute Current_Status in the Command winndow
% To restart, navigate to the indicated section (i.e. Motion Correction 
% finished, means start from CNMFE) and press 'Run and advance' in the editor 
% Section of Matlab

clear
clc
close

%%%%% Set Basic Parameters and select Folder for processing
% Specifies the initiallization folder for selection of Data and where the restart 
% variable is saved - If empty, files will be saved in directory that
% Matlab is currently initallized
Direc_Imaging = ''; % Don't forget last backward Slash !

% How many animals should be processed ? For each animal you can select as
% many Sessions as you want later. Every Animal will give you a prompt that
% allows you to select the Folder directories of the sessions you want to
% process for this animal
Animals2process = 1;

% Select the data type you want to process, currently the code is able to
% process tif, tiff, isxd, hdf5 and mat files. For extension requsts,
% please write me an email (julian.hinz@fmi.ch)
% If you select isxd, you will need to provide a Folder location with a
% valid ISDP program installed. The code utlizes the infrastructure
% provided by insxopix to load files
Data_Format_in = 'isxd'; % Specify if files were saved as hdf5, tif, tiff or isxd, or mat 
ISDP_Path = 'C:\Program Files\Inscopix\Data Processing\'; % Change to your local installation 

Recording_Speed = 20; % Recording speed of miniscope (in HZ)
length_identifier = 6; % Length of animal identifier e.g. 865347 -> 6 // The code uses this to generate a final Folder Output for the different Animals
Folder_Structure = 'FT'; % If the tifs are in a structure - Animal/Animal_Session/Tif - FFT, if Animal_Session/Tif - FT 

%%%%% Motion Correction
T_DS_factor = 1; % Used for Data DS
Spatial_Downsampling = 0.25; % Scale either 0.5 (half size), or 0.25(Quarter Size)-recommended // Improves Speed and RAM usage
max_reps_motioncorr = 6; % Number of maximally allowed repititions of motion corrections
Tolerance_motion = 0.01; % Number of tolerated motion before motion correction is stopped measure is = mean(shift of Frames) < Tolerance 
Indv_Template = true; % Do you want to generate a new template for each Session // Recommendation - leave it at true
Chunk_size = 3000; % How many frames are processed at the same time, decrease for lower RAM usage
Qualitycheck = true; % Creates a DS video and max intensity projection Images for quality control
Size_mov = 500; % Size of the displayed movie for quality control

%%%%% CNMFE
% This feature enables processing of large files. Expect that CNMFE
% will require ~20 GB of RAM for a 60 min recording @20Hz/core
NW2Test = [12 10 8 4]; % Select how many workers are supposed to be assigned during CNMFE - depends on emperical evidence, if empty standard will be used

%%%% Email parameters
% will send you a message when Matlab finished with all MC and CNMFE
Enable_Email = false; % Do you want to receive Email notifications, true/false
Email_sender = ''; % Set to gmail mail adress - need to enable third party program access
Email_receiver = ''; % Your email adress
PWD = ''; % The PWD to the GMail

%%%% Post processing
% These parameters will provide some automatic exclusion of CNMFE
% components. With the current parameter set expect the generation of many
% false positives, which proved most desirable for our use case, but
% requires rigourous exclusion. These parameters take care of all neurons
% that are not at least x pixel large (minPixel), that touch the outer rim
% of the recording area (Outer_area) and might therefore include Motion
% Correction artifacts, and (the most important) exclude cells that are low
% signal to noise and are overlapping with other components a lot. This
% feature works against the oversampling by selecting the lowest S/N
% components (1-Threshold_Low_SN) and checks whether they are overlapping
% more than indicated in Overlap and exclude them if they do.

Threshold_Low_SN = 0.6; % how many neurons are consider high S/N - all others are checked for overlap
Overlap = 0.6; % This is the Overlap neurons need to have to get kicked out
Outer_area = 6; % in pix, excludes all shapes that lie in this area 
minPixel = 10; % minimum Size of components

%%%%%%%  Registration over days
% Choose your registration Method. From experience Baifra 'Ahanonu's
% registration works well on datasets with both many and few components,
% while the method of Liron 'Sheintuch' does not work well with datasets
% with only few neurons (for example in subcortical regions). We are now
% using only the 'Ahanonu' method of centroid distance and component
% spatial correlation
Registration_Methods = 'Ahanonu'; % Sheintuch, 'Ahanonu'
reference_session_index = 1; % Reference Session for registration

%%% Sheintuch
% Registration 1
initial_registration_type = 'Centroid distance'; % 'best_model_string' , 'Spatial correlation', 'Centroid distance'
Model_Selection = false; % If true the code will use calculated threshold
Threshold_manual = 5; % Only relevant if above parameter is set to true default: 5 for Centroid and 0.65 for correlation

% Final Registration
registration_approach = 'Simple threshold'; % either 'Probabilistic' or 'Simple threshold', if Simple threshold - Registration type2 and threshold manual2 are important, otherwise threshold p and p_same

Registration_type2 = 'Centroid distance'; %'Spatial correlation' or 'Centroid distance' - if automatic selection is wanted , leave blank []
Threshold_manual2 = 5; % Only relevant if above parameter is set to true default: 5 for Centroid and 0.65 for correlation

Threshold_p = 0.95; %
p_same_threshold = 0.5; % Increasing the threshold will reduce the number of aligned cells as it will increase the threshold for 'same'

%%% Ahanonu
max_Distance = 5; % in pix, distance that is maximally allowed for still belonging to the same group
analysis_Type = 'pairwise'; % 'clustering'
runMotionCorrection = true; % refers to whether the code will try to shift components to optimize alignment
imageCorr = 0.5; % [0 - 1] minimum spatial correlation necessary for aligning components across days
imageCorrType = 'corr2'; % 'corr2', 'ochiai', 'jaccard' - type of spatial correlation used for assesing alignment quality

%%%%%% Post-Post process
% These are just visualization parameters that you can use to assess
% whether components were correctly aligned. We use both the overlap of
% spatial components, as well as the similarity of Transient structure to
% assess whether a neiron was corretly aligned
D_aligned = []; % if empty the code will show all components that are aligned over at least 2 sessions, otherwise the number will indicate how many days are allowed to be unaligned 
Display = 'Raw'; % Either Raw, Deconv or both - decides between the different plot options

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Start of Main functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(gcp('nocreate'))
    delete(gcp)
end

% Limit Matlab's ability to use paralel cores for single core computations
last_comp_threads = maxNumCompThreads(4);

% Add folder directory automatically
filePath = matlab.desktop.editor.getActiveFilename;
One_up = strfind(filePath, '\');
addpath(genpath(filePath(1:One_up(end))));

Fs = Recording_Speed/T_DS_factor; % recording speed (after Downsampling)

if strcmpi(Data_Format_in, 'isxd') 
    try
        if ~exist(ISDP_Path)
            error('Please provide valid path to the ISDP Folder !')
        else
            addpath(genpath(ISDP_Path));
        end
    catch
        error('Please provide valid path to the ISDP Folder !')
    end
end

% Save variables for easy restarting
Current_Status = 'Before MC';

Saved = false;
Counter_mat = 1;

while ~Saved  
    if ~exist([Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat'])    
        save([Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat'], '-v7.3');
        disp(['If you want to Restart, please load following matfile ' ...
            [Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat']])
        Saved = true;
    else
        Counter_mat = Counter_mat + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Motion Correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://www.osapublishing.org/ol/abstract.cfm?uri=ol-33-2-156
% Code: https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation

% Interactive Selection
All_nam = cell(Animals2process, 1);
for i = 1:size(All_nam, 1)
    All_nam{i} = uipickfiles('FilterSpec', Direc_Imaging);
end

All_nam(cell2mat(cellfun(@isempty, All_nam, 'UniformOutput', false)) | ...
    cellfun(@(x) isequal(x, 0), All_nam)) = [];

% Select Regions for motion correction in all Sessions that are selected
BW_all = cell(numel(All_nam), 1);

for i = 1:numel(All_nam)
        
    % Find all .tif files in the folder
    Current_folder = [All_nam{i}{1} '\'];

    searchString = fullfile(Current_folder, ['*.' Data_Format_in]);
    Filenames = dir(searchString);
    [~, idx] = sort_nat({Filenames.name}); % Order the indices according to recording time

    if strcmpi(Data_Format_in, 'hdf5') && size(Filenames, 1) == 1
        Info = h5info([Current_folder Filenames(idx(1)).name], '/images');
        imageData = mat2gray(h5read([Current_folder Filenames(idx(1)).name], '/images', [1 1 1], ...
            [Info.Dataspace.Size(1) Info.Dataspace.Size(2) 1]));
    elseif strcmpi(Data_Format_in, 'tif') || strcmpi(Data_Format_in, 'tiff')     
        t = Tiff([Current_folder Filenames(idx(1)).name]);
        imageData = mat2gray(read(t)); 
    elseif strcmpi(Data_Format_in, 'isxd')
        movie = isx.Movie.read([Current_folder Filenames(idx(1)).name]);
        imageData = mat2gray(movie.get_frame_data(5));
    elseif strcmpi(Data_Format_in, 'mat')
        HandleDat = matfile([Current_folder 'processed_data\dataset.mat']);
        imageData = HandleDat.data(:,:,5);
    end
    
    % Now select regions for Motion Correction
    if strcmpi(Data_Format_in, 'mat')
        [Position_all] = Define_ROIS_mini(double(imageData)./imgaussfilt(...
            double(imageData), 15), 1);
    else
        Resized_Image = imresize(imageData, [size(imageData, 1)*Spatial_Downsampling ...
            size(imageData, 2)*Spatial_Downsampling], ...
                    'bicubic');
        [Position_all] = Define_ROIS_mini(Resized_Image./imgaussfilt(Resized_Image, 15), 1);
    end
    
    % Save the regions for later processing
    BW_all{i} = Position_all;
end

% Motion Correction and Alingment
% Create Filter used for Motion correction - This filter is taken from the
% CNMFE initiallization procedure and works well for MC

gSig = 4; % Determines the size of the used gaussian kernel
psf = fspecial('gaussian', ceil(gSig*4+1), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;

[Movie_all, CN_all] = Motion_correction(All_nam, max_reps_motioncorr, ...
    Tolerance_motion, T_DS_factor, Spatial_Downsampling, BW_all, psf, Data_Format_in,...
    Chunk_size, Recording_Speed, Indv_Template);

% Save variables for easy restarting
Current_Status = 'MC finished';
save([Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat'], '-v7.3');

% Send Email that Motion Correction is finished
if Enable_Email
    try
        send_mail_message(Email_receiver, Email_sender, ['Finished the Motion' ...
            'Correction - Input required !'], PWD, [], []);
    catch
        disp(['Email Couldnt be send, please double check parameters and whether' ...
            ' third party access is enabled in Gmail!'])
    end
end

% Option to check Motion Correction before proceeding
if Qualitycheck 
    ScreenSize = get(0,'ScreenSize');   
    for i = 1:size(All_nam, 1)
        Temp_Ses = All_nam{i};
        Handles = cell(size(Temp_Ses, 2), 1);
        Reps = floor((ScreenSize(3)- 100)/(Size_mov*1.1));
        Index_Positions_s = unique([1:Reps:size(Temp_Ses, 2) size(Temp_Ses, 2)]);
        Index_Positions_e = unique([1:Reps:size(Temp_Ses, 2) size(Temp_Ses, 2)]);
       
        for iii = 1:ceil(size(Temp_Ses, 2)/Reps)
            for ii = 1:Reps
                try
                    [~,namesP,~] = fileparts(Temp_Ses{ii+((iii-1)*Reps)});
                    figure('Position', [100 + (Size_mov*1.1*(ii-1)) 200 + ...
                        Size_mov Size_mov Size_mov])
                    imagesc(CN_all{i}{ii+((iii-1)*Reps)}); axis off; ...
                        colormap bone; box off; title(namesP)
                    Handles{ii} = implay(Movie_all{i}{ii+((iii-1)*Reps)});
                    Handles{ii}.Parent.Position = [100 + (Size_mov*1.1*(ii-1)) ...
                        100 Size_mov Size_mov];
                catch
                    continue
                end 
            end
            disp('Please press enter if you want to see the next batch of movies.')
            pause; % Pause execution until key stroke
            
            for ii = 1:Reps % close all windows
                try
                    close(Handles{ii})
                catch
                    continue
                end
            end
            close all                
        end
    end
    
    % Check if the Motion Correction was satisfying 
    answer = questdlg('Would you like to proceed with processing ?', ...
        'MC', ...
        'Yes', 'No', 'Yes');

    switch answer
        case 'No'
            error('You aborted the process, adjust the MC parameters and run the process again!')
    end 
    
end

% If Yes is selected, then create Aligment Folder 
[Savepath, ~ ,~] = fileparts(All_nam{1}{1});

switch Folder_Structure
    case 'FFT'
        Indices_ups = strfind(Savepath, '\');
        Alignment_path = [Savepath(1:Indices_ups(end)-1) '\Alignment\'];
    case 'FT'
        Alignment_path = [Savepath '\Alignment\'];
end

if ~exist(Alignment_path)
    mkdir(Alignment_path);
end

% Now save the MC read outs in the Alignment folder

Motion_Correction_QC{1} = Movie_all;
Motion_Correction_QC{2} = CN_all;

Saved = false;
Counter = 1;

while ~Saved     
    if ~exist([Alignment_path 'MC_Control_' num2str(Counter) '.mat'])    
        savefast([Alignment_path 'MC_Control_' num2str(Counter) '.mat'], 'Motion_Correction_QC');
        disp(['If you want to check the Motion Correction, please refer to ' [Alignment_path 'MC_Control_' num2str(Counter) '.mat']])
        Saved = true;
    else
        Counter = Counter + 1;
    end
end
clear Movie_all CN_all Motion_Correction_QC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CNMFE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://elifesciences.org/articles/28728
% Code: https://github.com/zhoupc/CNMF_E

% Run CNMFE for all selected files
for p = 1:size(All_nam, 1)
    CNMFE(All_nam{p}, Fs, NW2Test);
    vars = who('-regexp', 'mat_data');
    clear(vars{:})
end
    
close all

if Enable_Email
    try
    	send_mail_message(Email_receiver, Email_sender, ['Finished the CNMFE' ...
            ' - Input required !'], PWD, [], []);
    catch
        disp(['Email Couldnt be send, please double check parameters and whether' ...
            ' third party access is enabled in Gmail!'])
    end
end

% Save variables for easy restarting
Current_Status = 'CNMFE finished';
save([Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat'], '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post-process CNMFE results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adapts a small sub part of the CNMFE to call all components and resave
% them after post processing - deletion of components
% based on code found in: https://github.com/zhoupc/CNMF_E

post_process_CNMFE(All_nam, Threshold_Low_SN, Overlap, Outer_area, minPixel)

Current_Status = 'PostProcess finished';
save([Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat'], '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Register Cells over days 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% https://www.cell.com/cell-reports/abstract/S2211-1247(17)31430-4
% Code: https://github.com/zivlab/CellReg
% http://science.sciencemag.org/content/363/6424/276
% Code: https://github.com/bahanonu/ciatah

Imaging_Method = 'one_photon';

switch Registration_Methods
    case 'Sheintuch'
        [Output_Direc] = Sheintuch_registration_over_days(All_nam, Imaging_Method, ...
            length_identifier, Folder_Structure, Threshold_p, p_same_threshold, ...
            registration_approach, initial_registration_type, Threshold_manual, ...
            Registration_type2, Threshold_manual2, Model_Selection, ...
            reference_session_index);
    case 'Ahanonu'
        [Output_Direc] = Ahanonu_registration_over_days(All_nam, length_identifier, ...
            Folder_Structure, max_Distance, reference_session_index, analysis_Type, ...
            runMotionCorrection, imageCorr, imageCorrType);
end

% Save variables for easy restarting
Current_Status = 'Alignment finished';
save([Direc_Imaging 'Restart_' num2str(Counter_mat) '.mat'], '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check alignment of cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Alignment_post_Process(Output_Direc, D_aligned, Display, Fs, ...
    Registration_Methods);

