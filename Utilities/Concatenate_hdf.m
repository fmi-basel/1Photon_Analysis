%% Small utility to concatenate hdf5 files
% Julian Hinz @Luthi lab 2018


% Set parameters
% How many sessions have hdf5 files that need to be concatenated - needed
Sessions = 1; 

% Directory of hdf5 files - Just to help not to click through everything 
% set it to a directory that is close to all sessions that you want to concatenate - optional
Directory_Files = [];

% Set it to a generic directory, all files will be saved there with a number code
% You ca and copy the files later - needed
Saving_Directory = ''; 

%%%%%%%
% Code starts from here

% Select which Files you want to concatenate
All_nam = cell(Sessions, 1);
for i = 1:Sessions
    All_nam{i} = uipickfiles('FilterSpec', [Directory_Files '*.hdf5']);
end

% Run through Sessions and load files, concatenate and save

for i = 1:Sessions
    Temp = All_nam{i};
    Temp_hdf = cell(size(Temp, 2), 1);
    
    for ii = 1:size(Temp, 2)
        Temp_hdf{ii} = h5read(Temp{ii}, '/images');
    end
    
    % Now concatenate and write it back
    C = cat(3, Temp_hdf{:});
    hdf5write([Saving_Directory 'HDFconc_' num2str(i)], '/images', C);
    clear C Temp_hdf
end

disp('Finished with concatenation !')

