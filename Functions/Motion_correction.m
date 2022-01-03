%% Workflow Motion Correction on hdf5 and tif files
% written by Julian Hinz @Luthi lab 2018

function [Movie_all, Max_all] = Motion_correction(All_nam, max_reps_motioncorr, ...
    Tolerance_motion, T_DS_factor, Resize_factor, BW_all, psf, Data_Format_in,...
    Chunk_size, Recording_Speed, Indv_Template)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Import Images - Resizing - Motion registration - Individual
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Check if the user specified paticular sessions and or animals to be 
    % processed again
    
   	% Correlation map / movie
    Max_all = cell(size(All_nam, 1), 1);
    Movie_all = cell(size(All_nam, 1), 1);
    
    for p = 1:size(All_nam, 1)
        
        % Identify Sessions of the same Animal
        subfolder_extracted = All_nam{p};
        
        % Correlation map
        Max = cell(size(subfolder_extracted, 2), 1);
        Movies = cell(size(subfolder_extracted, 2), 1);

        for pp = 1:size(subfolder_extracted, 2)
           
            % Time Loading and MC process
            StartTimer = tic;
            
            % Find all files in the folder, either tif or h5
            Current_folder = [subfolder_extracted{pp} '\'];
             
            % Get Folder Name
            [~,name_R,~] = fileparts(subfolder_extracted{pp});
            name = ['S' name_R];

            if ~exist([Current_folder 'processed_data'])
                mkdir([Current_folder 'processed_data'])
            else
                disp('File was already processed - Results will be overwritten !')
            end
 
            searchString = fullfile(Current_folder, ['*.' Data_Format_in]);
            Filenames = dir(searchString);
            
            % Define Region for motion correction
            BW_x = BW_all{p}{1, 1};
            BW_y = BW_all{p}{1, 2};
            
            % Here deal with preallocation of resources on disk to improve
            % RAM efficiency and also chunk hdf5 file and create mem mapped
            % file
            
            if strcmpi(Data_Format_in, 'hdf5') && size(Filenames, 1) == 1
                
                Info = h5info([Current_folder Filenames.name], '/images');
                
                chunks_orig_s = [1:Chunk_size:ceil(Info.Dataspace.Size(3)/Chunk_size)*Chunk_size];
                chunks_orig_e = [chunks_orig_s(2:end) - 1 Info.Dataspace.Size(3)];
                chunks_orig_s(2:end) = chunks_orig_s(2:end) ;
                
                chunks_MCs = ceil(chunks_orig_s/T_DS_factor);                 
                chunks_MCe = ceil(chunks_orig_e/T_DS_factor); 
                
                numFiles = numel(chunks_orig_s);
                
                data = imresize(h5read([Current_folder Filenames.name], '/images'), ...
                        [Info.Dataspace.Size(1)*Resize_factor, Info.Dataspace.Size(2)*Resize_factor], 'bicubic');
                Class_M = class(data);

                switch Class_M
                    case 'uint16'
                        Max_Value = 2^16 - 1;
                    case 'uint8'
                        Max_Value = 2^8 - 1;
                    otherwise
                         error('The data format is not uint, please check the format and run this script again !')
                end
                
                savefast([Current_folder 'processed_data\dataset'], 'data');
                Data_handle = matfile([Current_folder 'processed_data\dataset']);
                clear data
                
                MC = zeros(Info.Dataspace.Size(1)*Resize_factor, ...
                    Info.Dataspace.Size(2)*Resize_factor, numel(1:T_DS_factor:Info.Dataspace.Size(3)), Class_M);
                savefast([Current_folder 'processed_data\MC'], 'MC');
                clear MC
                MC_f = matfile([Current_folder 'processed_data\MC'], 'Writable', true);

            elseif strcmpi(Data_Format_in, 'tif') || strcmpi(Data_Format_in, 'tiff')
                
                [~, idx] = sort_nat({Filenames.date}); % Order the indices according to recording time
                
                num = size(Filenames, 1);
                A = cell(num, 1);
                t = Tiff([Current_folder Filenames(idx(1)).name]);
                Size = size(read(t));
                
                % This is taking care of the case that the miniscope
                % produces an error and does't produce a tif containing
                % only 1080x1080 frames, to solve this I catch errors and
                % replace the broken frame with the average of the 2 frames
                % immediatly adjacent to it
                try
                    parfor z = 1:num
                        A{z} = imresize(loadtiff([Current_folder Filenames(idx(z)).name]), ...
                            [Size(1)*Resize_factor, Size(2)*Resize_factor], 'bicubic');
                    end
                catch
                    for z = 1:num
                        Temp_tif = loadtiff([Current_folder Filenames(idx(z)).name]);
                        
                        try
                            A{z} = imresize(Temp_tif, [Size(1)*Resize_factor, Size(2)*Resize_factor], 'bicubic');
                        catch
                            if iscell(Temp_tif) && numel(Temp_tif) == 3
                                Piece1 = Temp_tif{1};
                                Piece2 = Temp_tif{3};                            
                                Temp_tif{2} = mean(cat(3, Piece1(:,:,end), Piece2(:,:,1)), 3, 'native');
                                A{z} = imresize(cat(3, Temp_tif{:}), [Size(1)*Resize_factor, Size(2)*Resize_factor], 'bicubic');
                            else
                               error('This case has not occured before, please contact julian.hinz@fmi.ch if you cant fix it.')        
                            end
                        end
                    end
                end
                                
                data = cat(3 , A{:});
                clear A
                Size = size(data);
                
                chunks_orig_s = [1:Chunk_size:ceil(Size(3)/Chunk_size)*Chunk_size];
                chunks_orig_e = [chunks_orig_s(2:end) - 1 Size(3)];
                chunks_orig_s(2:end) = chunks_orig_s(2:end) ;
                
                chunks_MCs = ceil(chunks_orig_s/T_DS_factor);                 
                chunks_MCe = ceil(chunks_orig_e/T_DS_factor); 
                
                numFiles = numel(chunks_orig_s);
                                
                % Check the class of the tif files to make the conversion
                % into gray scale accurate
                Class_M = class(data);

                switch Class_M
                    case 'uint16'
                        Max_Value = 2^16 - 1;
                    case 'uint8'
                        Max_Value = 2^8 - 1;
                    otherwise
                         error('The data format is not uint, please check the format and run this script again !')
                end
                
                savefast([Current_folder 'processed_data\dataset'], 'data');
                Data_handle = matfile([Current_folder 'processed_data\dataset']);
                clear data
                
                MC = zeros(Size(1), Size(2), numel(1:T_DS_factor:Size(3)), Class_M);
                savefast([Current_folder 'processed_data\MC'], 'MC');
                clear MC
                MC_f = matfile([Current_folder 'processed_data\MC'], 'Writable', true);
            elseif strcmpi(Data_Format_in, 'isxd')
                   
                func_name = 'isx_movie_get_frame_data_u16';
                
                if size(Filenames, 1) ~= 1
                   
                    [~, idx] = sort_nat({Filenames.name}); % Order the indices according to recording time

                    % First define general dataformat
                    movie = isx.Movie.read([Current_folder Filenames(idx(1)).name]);
                    Class_M = movie.data_type;

                    switch Class_M
                        case 'uint16'
                            Max_Value = 2^16 - 1;
                        otherwise
                             error('The data format is not uint, please check the format or contact julian.hinz@fmi.ch and run this script again !')
                    end
                    
                    Data_Collection = cell(size(Filenames, 1), 1);
                    
                    parfor z = 1:size(Filenames, 1)
                        movie = isx.Movie.read([Current_folder Filenames(idx(z)).name]);
                        
                        data_tmp = cell(movie.timing.num_samples, 1);
                        frame_data_ptr = libpointer([movie.data_type, 'Ptr'], zeros(fliplr(movie.spacing.num_pixels), movie.data_type));

                        for zz = 1:movie.timing.num_samples
                            try
                                isx.internal.call_c_lib(func_name, movie.ptr, zz, frame_data_ptr);
                                data_tmp{zz} = imresize(frame_data_ptr.Value', [movie.spacing.num_pixels(1)*Resize_factor, ...
                                    movie.spacing.num_pixels(2)*Resize_factor], 'bicubic');
                            catch
                                disp('Movie dropped a frame !')
                            end
                        end
                        
                        data_tmp(cellfun(@isempty, data_tmp)) = [];
                        Data_Collection{z} = cat(3, data_tmp{:});
                    end
                    
                    data = cat(3, Data_Collection{:});
                    clear Data_Collection data_tmp
                else
                    movie = isx.Movie.read([Current_folder Filenames.name]);
                    Class_M = movie.data_type;

                    switch Class_M
                        case 'uint16'
                            Max_Value = 2^16 - 1;
                        otherwise
                             error('The data format is not uint, please check the format or contact julian.hinz@fmi.ch and run this script again !')
                    end

                    frame_data_ptr = libpointer([movie.data_type, 'Ptr'], zeros(fliplr(movie.spacing.num_pixels), movie.data_type));
                    
                    Num_Reps = ceil(movie.timing.num_samples/Chunk_size);
                    data_tmp = cell(Num_Reps, 1);
                    
                    for z = 1:Num_Reps
                        
                        start = 1+Chunk_size*(z-1);
                        
                        if z == Num_Reps
                           stop = movie.timing.num_samples;
                        else
                           stop = Chunk_size*z;
                        end
                        
                        data_tmpTmp = cell(numel(start:stop), 1);
                        
                        for zz = start:stop
                            try
                                isx.internal.call_c_lib(func_name, movie.ptr, zz, frame_data_ptr);
                                data_tmpTmp{zz} = imresize(frame_data_ptr.Value', [movie.spacing.num_pixels(1)*Resize_factor, ...
                                    movie.spacing.num_pixels(2)*Resize_factor], 'bicubic');
                            catch
                                disp('Movie dropped a frame !')
                            end
                        end
                        
                        data_tmp(cellfun(@isempty, data_tmp)) = [];
                        data_tmp{z} = cat(3, data_tmpTmp{:});
                    end
                    data = cat(3 , data_tmp{:});
                end

                Size_dat = size(data);
                savefast([Current_folder 'processed_data\dataset'], 'data');
                clear data
                
                MC = zeros(Size_dat(1), Size_dat(2), numel(1:T_DS_factor:Size_dat(3)), Class_M);
                savefast([Current_folder 'processed_data\MC'], 'MC');
                clear MC
                
                MC_f = matfile([Current_folder 'processed_data\MC'], 'Writable', true);
                Data_handle = matfile([Current_folder 'processed_data\dataset']);
                
                
                chunks_orig_s = [1:Chunk_size:ceil(Size_dat(3)/Chunk_size)*Chunk_size];
                chunks_orig_e = [chunks_orig_s(2:end) - 1 Size_dat(3)];
                chunks_orig_s(2:end) = chunks_orig_s(2:end) ;
                
                chunks_MCs = ceil(chunks_orig_s/T_DS_factor);                 
                chunks_MCe = ceil(chunks_orig_e/T_DS_factor); 
                
                numFiles = numel(chunks_orig_s);
                
            elseif strcmpi(Data_Format_in, 'mat')
                Data_handle = matfile([Current_folder 'processed_data\dataset']);
                Size_dat = size(Data_handle.data);
                Class_M = class(Data_handle.data);
                
                MC = zeros(Size_dat(1), Size_dat(2), numel(1:T_DS_factor:Size_dat(3)), Class_M);
                savefast([Current_folder 'processed_data\MC'], 'MC');
                MC_f = matfile([Current_folder 'processed_data\MC'], 'Writable', true);
                clear MC
                
                chunks_orig_s = [1:Chunk_size:ceil(Size_dat(3)/Chunk_size)*Chunk_size];
                chunks_orig_e = [chunks_orig_s(2:end) - 1 Size_dat(3)];
                chunks_orig_s(2:end) = chunks_orig_s(2:end) ;
                
                chunks_MCs = ceil(chunks_orig_s/T_DS_factor);                 
                chunks_MCe = ceil(chunks_orig_e/T_DS_factor); 
                
                numFiles = numel(chunks_orig_s);
                
                 switch Class_M
                        case 'uint16'
                            Max_Value = 2^16 - 1;
                        otherwise
                             error('The data format is not uint, please check the format or contact julian.hinz@fmi.ch and run this script again !')
                 end
            else
                error('So far only hdf5, tif, tiff, isxd and mat files are supported.')
            end
           
            % Save applied shifts
            Shift_collection = cell(numFiles, 1);

            %%%%
            % Main loop for Motion correction with first part of image
            for ppp = 1:numFiles
                               
                % load saved mat file
                A = Data_handle.data(:,:,chunks_orig_s(ppp):chunks_orig_e(ppp));

                % estimate Template only during first iteration of that
                % paticular Session and animal or all Sessions, if
                % individual template was selected
                
                if Indv_Template && ppp == 1 || ppp == 1 && pp == 1
                    Resized_Image = mat2gray(A, [0 Max_Value]);
                    clear A
                    
                    Shifts = 1; Start_time = tic; count = 0;
                    Registered_Images_template = Resized_Image(:,:,1:100);
                    while mean(abs(Shifts(:))) > Tolerance_motion && count < max_reps_motioncorr                 
                        count = count + 1;
                        [Registered_Images_template, Shifts] = Image_Registration_gray...
                            (Registered_Images_template, Registered_Images_template(:,:,5), 2, psf, BW_x, BW_y);                
                    end
                    Time_final = toc(Start_time);
                    fprintf('The Template was created within %d s\n', Time_final);            

                    Template = median(Registered_Images_template, 3); % create new template
                    clear Registered_Images_template

                    if T_DS_factor > 1
                        Resized_Image = Resized_Image(:,:,1:T_DS_factor:end);
                    end
                else
                    % Spatial and Temporal DS
                    if T_DS_factor > 1
                        Resized_Image = A(:,:,1:T_DS_factor:end);
                        clear A
                    else
                        Resized_Image = A;
                        clear A
                    end

                    % convert data type
                    Resized_Image = mat2gray(Resized_Image, [0 Max_Value]);
                end

                %%% Motion correction 
                MRDFT = Resized_Image;        
                Start_time = tic; count = 0; Shifts = 1;
                Shifts_all = nan(size(MRDFT, 3), 2, max_reps_motioncorr);

                % Run either as long as it takes to get to the motion tolerance
                % or stop after X iterations
                while mean(abs(Shifts(:))) > Tolerance_motion && count < max_reps_motioncorr                 
                    count = count + 1;
                    [MRDFT, Shifts] = Image_Registration_gray(MRDFT, Template, 2, psf, BW_x, BW_y);
                    Shifts_all(:,:,count) = Shifts;
                end
                Time_final = toc(Start_time);
                
                if count < max_reps_motioncorr  
                    fprintf(['Motion correction part %d out of %d converged after' ...
                        ' %d repetitions and took %e s\n'], ppp, numFiles, count, Time_final);
                else
                    disp(['Motion Correction didnt converge, within allocated repetitions,' ...
                        ' try to increase the size of the MC ROI withing the lens.'])
                end
                
                Shifts_all = nansum(Shifts_all, 3);
                Shift_collection{ppp} = Shifts_all;
                
                MC_f.MC(:,:,chunks_MCs(ppp):chunks_MCe(ppp)) = im2uint16(MRDFT);
            end

            % Generate the Quality Control metrics
            Indices = 1:20:size(MC_f.MC, 3);
            Sub_dat = mat2gray(MC_f.MC(:,:,Indices), [0 Max_Value]);

            Max{pp} = max(imfilter(Sub_dat, psf, 'replicate'), [], 3);

            F(numel(Indices)) = struct('cdata',[],'colormap',[]);
            
            g = figure;
            for j = 1:numel(Indices)
                imagesc(Sub_dat(:,:,j)); colormap bone; axis off
                title(['Time: ' num2str(j/Recording_Speed * 20) 's'])
                drawnow
                F(j) = getframe;
            end
            
            close(g);
            Movies{pp} = F;

            savefast([Current_folder 'processed_data\MC_Shifts.mat'], 'Shift_collection');
            
            if Indv_Template || pp == 1
                Template = im2uint16(Template);
                savefast([Current_folder 'processed_data\Template.mat'], 'Template');
            end
            
            clear Full_Session_MRDFT MRDFT Resized_Image F
            delete(gcp);    
            
            TTime = toc(StartTimer);
            
            disp('******************************************************')
            fprintf(['Finished with folder ' name '\n'])
            fprintf(['Motion correction and loading took ' num2str(TTime) ' s \n'])
            disp('******************************************************')
        end
        
        % Save Projection Image and DS Movie for Visualization
        Max_all{p} = Max;
        Movie_all{p} = Movies;
    end
end