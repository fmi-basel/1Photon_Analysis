%% New registration, based on Schnitzer publication
% Wrapper written by Julian Hinz, 2019
% Original function for registration from http://science.sciencemag.org/content/363/6424/276.
% Taken from the Github repository:
% https://github.com/bahanonu/calciumImagingAnalysis 
% On 09.04.2019

function [Output_Direc] = Ahanonu_registration_over_days(All_nam, length_identifier, ...
    Folder_Structure, max_Distance, reference_Session, analysis_Type, ...
    runMotionCorrection, imageCorr, imageCorrType)

Output_Direc = cell(size(All_nam, 1), 1);

    % First Sort out which Animals/Sessions the user wants to process
    for i = 1:size(All_nam, 1)

       Sessions = All_nam{i};

        spatial_footprints = cell(1, size(Sessions, 2));
        Temporal_Traces =  cell(size(Sessions, 2), 3);
        [Savepath, Foldern ,~] = fileparts(All_nam{i}{1});
        Alignment = cell(5, 1);
        Corr_Img = cell(size(Sessions, 2), 1);

        for ii = 1:size(Sessions, 2)

            try
                load([All_nam{i}{ii} '\processed_data\CNMFE_Data.mat']); 
                load([All_nam{i}{ii} '\processed_data\CN.mat']);
                
            catch
                if exist([All_nam{i}{ii} '\processed_data\Result_CNMFE_noResults.mat']) ~= 0 
                    disp('No Neurons detected in this session - Registration of this session not executed.')
                    continue
                else
                    [~,ses, ~] = fileparts(All_nam{i}{ii});
                    disp(['Session ' ses ' could not be processed, check why and return to the Script.'])
                    continue
                end
            end
            
            try
                spatial_footprints{ii} =  reshape(full(Data.A), [size(Data.A_reshaped, 1) ...
                    size(Data.A_reshaped, 2) size(Data.A, 2)]);
            catch  
                Spatial_Components = full(Data.A);
                Weighted_Components = double(Spatial_Components./max(...
                    Spatial_Components) > 0.5);
                Data.A_reshaped = reshape(Weighted_Components, ...
                    [size(Cn, 1) size(Cn, 2), size(Weighted_Components, 2)]);         
                spatial_footprints{ii} =  reshape(full(Data.A), [size(Data.A_reshaped, 1) ...
                    size(Data.A_reshaped, 2) size(Data.A, 2)]);
            end
            
            Corr_Img{ii} = Cn;
            Temporal_Traces{ii, 1} = Data.C;
            Temporal_Traces{ii, 2} = Data.S;
            Temporal_Traces{ii, 3} = Data.C_raw;
        end

        clear Data Cn
        Non_registered_Ses = find(cellfun(@isempty, spatial_footprints));
        Temporal_Traces(Non_registered_Ses, :) = [];
        spatial_footprints(Non_registered_Ses) = [];

        % Find proper naming for the results directory  
        Find_pos = isstrprop(Foldern, 'digit');
        Index = strfind(Find_pos, ones(1, length_identifier));

        % Defining the results_directory and creating the figures_directory:

        switch Folder_Structure
            case 'FFT'
                Indices_ups = strfind(Savepath, '\');
                results_directory = [Savepath(1:Indices_ups(end)-1) '\Alignment\S_' Foldern(Index:Index+length_identifier - 1) '\'];
                Alignment_path = [Savepath(1:Indices_ups(end)-1) '\Alignment\'];
                Output_Direc{i} = results_directory;

            case 'FT'
                results_directory = [Savepath '\Alignment\S_' Foldern(Index...
                    :Index+length_identifier - 1) '\'];
                Alignment_path = [Savepath '\Alignment\'];
                Output_Direc{i} = results_directory;
        end

        if ~exist(results_directory)
            mkdir(results_directory);
        end


        figures_directory=fullfile(results_directory,'Figures');

        if exist(figures_directory,'dir')~=7
            mkdir(figures_directory);
        end

        figures_visibility = 'on'; % either 'on' or 'off' (in any case figures are saved)

        if size(All_nam{i}, 2) == 1 
            if exist([All_nam{i}{1} '\processed_data\Result_CNMFE_noResults.mat']) ~= 0
                continue
            else
                Data.S = Temporal_Traces{1, 2};
                Data.C = Temporal_Traces{1, 1};
                Data.A = spatial_footprints;
                Data.C_Raw = Temporal_Traces{1, 3}; 
                savefast([results_directory '\Data_Miniscope.mat'], 'Data')
                disp('Only one Session provided, no Alignment is performed !')
                continue
            end
        end
    
        % aligning the Sessions - powered by the CIAtah
        [OutStruct] = matchObjBtwnTrials(spatial_footprints, max_Distance, ...
            reference_Session, analysis_Type, runMotionCorrection, imageCorr, ...
            imageCorrType);
        
        close all
        
        %%% Post Processing 
        try       
            Num_Ses_aligned = cumsum(transpose(OutStruct.globalIDs~=0))';

            b = hist(Num_Ses_aligned(:,end), [1:size(Num_Ses_aligned, 2)]);
            g = figure;
            bar(b, 'k'); box off
            xlabel('Number of aligned Sessions'); ylabel('# Cells')
            saveas(gcf,[results_directory '\Alignment_Result.png'])
            close(g)

            [Ses, index] = sortrows(Num_Ses_aligned, 'descend');
            All_Ses = max(find(Ses(:,end) == size(Ses, 2)));
            Length_all_s = mean(cellfun(@(x)x(2), cellfun(@size, Temporal_Traces, 'UniformOutput', false)), 2);

            % Create Data structure containing aligned cells
            Indices_loop = [0; cumsum(Length_all_s)];
            Data = struct;
            Data.Start_Frame_Session = Indices_loop(1:end-1) + 1;

            Neural_S = nan(All_Ses, sum(Length_all_s));
            Neural_C = nan(All_Ses, sum(Length_all_s));
            Neural_C_raw = nan(All_Ses, sum(Length_all_s));
            Neural_A = nan(size(OutStruct.objectMapTurboreg, 1), size(OutStruct.objectMapTurboreg, 2), All_Ses, size(All_nam{i}, 2));

            Indices_aligned = OutStruct.globalIDs(Num_Ses_aligned(:,end) == size(OutStruct.globalIDs, 2), :);

            for j = 1:size(Temporal_Traces, 1)
                Neural_C(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                    (Temporal_Traces{j, 1}(Indices_aligned(:,j), :));
                Neural_S(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                    (Temporal_Traces{j, 2}(Indices_aligned(:,j), :));        
                Neural_C_raw(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                    (Temporal_Traces{j, 3}(Indices_aligned(:,j), :));
                Neural_A(:, :, :, j) = spatial_footprints{j}(:, :, Indices_aligned(:,j));
            end

            Data.S = Neural_S;
            Data.C = Neural_C;
            Data.A = squeeze(mean(Neural_A, 4));
            Data.C_Raw = Neural_C_raw;  
            Data.CN = Corr_Img;

            Neural_S_all = nan(size(index, 1), sum(Length_all_s));
            Neural_C_all = nan(size(index, 1), sum(Length_all_s));
            Neural_C_raw_all = nan(size(index, 1), sum(Length_all_s));

            for j = 1:size(Temporal_Traces, 1)
                Temp_idx = OutStruct.globalIDs(OutStruct.globalIDs(:,j) ~=0, j);
                Temp_idx2 = OutStruct.globalIDs(:,j) ~=0;
                Neural_C_all(Temp_idx2, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                    (Temporal_Traces{j, 1}(Temp_idx, :));
                Neural_S_all(Temp_idx2, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                    (Temporal_Traces{j, 2}(Temp_idx, :));
                Neural_C_raw_all(Temp_idx2, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                    (Temporal_Traces{j, 3}(Temp_idx, :));
            end

            Data.S_all = Neural_S_all;
            Data.C_all = Neural_C_all;
            Data.C_Raw_all = Neural_C_raw_all; 

            Data.C_unsorted = {Temporal_Traces{:, 1}};
            Data.S_unsorted = {Temporal_Traces{:, 2}};
            Data.C_Raw_unsorted = {Temporal_Traces{:, 3}};
            spatial_footprints_resorted = cell(size(spatial_footprints, 2), 1);

            for zz = 1:size(spatial_footprints, 2)
                spatial_footprints_resorted{zz} = permute(spatial_footprints{zz}, [3, 1, 2]);
            end

            Data.A_unsorted = spatial_footprints_resorted;
            Data.Order = OutStruct.globalIDs;
            Data.Non_registered_Ses = Non_registered_Ses;
            Data.SessionNames = All_nam{i};
            
            % Plot all the cells that were selected and Alignment

            g = figure('Position', [500 500 1500 800]);
            suptitle('All Cells')
            subplot(1,2,1)
            imagesc(Neural_C_all, [0 10]); colormap bone
            set(gca, 'YTick', [], 'XTick', []); xlabel('Time'); ylabel('# Cells')

            subplot(1,2,2)
            hold on;
            for j = 1:size(Neural_C_raw_all, 1)
                plot(Neural_C_raw_all(j,:) + 3*j, 'k')
            end
            set(gca, 'YTick', [], 'XTick', []); xlabel('Time'); 
            xlim([0 size(Neural_C_all, 2)]); ylim([0 3*j+5]);
            saveas(gcf,[results_directory '\All_Cells.png'])
            close(g)
            
            % Plot Aligned Cells
            if ~isempty(Neural_C_raw)         
                g = figure;
                suptitle('Aligned Cells')
                subplot(1,2,1)        
                imagesc(Neural_C_raw, [0 25]);colormap bone
                set(gca, 'YTick', [], 'XTick', []); xlabel('Time'); ylabel('# Cells')

                subplot(1,2,2)
                hold on
                for j = 1:size(Neural_C_raw, 1)
                    plot(Neural_C_raw(j,:) + 15*j, 'k')
                end
                xlim([0 size(Neural_C_raw, 2)]); ylim([0 15*j+15]);
                set(gca, 'YTick', [], 'XTick', []); xlabel('Time');
                saveas(gcf,[results_directory '\Aligned_Cells.png'])
                close(g)
            end
            
            % Plot Spatial footptints of aligned Cells - Taken from Sheintuch
            % et al.
            figures_visibility = 'on'; % either 'on' or 'off' (in any case figures are saved)
            optimal_cell_to_index_map = OutStruct.globalIDs;
            Permuted_Spatial = cellfun(@(x) permute(x, [3, 1, 2]), spatial_footprints, 'UniformOutput', false);
            plot_all_registered_projections(Permuted_Spatial,optimal_cell_to_index_map,figures_directory,figures_visibility)

            Alignment{1} = Data;
            Alignment{2} = OutStruct.globalIDs;
            Alignment{3} = OutStruct;
            Alignment{4} = Temporal_Traces;
            Alignment{5} = Permuted_Spatial;

            savefast([results_directory '\Data_Miniscope.mat'], 'Data')
            savefast([results_directory '\Alignment.mat'], 'Alignment')
            close all
        catch
            disp('Alignment failed, probably due to low cell yield and insufficient data to model.') 

            Data.C_unsorted = {Temporal_Traces{:, 1}};
            Data.S_unsorted = {Temporal_Traces{:, 2}};
            Data.C_Raw_unsorted = {Temporal_Traces{:, 3}};
            Data.A_unsorted = spatial_footprints;
            savefast([results_directory '\Data_Miniscope.mat'], 'Data')
            close all
        end        
    end
end








