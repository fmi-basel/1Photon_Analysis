%% Alignment over days
function [Output_Direc] = Sheintuch_registration_over_days(All_nam, imaging_technique, length_identifier, Folder_Structure, ...
    p_same_certainty_threshold, p_same_threshold, registration_approach, initial_registration_type, Threshold_manual, ...
    Registration_type2, Threshold_manual2, Model_Selection, reference_session_index)

% Method overview:
% This is an implementation of a probabilistic approach for the
% identification of the same neurons (cell registration) across multiple sessions
% in Ca2+ imaging data, developed by Sheintuch et al., 2017.

% Input: The inputs for the cell registration method are the spatial footprints of
% cellular activity (weighted ROIs) of the cells that were detected in the different
% sessions. Each spatial footprint is a matrix the size of the frame and
% each pixel's value represents its contribution to the
% cell's fluorescence.

% Output: The main output for the cell registration method is the obtained mapping of
% cell identity across all registered sessions. It is a matrix the size of
% the final number of registered cells by the number of registered
% sessions. Each entry holds the index for the cell in a given session.
% Other outputs include:
% 1. register scores - providing with the registration quality of each cell register
% 2. log file - with all the relevant information regarding the data, registration
% configurations, and a summary of the registration results and quality.
% 3. figures - saved automatically in a designated folder.

% The code includes the following stages:
% 1. Loading the spatial footprints of cellular activity from the different sessions.
% 2. Aligning all the sessions according to a reference coordinate system
% 3. Computing a probabilistic model of the spatial footprints similarities
% 4. Obtaining an initial cell registration according to an optimized registration threshold.
% 5. Obtaining the final cell registration based on a clustering algorithm.

% Setting paths for the cell registration procedure:
Output_Direc = cell(size(All_nam, 1), 1);

% Defining the parameters for image alignment:
alignment_type = 'Translations and Rotations'; % either 'Translations' or 'Translations and Rotations'
use_parallel_processing = true; % either true or false
maximal_rotation = 10; % in degrees - only relevant if 'Translations and Rotations' is used    

for i = 1:size(All_nam, 1)
    
    Sessions = All_nam{i};

    spatial_footprints = cell(1, size(Sessions, 2));
    Temporal_Traces =  cell(size(Sessions, 2), 3);
    [Savepath, Foldern ,~] = fileparts(All_nam{i}{1});
    Alignment = cell(6);
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
                break
            end
        end
               
        try
            spatial_footprints{ii} =  permute(reshape(full(Data.A), [size(Data.A_reshaped, 1) ...
                size(Data.A_reshaped, 2) size(Data.A, 2)]), [3 1 2]);
        catch  
            Spatial_Components = full(Data.A);
            Weighted_Components = double(Spatial_Components./max(...
                Spatial_Components) > 0.5);
            Data.A_reshaped = reshape(Weighted_Components, ...
                [size(Cn, 1) size(Cn, 2), size(Weighted_Components, 2)]);         
            spatial_footprints{ii} =  permute(reshape(full(Data.A), [size(Data.A_reshaped, 1) ...
                size(Data.A_reshaped, 2) size(Data.A, 2)]), [3 1 2]);
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
            results_directory = [Savepath '\Alignment\S_' Foldern(Index:Index+length_identifier - 1) '\'];
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

    %% Stage 1 - Loading the spatial footprints of cellular activity:
    % This stage loads a new data set which includes several sessions with the
    % identified spatial footprints.

    % Defining the parameters:
    microns_per_pixel=2.4; % Parameter to be determined !

    % Loading the data:
    disp('Stage 1 - Loading sessions')
    % [spatial_footprints,number_of_sessions]=load_multiple_sessions(file_names);
    [footprints_projections]=compute_footprints_projections(spatial_footprints);
    plot_all_sessions_projections(footprints_projections,figures_directory,figures_visibility)
    disp('Done')
    
    %% Stage 2 - Aligning all the sessions to a reference coordinate system:
    % A rigid-body transfomration is applied to all the sessions
    % according to a chosen reference ssseion. The alignment includes:
    % 1. Preparing the data for alignment
    % 2. Aligning all the sessions according to a reference coordinate system
    % 3. Evaluating how suitable the data is for longitudinal analysis


    try
        % Preparing the data for alignment:
        disp('Stage 2 - Aligning sessions')
        [normalized_spatial_footprints]=normalize_spatial_footprints(spatial_footprints);
        [adjusted_spatial_footprints,adjusted_FOV,adjusted_x_size,adjusted_y_size,adjustment_zero_padding]=...
            adjust_FOV_size(normalized_spatial_footprints);
        [adjusted_footprints_projections]=compute_footprints_projections(adjusted_spatial_footprints);
        [centroid_locations]=compute_centroid_locations(adjusted_spatial_footprints,microns_per_pixel); 
        [centroid_projections]=compute_centroids_projections(centroid_locations,adjusted_spatial_footprints);

        % Aligning the cells according to the tranlations/rotations that maximize their similarity:
        if strcmp(alignment_type,'Translations and Rotations')
            [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV]=...
                align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,use_parallel_processing,maximal_rotation);
        elseif strcmp(alignment_type,'Translations')
            [spatial_footprints_corrected,centroid_locations_corrected,footprints_projections_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,overlapping_FOV]=...
                align_images(adjusted_spatial_footprints,centroid_locations,adjusted_footprints_projections,centroid_projections,adjusted_FOV,microns_per_pixel,reference_session_index,alignment_type,use_parallel_processing);
        end

        % Evaluating data quality:
        [all_centroid_projections_correlations,number_of_cells_per_session]=...
            evaluate_data_quality(spatial_footprints_corrected,centroid_projections_corrected,maximal_cross_correlation,alignment_translations,reference_session_index,alignment_type);
        plot_alignment_results(adjusted_spatial_footprints,centroid_locations,spatial_footprints_corrected,centroid_locations_corrected,adjusted_footprints_projections,footprints_projections_corrected,reference_session_index,all_centroid_projections_correlations,maximal_cross_correlation,alignment_translations,overlapping_FOV,alignment_type,number_of_cells_per_session,figures_directory,figures_visibility)

        disp('Done')

        %% Stage 3 (part a) - Calculating the similarities distributions form the data:
        % This stage uses the ditribtuions of centroid distance and spatial correlations
        % to compute the probabilities of neighboring cell-pairs to be the same cell (P_same).

        % part a includes the calculation of the distributions of centroid distances and spatial
        % correlations from the data.

        % Defining the parameters for the probabilstic modeling:
        maximal_distance = 12; % cell-pairs that are more than 12 micrometers apart are assumed to be different cells
        normalized_maximal_distance = maximal_distance/microns_per_pixel;

        % Computing correlations and distances across days:
        disp('Stage 3 - Calculating a probabilistic model of the data')
        [number_of_bins,centers_of_bins]=estimate_number_of_bins(spatial_footprints,normalized_maximal_distance);
        [all_to_all_indexes,all_to_all_spatial_correlations,all_to_all_centroid_distances,neighbors_spatial_correlations,neighbors_centroid_distances,neighbors_x_displacements,neighbors_y_displacements,NN_spatial_correlations,NNN_spatial_correlations,NN_centroid_distances,NNN_centroid_distances]=...
            compute_data_distribution(spatial_footprints_corrected,centroid_locations_corrected,normalized_maximal_distance);

        % Plotting the (x,y) displacements:
        plot_x_y_displacements(neighbors_x_displacements,neighbors_y_displacements,microns_per_pixel,normalized_maximal_distance,number_of_bins,centers_of_bins,figures_directory,figures_visibility);
        disp('Part a done')

        %% Stage 3 (part b) - Compute a probabilistic model:
        % Modeling the data as a weighted sum of same cells and different cells,
        % and estimating the attainable registration accuracy:

        disp('Calculating a probabilistic model of the data')
        % Modeling the distribution of centroid distances:
    
        [centroid_distances_model_parameters,p_same_given_centroid_distance,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,MSE_centroid_distances_model,centroid_distance_intersection]=...
            compute_centroid_distances_model(neighbors_centroid_distances,microns_per_pixel,centers_of_bins);

 
        % Modeling the distribution of spatial correlations:
        if strcmp(imaging_technique,'one_photon')
            [spatial_correlations_model_parameters,p_same_given_spatial_correlation,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,MSE_spatial_correlations_model,spatial_correlation_intersection]=...
                compute_spatial_correlations_model(neighbors_spatial_correlations,centers_of_bins);
        end

        % estimating registration accuracy:
        if strcmp(imaging_technique,'one_photon')
            [p_same_centers_of_bins,uncertain_fraction_centroid_distances,cdf_p_same_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold,uncertain_fraction_spatial_correlations,cdf_p_same_spatial_correlations,false_positive_per_correlation_threshold,true_positive_per_correlation_threshold]=...
                estimate_registration_accuracy(p_same_certainty_threshold,neighbors_centroid_distances,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,centers_of_bins,neighbors_spatial_correlations,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,p_same_given_spatial_correlation);
            % Checking which model is better according to a defined cost function:
            [best_model_string]=choose_best_model(uncertain_fraction_centroid_distances,MSE_centroid_distances_model,imaging_technique,uncertain_fraction_spatial_correlations,MSE_spatial_correlations_model);
        else
            [p_same_centers_of_bins,uncertain_fraction_centroid_distances,cdf_p_same_centroid_distances,false_positive_per_distance_threshold,true_positive_per_distance_threshold]=...
                estimate_registration_accuracy(p_same_certainty_threshold,neighbors_centroid_distances,centroid_distances_model_same_cells,centroid_distances_model_different_cells,p_same_given_centroid_distance,centers_of_bins);
            [best_model_string]=choose_best_model(uncertain_fraction_centroid_distances,MSE_centroid_distances_model,imaging_technique);
        end

        % Plotting the probabilistic models and estimated registration accuracy:
        if strcmp(imaging_technique,'one_photon')
            plot_models(centroid_distances_model_parameters,NN_centroid_distances,NNN_centroid_distances,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,centroid_distance_intersection,centers_of_bins,microns_per_pixel,normalized_maximal_distance,figures_directory,figures_visibility,spatial_correlations_model_parameters,NN_spatial_correlations,NNN_spatial_correlations,spatial_correlations_distribution,spatial_correlations_model_same_cells,spatial_correlations_model_different_cells,spatial_correlations_model_weighted_sum,spatial_correlation_intersection)
            plot_estimated_registration_accuracy(p_same_centers_of_bins,p_same_certainty_threshold,p_same_given_centroid_distance,centroid_distances_distribution,cdf_p_same_centroid_distances,uncertain_fraction_centroid_distances,true_positive_per_distance_threshold,false_positive_per_distance_threshold,centers_of_bins,normalized_maximal_distance,microns_per_pixel,imaging_technique,figures_directory,figures_visibility,p_same_given_spatial_correlation,spatial_correlations_distribution,cdf_p_same_spatial_correlations,uncertain_fraction_spatial_correlations,true_positive_per_correlation_threshold,false_positive_per_correlation_threshold)
        else
            plot_models(centroid_distances_model_parameters,NN_centroid_distances,NNN_centroid_distances,centroid_distances_distribution,centroid_distances_model_same_cells,centroid_distances_model_different_cells,centroid_distances_model_weighted_sum,centroid_distance_intersection,centers_of_bins,microns_per_pixel,normalized_maximal_distance,figures_directory,figures_visibility)
            plot_estimated_registration_accuracy(p_same_centers_of_bins,p_same_certainty_threshold,p_same_given_centroid_distance,centroid_distances_distribution,cdf_p_same_centroid_distances,uncertain_fraction_centroid_distances,true_positive_per_distance_threshold,false_positive_per_distance_threshold,centers_of_bins,normalized_maximal_distance,microns_per_pixel,imaging_technique,figures_directory,figures_visibility)
        end

        % Computing the P_same for each neighboring cell-pair according to the different models:
        if strcmp(imaging_technique,'one_photon')
            [all_to_all_p_same_centroid_distance_model,all_to_all_p_same_spatial_correlation_model]=...
                compute_p_same(all_to_all_centroid_distances,p_same_given_centroid_distance,centers_of_bins,imaging_technique,all_to_all_spatial_correlations,p_same_given_spatial_correlation);
        else
            [all_to_all_p_same_centroid_distance_model]=...
                compute_p_same(all_to_all_centroid_distances,p_same_given_centroid_distance,centers_of_bins,imaging_technique);
        end
        disp('Done')

        %% Stage 4 - Initial cell registration
        % This stage performs an initial cell registration according to an
        % optimized threshold of either spatial correlations or centroid distances.

        % Defining the parameters for initial registration:
        
        % initial_registration_type=best_model_string; % either 'Spatial correlation', 'Centroid distance', or 'best_model_string';
        
        % The threshold that corresponds to p_same=0.5 is automatically chosen.
        % if a specific distance/correlation threshold is to be used - change the
        % initial threshold manually in the next few lines.
        
        if Model_Selection
            if exist('spatial_correlation_intersection','var')
                initial_threshold=spatial_correlation_intersection; % the threshold for p_same=0.5;
            end
            if exist('centroid_distance_intersection','var')
                initial_threshold=centroid_distance_intersection; % the threshold for p_same=0.5;
            end
        else
            initial_threshold = Threshold_manual;
        end

        % Computing the initial registration according to a simple threshold:
        disp('Stage 4 - Performing initial registration')
        if strcmp(initial_registration_type,'Spatial correlation') % if spatial correlations are used
%             if exist('spatial_correlation_intersection','var')
%                 initial_threshold=spatial_correlation_intersection; % the threshold for p_same=0.5;
%             else
%                 initial_threshold=0.65; % a fixed correlation threshold not based on the model
%             end
            if strcmp(imaging_technique,'two_photon')
                error('The spatial correlations model is only applied for 1-photon imaging data')
            else
                [cell_to_index_map,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations]=...
                    initial_registration_spatial_correlations(normalized_maximal_distance,initial_threshold,spatial_footprints_corrected,centroid_locations_corrected);
                plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints_corrected,initial_registration_type,figures_directory,figures_visibility,registered_cells_spatial_correlations,non_registered_cells_spatial_correlations)
            end
        else % if centroid distances are used
%             if exist('centroid_distance_intersection','var')
%                 initial_threshold=centroid_distance_intersection; % the threshold for p_same=0.5;
%             else
%                 initial_threshold=5; % a fixed distance threshold not based on the model
%             end
            normalized_distance_threshold=initial_threshold/microns_per_pixel;
            [cell_to_index_map,registered_cells_centroid_distances,non_registered_cells_centroid_distances]=...
                initial_registration_centroid_distances(normalized_maximal_distance,normalized_distance_threshold,centroid_locations_corrected);
            plot_initial_registration(cell_to_index_map,number_of_bins,spatial_footprints_corrected,initial_registration_type,figures_directory,figures_visibility,registered_cells_centroid_distances,non_registered_cells_centroid_distances,microns_per_pixel,normalized_maximal_distance)
        end

        disp([num2str(size(cell_to_index_map,1)) ' cells were found'])
        disp('Done')

        %% Stage 5 - Final cell registration:
        % This stage performs the final cell registration with a clustering algorithm 
        % that is based on the probability model for same cells and different cells. 
        % P_same can be either according to centroid distances or spatial
        % correlations.

        % Defining the parameters for final registration:
        if isempty(Registration_type2)
            model_type=best_model_string; % either 'Spatial correlation' or 'Centroid distance'
        else
            model_type = Registration_type2;            
        end
        
        % Deciding on the registration threshold:
        transform_data=false;
        if strcmp(registration_approach,'Simple threshold') % only relevant if a simple threshold is used
            if strcmp(model_type,'Spatial correlation')
                
                if Model_Selection
                    if exist('spatial_correlation_intersection','var')
                        final_threshold=spatial_correlation_intersection; % the threshold for p_same=0.5;
                    else
                        final_threshold=0.65; % a fixed correlation threshold not based on the model
                    end
                else
                   final_threshold = Threshold_manual2;  
                end
                
            elseif strcmp(model_type,'Centroid distance')
                
                if Model_Selection
                    if exist('centroid_distance_intersection','var')
                        final_threshold=centroid_distance_intersection; % the threshold for p_same=0.5;
                    else
                        final_threshold=5; % a fixed distance threshold not based on the model
                    end
                else
                    final_threshold = Threshold_manual2;
                end
                
                normalized_distance_threshold=(maximal_distance-final_threshold)/maximal_distance;
                transform_data=true;
            end
        else
            final_threshold=p_same_threshold;
        end

        % Registering the cells with the clustering algorithm:
        disp('Stage 5 - Performing final registration')
        if strcmp(registration_approach,'Probabilistic')    
            if strcmp(model_type,'Spatial correlation')
                [optimal_cell_to_index_map,registered_cells_centroids,cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=...
                    cluster_cells(cell_to_index_map,all_to_all_p_same_spatial_correlation_model,all_to_all_indexes,normalized_maximal_distance,p_same_threshold,centroid_locations_corrected,registration_approach,transform_data);
            elseif strcmp(model_type,'Centroid distance')
                [optimal_cell_to_index_map,registered_cells_centroids,cell_scores,cell_scores_positive,cell_scores_negative,cell_scores_exclusive,p_same_registered_pairs]=...
                    cluster_cells(cell_to_index_map,all_to_all_p_same_centroid_distance_model,all_to_all_indexes,normalized_maximal_distance,p_same_threshold,centroid_locations_corrected,registration_approach,transform_data);
            end
            plot_cell_scores(cell_scores_positive,cell_scores_negative,cell_scores_exclusive,cell_scores,p_same_registered_pairs,figures_directory,figures_visibility)
        elseif strcmp(registration_approach,'Simple threshold')
            if strcmp(model_type,'Spatial correlation')
                [optimal_cell_to_index_map,registered_cells_centroids]=...
                    cluster_cells(cell_to_index_map,all_to_all_spatial_correlations,all_to_all_indexes,normalized_maximal_distance,final_threshold,centroid_locations_corrected,registration_approach,transform_data);
            elseif strcmp(model_type,'Centroid distance')
                [optimal_cell_to_index_map,registered_cells_centroids]=...
                    cluster_cells(cell_to_index_map,all_to_all_centroid_distances,all_to_all_indexes,normalized_maximal_distance,normalized_distance_threshold,centroid_locations_corrected,registration_approach,transform_data);
            end
        end
        [is_in_overlapping_FOV]=check_if_in_overlapping_FOV(registered_cells_centroids,overlapping_FOV);

        % Plotting the registration results with the cell maps from all sessions:
        plot_all_registered_projections(spatial_footprints_corrected,optimal_cell_to_index_map,figures_directory,figures_visibility)

        % saving the final registration results:
        disp('Saving the results')
        cell_registered_struct=struct;
        cell_registered_struct.cell_to_index_map=optimal_cell_to_index_map;
        if strcmp(registration_approach,'Probabilistic')
            cell_registered_struct.cell_scores=cell_scores';
            cell_registered_struct.true_positive_scores=cell_scores_positive';
            cell_registered_struct.true_negative_scores=cell_scores_negative';
            cell_registered_struct.exclusivity_scores=cell_scores_exclusive';
            cell_registered_struct.p_same_registered_pairs=p_same_registered_pairs';
        end
        cell_registered_struct.is_cell_in_overlapping_FOV=is_in_overlapping_FOV';
        cell_registered_struct.registered_cells_centroids=registered_cells_centroids';
        cell_registered_struct.centroid_locations_corrected=centroid_locations_corrected';
        cell_registered_struct.spatial_footprints_corrected=spatial_footprints_corrected';
        cell_registered_struct.spatial_footprints_corrected=spatial_footprints_corrected';
        cell_registered_struct.alignment_x_translations=alignment_translations(1,:);
        cell_registered_struct.alignment_y_translations=alignment_translations(2,:);
        if strcmp(alignment_type,'Translations and Rotations')
            cell_registered_struct.alignment_rotations=alignment_translations(3,:);
        end
        cell_registered_struct.adjustment_x_zero_padding=adjustment_zero_padding(1,:);
        cell_registered_struct.adjustment_y_zero_padding=adjustment_zero_padding(2,:);
        savefast(fullfile(results_directory,['cellRegistered_' datestr(clock,'yyyymmdd_HHMMss') '.mat']),'cell_registered_struct')

        % Saving a log file with all the chosen parameters:

        disp([num2str(size(optimal_cell_to_index_map,1)) ' cells were found'])
        disp('End of cell registration procedure')

         %%% Post Processing 
        Num_Ses_aligned = cumsum(transpose(cell_to_index_map~=0))';

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
        Neural_A = nan(All_Ses, size(overlapping_FOV, 1), size(overlapping_FOV, 2), size(All_nam{i}, 2));

        Indices_aligned = cell_to_index_map(Num_Ses_aligned(:,end) == size(cell_to_index_map, 2), :);

        for j = 1:size(Temporal_Traces, 1)
            Neural_C(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Temporal_Traces{j, 1}(Indices_aligned(:,j), :));
            Neural_S(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Temporal_Traces{j, 2}(Indices_aligned(:,j), :));        
            Neural_C_raw(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Temporal_Traces{j, 3}(Indices_aligned(:,j), :));
            Neural_A(:, :, :, j) = spatial_footprints{j}(Indices_aligned(:,j), :, :);
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
            Temp_idx = cell_to_index_map(cell_to_index_map(:,j) ~=0, j);
            Temp_idx2 = cell_to_index_map(:,j) ~=0;
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
        Data.A_unsorted = spatial_footprints;
        Data.Order = cell_to_index_map;
        Data.Non_registered_Ses = Non_registered_Ses;

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
        xlim([0 size(Neural_C_raw, 2)]); ylim([0 3*j+5]);
        saveas(gcf,[results_directory '\All_Cells.png'])
        close(g)
        
        % Plot Aligned Cells
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
        
        Alignment{1} = Data;
        Alignment{2} = cell_to_index_map;
        Alignment{3} = cell_registered_struct;
        Alignment{4} = adjusted_FOV;
        Alignment{5} = Temporal_Traces;
        Alignment{6} = spatial_footprints;

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

if use_parallel_processing
    delete(gcp);
end

end