%% Extract neural signal with CNMFE

function CNMFE(All_nam, Fs, NW2Test)
    

    for i = 1:size(All_nam, 2)

        nam_temp = All_nam{i};         
        Filenames = dir(fullfile([nam_temp '\processed_data'], '*.mat'));
        index2 = cellfun(@(x) any(strfind(x, 'MC.mat')), {Filenames.name});

        %% choose multiple datasets or just one  
        neuron = Sources2D(); 
        nam = [nam_temp '\processed_data\' Filenames(index2).name];         % you can put all file names into a cell array; when it's empty, manually select files 
        nam = neuron.select_data(nam);  %if nam is [], then select data interactively
        
        pars_envs = struct('memory_size_to_use', 150 / 12, ...   % GB, memory space you allow to use in MATLAB
            'memory_size_per_patch', 4, ...   % GB, space for loading data within one patch
            'patch_dims', [64, 64]);  %GB, patch size
        
        %% parameters
        % -------------------------    COMPUTATION    -------------------------  %

        % -------------------------      SPATIAL      -------------------------  %
        gSig = 4;           % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
        gSiz = 16;          % pixel, neuron diameter
        ssub = 1;           % spatial downsampling factor
        with_dendrites = true;   % with dendrites or not
        if with_dendrites
            % determine the search locations by dilating the current neuron shapes
            updateA_search_method = 'dilate';  % #ok<UNRCH>
            updateA_bSiz = 5;
            updateA_dist = neuron.options.dist;
        else
            % determine the search locations by selecting a round area
            updateA_search_method = 'ellipse'; %#ok<UNRCH>
            updateA_dist = 5;
            updateA_bSiz = neuron.options.dist;
        end
        spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'
        spatial_algorithm = 'hals_thresh';

        % -------------------------      TEMPORAL     -------------------------  %
        tsub = 1;           % temporal downsampling factor
        deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
            'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
            'smin', -5, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level
            'optimize_pars', true, ...  % optimize AR coefficients
            'optimize_b', true, ...% optimize the baseline);
            'max_tau', 100);    % maximum decay time (unit: frame);

        nk = 3;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
        % when changed, try some integers smaller than total_frame/(Fs*30)
        detrend_method = 'spline';  % compute the local minimum as an estimation of trend.

        % -------------------------     BACKGROUND    -------------------------  %
        bg_model = 'ring';  % model of the background {'ring', 'svd'(default), 'nmf'}
        nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
        ring_radius = 18;  % when the ring model used, it is the radius of the ring used in the background model.
        %otherwise, it's just the width of the overlapping area
        num_neighbors = []; % number of neighbors for each neuron

        % -------------------------      MERGING      -------------------------  %
        show_merge = false;  % if true, manually verify the merging step
        merge_thr = 0.5;  % old: 0. 65, 0.5   % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
        method_dist = 'mean';   % method for computing neuron distances {'mean', 'max'}
        dmin = 5; % old: 5  % minimum distances between two neurons. it is used together with merge_thr
        dmin_only = 4; % old: 2, 4, 6 % merge neurons if their distances are smaller than dmin_only.
        merge_thr_spatial = [0.6, 0.3, -inf]; % old: 0.8, 0.6 % merge components with highly correlated spatial shapes (corr=0.8) and small temporal correlations (corr=0.1)

        % -------------------------  INITIALIZATION   -------------------------  %
        K = [];             % maximum number of neurons per patch. when K=[], take as many as possible.
        min_corr = 0.8;     % minimum local correlation for a seeding pixel
        min_pnr = 8;       % minimum peak-to-noise ratio for a seeding pixel
        min_pixel = gSig^2;      % minimum number of nonzero pixels for each neuron
        bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
        frame_range = [];   % when [], uses all frames
        save_initialization = false;    % save the initialization procedure as a video.
            % use parallel computation for parallel computing
        show_init = false;   % show initialization results
        choose_params = false; % manually choose parameters
        center_psf = true;  % set the value as true when the background fluctuation is large (usually 1p data)
        % set the value as false when the background fluctuation is small (2p)

        % -------------------------  Residual   -------------------------  %
        min_corr_res = 0.7;
        min_pnr_res = 6;
        seed_method_res = 'auto';  % method for initializing neurons from the residual
        update_sn = true;

        % ----------------------  WITH MANUAL INTERVENTION  --------------------  %
        with_manual_intervention = false;

        % -------------------------  FINAL RESULTS   -------------------------  %
        save_demixed = true;    % save the demixed file or not
        kt = 3;                 % frame intervals

        % -------------------------    UPDATE ALL    -------------------------  %
        neuron.updateParams('gSig', gSig, ...       % -------- spatial --------
            'gSiz', gSiz, ...
            'ring_radius', ring_radius, ...
            'ssub', ssub, ...
            'search_method', updateA_search_method, ...
            'bSiz', updateA_bSiz, ...
            'dist', updateA_bSiz, ...
            'spatial_constraints', spatial_constraints, ...
            'spatial_algorithm', spatial_algorithm, ...
            'tsub', tsub, ...                       % -------- temporal --------
            'deconv_options', deconv_options, ...
            'nk', nk, ...
            'detrend_method', detrend_method, ...
            'background_model', bg_model, ...       % -------- background --------
            'nb', nb, ...
            'ring_radius', ring_radius, ...
            'num_neighbors', num_neighbors, ...
            'merge_thr', merge_thr, ...             % -------- merging ---------
            'dmin', dmin, ...
            'method_dist', method_dist, ...
            'min_corr', min_corr, ...               % ----- initialization -----
            'min_pnr', min_pnr, ...
            'min_pixel', min_pixel, ...
            'bd', bd, ...
            'center_psf', center_psf);
            neuron.Fs = Fs;

         %% distribute data and be ready to run source extraction
        neuron.getReady(pars_envs);
        
        % Deal with big files, first determine how many cores are available
        % There are currently two places that the Code will be run, Rig1 or
        % the VM. If it runs on the VM, select different numbers of cores
        % than with Rig 1. Then catch potential problems with RAM and 
        % reinitiate with fewer cores.
        
        c = parcluster('local'); % Create cluster
        nw = c.NumWorkers; % get Number of workers    
        
        % Set Big file processing true if file is above certain size
        b = matfile(nam);
        [d1, d2, d3] = size(b.MC);
        Size_mov = d1*d2*d3;

        if Size_mov <  3.5581e+9
            Big_Files = true;
        else
            Big_Files = false;
        end  

        % Pre specified number 
        if ~isempty(NW2Test)
            NW2Testa = NW2Test;
        elseif nw > 20 && ~Big_Files
            NW2Testa = [20 10 4 0];
        elseif nw < 20 && nw > 10 && ~Big_Files
            NW2Testa = [10 6 4 0];
        elseif nw < 10 && ~Big_Files
            NW2Testa = [4 2 0];
        elseif Big_Files
            NW2Testa = [6 4 0];
        else
            NW2Testa = [2 0];
        end

        % Now initiate a while loop that will iterate 4 times with
        % different number of cores
        
        Processing_finished = false;
        Runs = 1;
        
        while ~Processing_finished && Runs <= size(NW2Testa, 2)                       
            
            if NW2Testa(Runs) > 0
                
               if ~isempty(gcp('nocreate'))
                    delete(gcp)
               end
                
               parpool(NW2Testa(Runs));
               use_parallel = true;
            else
               if ~isempty(gcp('nocreate'))
                    delete(gcp)
               end
               
               use_parallel = false;               
            end
            
            try                 
                 [center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel);
                 Processing_finished = true;
            catch
                disp('Initiallization didnt work, lets try with fewer cores.')
                Runs = Runs + 1;
                continue
            end
            
        end
        
        
        % Continue with script
        neuron.compactSpatial();
        
        if show_init
            figure();
            ax_init= axes();
            imagesc(Cn, [0, 1]); colormap gray;
            hold on;
            plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
        end

        %% estimate the background components
        neuron.update_background_parallel(use_parallel);
        neuron_init = neuron.copy();

        %%  merge neurons and update spatial/temporal components
        neuron.merge_neurons_dist_corr(show_merge);
        neuron.merge_high_corr(show_merge, merge_thr_spatial);

        %% udpate spatial&temporal components, delete false positives and merge neurons
        % update spatial
        if update_sn
            neuron.update_spatial_parallel(use_parallel, true);
            udpate_sn = false;
        else
            neuron.update_spatial_parallel(use_parallel);
        end
        % merge neurons based on correlations 
        neuron.merge_high_corr(show_merge, merge_thr_spatial);

        for m=1:2
            % update temporal
            neuron.update_temporal_parallel(use_parallel);

            % delete bad neurons
            neuron.remove_false_positives();

            % merge neurons based on temporal correlation + distances 
            neuron.merge_neurons_dist_corr(show_merge);
        end


        %% run more iterations
        neuron.update_background_parallel(use_parallel);
        neuron.update_spatial_parallel(use_parallel);
        neuron.update_temporal_parallel(use_parallel);

        K = size(neuron.A,2);
        tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
        neuron.remove_false_positives();
        neuron.merge_neurons_dist_corr(show_merge);
        neuron.merge_high_corr(show_merge, merge_thr_spatial);

        if K~=size(neuron.A,2)
            neuron.update_spatial_parallel(use_parallel);
            neuron.update_temporal_parallel(use_parallel);
            neuron.remove_false_positives();
        end

        %% save the workspace for future analysis
        neuron.orderROIs('snr');
        cnmfe_path = neuron.save_workspace();

        %% show neuron contours
        Coor = neuron.show_contours(0.6);

        %% save neurons shapes
        neuron.save_neurons();

        savefast([nam_temp '\processed_data\Result_CNMFE.mat'], 'neuron');
        savefast([nam_temp '\processed_data\Outline.mat'], 'Coor');
        savefast([nam_temp '\processed_data\CN.mat'], 'Cn');


    end

end


%         switch Data_Format
%             case 'mat'
% 
%                 if Size_mov <  3.5581e+10
% 
%                      parpool(12);
%                      use_parallel = true;
%                 else
%                      pars_envs = struct('memory_size_to_use', 150 / 12, ...   % GB, memory space you allow to use in MATLAB
%                         'memory_size_per_patch', 4, ...   % GB, space for loading data within one patch
%                         'patch_dims', [64, 64]);  %GB, patch size
%                         use_parallel = false;
%                 end
%                 
%             case 'tif'      
%                 pars_envs = struct('memory_size_to_use', 150 / 12, ...   % GB, memory space you allow to use in MATLAB
%                     'memory_size_per_patch', 1, ...   % GB, space for loading data within one patch
%                     'patch_dims', [32, 32]);  %GB, patch size
%                 disp('Size not determined, might be too big and the RAM might blow up in your face !')
%         end
        

%% create a video for displaying the
%     amp_ac = 3000;
%     range_ac = 5+[0, amp_ac];
%     multi_factor = 20;
%     range_Y = 1e4+[0, amp_ac*multi_factor];
% 
%     avi_filename = neuron.show_demixed_video(save_demixed, kt, [], amp_ac, range_ac, range_Y, multi_factor);


%% update spatial components

%% pick neurons from the residual
% [center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel([], save_initialization, use_parallel, min_corr_res, min_pnr_res, seed_method_res);
% if show_init
%     axes(ax_init);
%     plot(center_res(:, 2), center_res(:, 1), '.g', 'markersize', 10);
% end
% neuron_init_res = neuron.copy();

%     %% initialize neurons from the video data within a selected temporal range
% %     if choose_params
% %         % change parameters for optimized initialization
% %         [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
% %     end
% 
%     %% add a manual intervention and run the whole procedure for a second time
%     neuron.options.spatial_algorithm = 'nnls';
%     if with_manual_intervention
%         show_merge = true;
%         neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
%         neuron.viewNeurons([], neuron.C_raw);
% 
%         % merge closeby neurons
%         neuron.merge_close_neighbors(true, dmin_only);
% 
%         % delete neurons
%         tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
%         ids = find(tags>0); 
%         if ~isempty(ids)
%             neuron.viewNeurons(ids, neuron.C_raw);
%         end
%     end
