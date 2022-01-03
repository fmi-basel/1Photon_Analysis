%% Post processing of CNMFE components
% Written by Julian Hinz @ Luthi lab 2018
% Based on scripts from Pengcheng Zhou, Columbia University
% Used to process results from the CNMFE algorithm
% https://elifesciences.org/articles/28728

function post_process_CNMFE(All_nam, Threshold_Low_SN, Overlap, Outer_area, ...
    minPixel)

    for i = 1:size(All_nam, 1)
        
        % EXtract Sessions 
        Sessions = All_nam{i};
        
        for ii = 1:size(Sessions, 2)

            load([All_nam{i}{ii} '\processed_data\Result_CNMFE.mat']);
            load([All_nam{i}{ii} '\processed_data\CN.mat']);        
            
            % find ROIs that lie exclusively in the outer area of the FOV
            % and exclude any small neurons
            nums = size(neuron.A, 2);

            ROI_Shape = reshape(full(neuron.A), size(Cn, 1), size(Cn, 2), nums);

            Exclude = any([reshape(ROI_Shape(1:Outer_area,: , :), Outer_area*size(Cn, 2), nums); ...
                reshape(ROI_Shape(end-Outer_area+1:end,: , :), Outer_area*size(Cn, 2), nums); ...
                reshape(ROI_Shape(:, 1:Outer_area, :), Outer_area*size(Cn, 1), nums); ...
                reshape(ROI_Shape(:, end-Outer_area+1:end, :), Outer_area*size(Cn, 1), nums)] > 0) ...
                | sum(full(neuron.A) > 1) <= minPixel;
            
            %%%%%% Exclude Cells that are overlapping
            % Calculate Overlap of components for the low SNR neurons
            Num_Neurons = size(neuron.A, 2);
            Spatial_Components = full(neuron.A);
            Weighted_Components = double(Spatial_Components./max(Spatial_Components) > 0.5);

            Delete = zeros(Num_Neurons, 1);
            % run only through low Signal to noise neurons

            for k = floor(Num_Neurons*Threshold_Low_SN):Num_Neurons

                Size_Comp = sum(Weighted_Components(:, k)); % Size of component
                Temp_comp = Weighted_Components - repmat(Weighted_Components(:, k), 1, Num_Neurons); % subtract location of the component in question from all the others
                              
                Overlaps = abs(sum(Temp_comp < 0) - Size_Comp); % Now 'normalize' for size    
                % If more than x, specified by input, overlap with other
                % components exclude it 
                if sum(Overlaps(1:k-1)) > ceil(Overlap * Size_Comp)
                    Delete(k) = 1;                    
                end
            end
    
            Exclude_comb = logical(Delete) | Exclude';
            savefast([All_nam{i}{ii} '\processed_data\Delete_Overlap_1.mat'], 'Exclude_comb');
            
            % Utilize the tool from the CNMFE with adaptations to manually
            % exclude neurons that don't fit selection criteria
            neuron.options.spatial_algorithm = 'nnls';
            neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
            neuron.viewNeurons2([], Cn, Exclude_comb, neuron.C_raw);
            
            if isempty(neuron.A) 
                Dummyvar = zeros(10,10);
                savefast([All_nam{i}{ii} '\processed_data\Result_CNMFE_noResults.mat'], 'Dummyvar');
                disp('All Neurons were deleted, a dummy variable saved. The alignment will skip this session later.')
            else    
                savefast([All_nam{i}{ii} '\processed_data\Result_CNMFE_postprocessed.mat'], 'neuron');
                
                % Calculate Overlap of components for the low SNR neurons
                Num_Neurons = size(neuron.A, 2);
                Spatial_Components = full(neuron.A);
                Weighted_Components = double(Spatial_Components./max(Spatial_Components) > 0.5);

                Delete = zeros(Num_Neurons, 1);
                
                % run only through low Signal to noise neurons
                for k = floor(Num_Neurons*Threshold_Low_SN):Num_Neurons

                    Size_Comp = sum(Weighted_Components(:, k)); % Size of component
                    Temp_comp = Weighted_Components - repmat(Weighted_Components(:, k), 1, Num_Neurons); % subtract location of the component in question from all the others

                    Overlaps = abs(sum(Temp_comp < 0) - Size_Comp); % Now 'normalize' for size    

                    % If more than x, specified by input, overlap with other
                    % components exclude it 
                    if sum(Overlaps(1:k-1)) > ceil(Overlap * Size_Comp)
                        Delete(k) = 1;                    
                    end
                end

                savefast([All_nam{i}{ii} '\processed_data\Delete_Overlap_2.mat'], 'Delete');
                
                % Summarize Data in new structure
                Data.A = neuron.A(:,~logical(Delete)); 
                Data.A_reshaped = reshape(Weighted_Components, ...
                    [size(Cn, 1) size(Cn, 2), size(Weighted_Components, 2)]);
                Data.A_weighted = sparse(Weighted_Components(:,~logical(Delete)));
                Data.C = neuron.C(~logical(Delete),:);
                Data.C_raw = neuron.C_raw(~logical(Delete),:);
                Data.S = neuron.S(~logical(Delete),:);

                savefast([All_nam{i}{ii} '\processed_data\CNMFE_Data.mat'], 'Data');

                % Get a new Img of the Outlines
                neuron.A = neuron.A(:,~logical(Delete));            
                Coor = neuron.show_contours(0.6);
                saveas(gcf,[All_nam{i}{ii} '\processed_data\Outlines_Neurons.png'])
            end
            
            close all
            clear neuron
            
            [~,Name,~] = fileparts(All_nam{i}{ii});
            
            disp('******************************************************')
            fprintf(['Finished with folder ' Name '\n'])
            fprintf('This is Session %d out of %d sessions of this animal.\n', ...
                ii, size(All_nam{i}, 2))
            disp('******************************************************')
        end
        disp('******************************************************')
        fprintf('You now finished animal %d out of %d.\n', i, size(All_nam, 1))
        disp('******************************************************')
    end

end
