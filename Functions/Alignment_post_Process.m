%% Post Processing of aligned traces
% Script to Postprocess aligned Sessions
% written by Julian Hinz @ Luthi lab 2018

function Alignment_post_Process(Output_Direc, D_aligned, Display, Fs, Animals, Registration_Methods)
    

    if isempty(Animals) || size(Output_Direc, 2) == 1
        Runs = size(Output_Direc, 1);
        A2P = 1:size(Output_Direc, 1);
    elseif ~isempty(Animals)
        A2P = Animals;
        Runs = numel(Animals);
    end

    
    for i = 1:Runs
        
        Ind_name = strfind(Output_Direc{A2P(i)}, '\');
        if exist([Output_Direc{A2P(i)} 'Alignment.mat']) ~= 0 
            load([Output_Direc{A2P(i)} 'Data_Miniscope.mat']); 
        else
            disp(['No Alignment for Folder ' Output_Direc{A2P(i)}(Ind_name(end-1)+1:Ind_name(end)-1)...
                ' took place, please check the aligment and return to this section !'])
            continue         
        end
        
        if isempty(D_aligned)
            mNum_Sessions = 2;
            Loops = size(Data.C_unsorted, 2) - mNum_Sessions + 1;
        else
            mNum_Sessions = size(Data.C_unsorted, 2) - D_aligned;
            Loops = size(Data.C_unsorted, 2) - mNum_Sessions + 1;
        end
        
        Num_aligned = fliplr(mNum_Sessions:size(Data.C_unsorted, 2));
        
        try 
            test = Data.Order;  
        catch
            disp('No Alignment !')
            continue
        end

        
        Num_Ses_aligned = cumsum(transpose(Data.Order~=0))';
        Num_Ses_aligned = Num_Ses_aligned(:,end);
        
        colors = distinguishable_colors(size(Data.C_unsorted, 2));
        ind = find(floor(sum(colors, 2)) == 0);
        
        for z = 1:numel(ind)
            colors(ind(z), :) = ones(3, 1) - 0.1 * z;
        end

        Boundries_all = cell(Loops, 1);
        Mean_Comp_all = cell(Loops, 1);
        Boundries_indv_all = cell(Loops, 1);
        
        for ii = 1:Loops
           
        Ind_2_alignedSes = Num_Ses_aligned == Num_aligned(ii);
        IND2 = find(Ind_2_alignedSes);
         
            if ~isempty(IND2)
                Boundries_loop = cell(numel(IND2), 1);
                Mean_Comp_loop = cell(numel(IND2), 1);
                Boundries_indv = cell(numel(IND2), 1);
                
                for iii = 1:numel(IND2)
                    Ord_tmp = Data.Order(IND2(iii), :);
                    ind_nonzer = find(Ord_tmp > 0, 1, 'first');
                    
                    BB = regionprops(squeeze(Data.A_unsorted{ind_nonzer}(Ord_tmp(ind_nonzer),:,:)) > 0, 'BoundingBox');
                    
                    lims_x = [BB.BoundingBox(1)-(BB.BoundingBox(3)) ...
                        BB.BoundingBox(1)+(BB.BoundingBox(3))+BB.BoundingBox(3)];
                    lims_y = [BB.BoundingBox(2)-(BB.BoundingBox(4)/2) ...
                        BB.BoundingBox(2)+(BB.BoundingBox(4))+BB.BoundingBox(4)];
                    
                    Boundries_loop{iii, :} = [lims_x lims_y]; 
                    
                    Circle = cell(size(Data.C_unsorted, 2), 1);
                    comp = cell(size(Data.C_unsorted, 2), 1);
                    
                    for p = 1:size(Data.C_unsorted, 2)

                        if ~any(Ord_tmp(p))
                            continue
                        end

                        B = bwboundaries(squeeze(Data.A_unsorted{p}(Ord_tmp(p),:,:)) > 0);
                        Circle{p} = B;
                        comp{p} = squeeze(Data.A_unsorted{p}(Ord_tmp(p),:,:));
                        

                    end
                    Boundries_indv{iii} = Circle;
                    Mean_Comp_loop{iii} = mean(cat(3, comp{:}), 3);
                end
                Mean_Comp_all{ii} = Mean_Comp_loop;
                Boundries_all{ii} = Boundries_loop;
                Boundries_indv_all{ii} = Boundries_indv;
            end         
        end        
        
        % Create Variable that can be edited
        Order_new = Data.Order;

        for ii = 1:Loops
 
            Ind_2_alignedSes = Num_Ses_aligned == Num_aligned(ii);
            IND2 = find(Ind_2_alignedSes);
            
            if ~isempty(IND2)
                iii = 1;
                
                while and(iii>=1, iii<=numel(IND2)) %for iii = 1:numel(IND2)
                    g = figure('position', [1300, 100, 1024, 1024]);
                    
                    Ord_tmp = Data.Order(IND2(iii), :);
                    ind_nonzer = find(Ord_tmp > 0, 1, 'first');
                    lims_x = Boundries_all{ii}{iii}(1:2);
                    lims_y = Boundries_all{ii}{iii}(3:4);
                    
                    for p = 1:size(Data.C_unsorted, 2)
                        
                        subplot(3, size(Data.C_unsorted, 2) + 1, p)
                        title(['Session ' num2str(p)])
                        
                        if ~any(Ord_tmp(p))
                            axis off
                            continue
                        end
                        
                        B = Boundries_indv_all{ii}{iii}{p};                        
                        hold on
                        imagesc(Data.CN{p})
                        % imagesc(squeeze(Data.A_unsorted{p}(Ord_tmp(p),:,:)))
                        colormap bone; axis off;
                        plot(B{1}(:,2), B{1}(:,1), 'color', colors(p, :))
                        xlim([0 size(squeeze(Data.A_unsorted{p}(Ord_tmp(p),:,:)), 1)]);
                        ylim([0 size(squeeze(Data.A_unsorted{p}(Ord_tmp(p),:,:)), 2)])
                        hold off
                        
                        subplot(3, size(Data.C_unsorted, 2) + 1, p + size(Data.C_unsorted, 2) + 1)
                        hold on
                        imagesc(squeeze(Data.A_unsorted{p}(Ord_tmp(p),:,:))); axis off;
                        colormap bone
                        plot(B{1}(:,2), B{1}(:,1), 'color', colors(p, :)); 
                        xlim([lims_x(1) lims_x(2)]);ylim([lims_y(1) lims_y(2)])
                        hold off
                        
                    end
                    
                    Circle = Boundries_indv_all{ii}{iii};
                    subplot(3, size(Data.C_unsorted, 2) + 1, size(Data.C_unsorted, 2) + 1)
                    title('Overlay')
                    hold on
                    % imagesc(Mean_Comp_all{ii}{iii}); 
                    imagesc(mean(cat(3,Data.CN{:}), 3)); 
                    axis off; colormap bone
                    Indices = find(~cell2mat(cellfun(@isempty, Circle, 'UniformOutput', false)));  
                    
                    for pp = 1:numel(Indices)                  
                        plot(Circle{Indices(pp)}{1}(:,2), Circle{Indices(pp)}{1}(:,1), 'color', colors(pp, :)); 
                    end
                    xlim([0 size(squeeze(Data.A_unsorted{ind_nonzer}(Ord_tmp(ind_nonzer),:,:)), 1)]);
                    ylim([0 size(squeeze(Data.A_unsorted{ind_nonzer}(Ord_tmp(ind_nonzer),:,:)), 2)])
                    hold off
                    
                    subplot(3, size(Data.C_unsorted, 2) + 1, size(Data.C_unsorted, 2)* 2 + 2)
                    hold on
                    imagesc(Mean_Comp_all{ii}{iii}); axis off; colormap bone
                        
                    for pp = 1:numel(Indices)                  
                        plot(Circle{Indices(pp)}{1}(:,2), Circle{Indices(pp)}{1}(:,1), 'color', colors(pp, :)); 
                    end
                    xlim([lims_x(1) lims_x(2)]);ylim([lims_y(1) lims_y(2)])
                    hold off
                    
                    Addit = (size(Data.C_unsorted, 2) + 1) * 2;
                    
                    subplot(3, size(Data.C_unsorted, 2) + 1, 1+Addit:size(Data.C_unsorted, 2) + Addit)
                    
                    switch Display
                        case 'Raw'
                            plot([1:size(Data.C_Raw_all, 2)]/Fs/60, Data.C_Raw_all(IND2(iii), :), 'k')
                        case 'Deconv'
                            plot([1:size(Data.C_all, 2)]/Fs/60, Data.C_all(IND2(iii), :), 'r')
                        case 'Both'
                            hold on
                            plot([1:size(Data.C_Raw_all, 2)]/Fs/60, Data.C_Raw_all(IND2(iii), :), 'k')                           
                            plot([1:size(Data.C_all, 2)]/Fs/60, Data.C_all(IND2(iii), :), 'r')
                    end
                         
                    xlabel('Time [min]'); ylabel('Activity'); xlim([0 size(Data.C_Raw_all, 2)/Fs/60 + 1]); box off
                    ylim([-0.5 nanmax(Data.C_Raw_all(IND2(iii), :)) + 1])
                    
                    vline(Data.Start_Frame_Session(2:end)/Fs/60, 'b', [], 2)
                    
                    % Now get User Input 
                    % Ask if Session is properly aligned - Yes means it is
                    % accepted as is, if it is completely wrong delete the
                    % whole order, otherwise just delete part of it based
                    % on input                  
                    
                    fprintf('Component %d, is aligned over %d of %d Sessions.\n', iii, Num_aligned(ii), size(Data.C_unsorted, 2))
                    fprintf('To accept alignment (k), to delete the whole alignment(d), to go back (b), \n to remove individual sessions type in the numbers and seperate them by spaces: ');
                    
                    temp = input('', 's');
                    
                    if any(isempty(temp) | strcmpi(temp, 'k')) 
                        iii = iii+1;
                    elseif strcmpi(temp, 'b')
                        if iii > 1
                            iii = iii-1;
                        else
                            iii = 1;
                        end
                    elseif strcmpi(temp, 'd')
                         Order_new(IND2(iii), :) = zeros(size(Data.C_unsorted, 2), 1);
                         iii = iii+1;
                    elseif strcmpi(temp, 'e')
                        break;
                    else
                        Numbers = cellfun(@str2num, regexp(temp,'\d*','Match'));
                        
                        if isempty(Numbers) || numel(Numbers) > size(Data.C_unsorted, 2)
                            disp('Unrecognized input, try again !')
                        else
                            Order_new(IND2(iii), Numbers) = zeros(numel(Numbers), 1);                                                      
                            iii = iii+1;
                        end                        
                    end
                    
                    try
                        close(g)
                    catch
                        disp('Please do not close the figure manually !')
                    end
                    
                end                
            else
                disp(['No component aligned over ' num2str(Num_aligned(ii)) ' Sessions.'])
                
            end
        end
        
        % Check if there is cells that are now not aligned over all
        % Sessions anymore
        Length_all_s = size(Data.C_all, 2); 
        
        Ind_2_alignedSes = Num_Ses_aligned == size(Data.C_unsorted, 2);
        Indices2delete = any(transpose(Order_new(Ind_2_alignedSes, :) == 0));
                       
        Data.S(Indices2delete, :) = [];
        Data.C(Indices2delete, :) = [];
        Data.C_Raw(Indices2delete, :) = [];
        Data.A(Indices2delete, :, :) = [];     
        
        % New Order File
        Order_old = Data.Order;
        Order_updated = nan(size(Order_new, 1)*2, size(Order_new, 2));
        Position = size(Order_new, 1) + 1;
        for x = 1:size(Order_new, 1)
            Changes_tmp = Order_new(x, :) ~= Order_old(x, :);
            
            if any(~Changes_tmp(~Order_old(x, :) == 0)) 
                Order_updated(x, :) = zeros(size(Order_updated, 2), 1);
                Order_updated(x, ~Changes_tmp) = Order_old(x, ~Changes_tmp);
            end
            
            if any(Changes_tmp) % Now add the single neurons to variable
                Changes = find(Changes_tmp);
                for xx = 1:numel(Changes)
                    Order_updated(Position, :) = zeros(size(Order_updated, 2), 1);
                    Order_updated(Position, Changes(xx)) = Order_old(x, Changes(xx));
                    Position = Position + 1;
                end
            end
        end
        
        Order_updated(all(isnan(Order_updated), 2), :) = [];
        Num_Ses_aligned = cumsum(transpose(Order_updated~=0))';
        [~, index] = sortrows(Num_Ses_aligned, 'descend');
        Order_updated = Order_updated(index, :);
        Num_Ses_aligned_aSort = cumsum(transpose(Order_updated~=0))';
        Data.Order = Order_updated;   
        
        % Delete the C_Raw all files etc.
        b = hist(Num_Ses_aligned(:,end), [1:size(Num_Ses_aligned, 2)]);
        g = figure;
        bar(b, 'k'); box off
        xlabel('Number of aligned Sessions'); ylabel('# Cells');
        saveas(gcf,[Output_Direc{A2P(i)} '\Figures\Alignment_Result_new.png']);
        close(g);

        [Ses, index] = sortrows(Num_Ses_aligned, 'descend');
        All_Ses = max(find(Ses(:,end) == size(Ses, 2)));

        % Create Data structure containing aligned cells
        Indices_loop = [Data.Start_Frame_Session-1; Length_all_s];
        Neural_S = nan(All_Ses, sum(Length_all_s));
        Neural_C = nan(All_Ses, sum(Length_all_s));
        Neural_C_raw = nan(All_Ses, sum(Length_all_s));
        Neural_A = nan(size(Data.A_unsorted{1}, 2), size(Data.A_unsorted{1}, 3), ...
            All_Ses, size(Data.A_unsorted, 1));

        Indices_aligned = Order_updated(Num_Ses_aligned_aSort(:,end) == size(Order_updated, 2), :);

        for j = 1:size(Data.A_unsorted, 1)
            Neural_C(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Data.C_unsorted{j}(Indices_aligned(:,j), :));
            Neural_S(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Data.S_unsorted{j}(Indices_aligned(:,j), :));        
            Neural_C_raw(:, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Data.C_Raw_unsorted{j}(Indices_aligned(:,j), :));
            Neural_A(:, :, :, j) = permute(Data.A_unsorted{j}(Indices_aligned(:,j), :, :), [2, 3, 1]);
        end

        Data.S = Neural_S;
        Data.C = Neural_C;
        Data.A = squeeze(mean(Neural_A, 4));
        Data.C_Raw = Neural_C_raw;  

        Neural_S_all = nan(size(index, 1), sum(Length_all_s));
        Neural_C_all = nan(size(index, 1), sum(Length_all_s));
        Neural_C_raw_all = nan(size(index, 1), sum(Length_all_s));
        Neural_A_all = nan(size(Data.A_unsorted{1}, 2), size(Data.A_unsorted{1}, 3), ...
            size(index, 1), size(Data.A_unsorted, 1));
        
        for j = 1:size(Data.A_unsorted, 1)
            Temp_idx = Order_updated(Order_updated(:,j) ~=0, j);
            Temp_idx2 = Order_updated(:,j) ~=0;
            Neural_C_all(Temp_idx2, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Data.C_unsorted{j}(Temp_idx, :));
            Neural_S_all(Temp_idx2, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Data.S_unsorted{j}(Temp_idx, :));
            Neural_C_raw_all(Temp_idx2, Indices_loop(j)+ 1:Indices_loop(j+1)) = ...
                (Data.C_Raw_unsorted{j}(Temp_idx, :));
            Neural_A_all(:, :, Temp_idx2, j) = permute(Data.A_unsorted{j}(Temp_idx,...
                :, :), [2, 3, 1]);
        end

        Data.S_all = Neural_S_all;
        Data.C_all = Neural_C_all;
        Data.C_Raw_all = Neural_C_raw_all; 
        Data.Neural_A_all = Neural_A_all;
        Data.A_unsorted = cellfun(@(x) permute(x, [2, 3, 1]), Data.A_unsorted, ...
            'UniformOutput', false);
        
        % Save the Output of the annotated Data
        savefast([Output_Direc{A2P(i)} 'Data_Miniscope_PP.mat'], 'Data');    
        
    end

end