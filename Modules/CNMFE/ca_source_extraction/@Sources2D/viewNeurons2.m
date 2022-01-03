function ind_del = viewNeurons2(obj, ind, Cn, del, C2, folder_nm)
%% view all components and delete components manually. it shows spatial
%   components in the full-frame and zoomed-in view. It also shows temporal
%   components
%% input:
%   ind: vector, indices of components to be displayed, no bigger than the maximum
%       number of neurons
%   C2:  K*T matrix, another temporal component to be displayed together
%       with the esitmated C. usually it is C without deconvolution.

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

if ~exist('ind', 'var') || isempty(ind)
    % display all neurons if ind is not specified.
    ind = 1:size(obj.A, 2);
elseif ind==-1
    ind = size(obj.A,2):-1:1;
end
if ~exist('C2', 'var'); C2=[]; end

if exist('folder_nm', 'var')&&(~isempty(folder_nm))
    % create a folder to save images
    save_img = true;
    cur_cd = cd();
    if ~exist(folder_nm, 'dir'); mkdir(folder_nm);
    else
        fprintf('The folder has been created and old results will be overwritten. \n');
    end
    cd(folder_nm);
else
    save_img = false;
end

% obj.delete(sum(obj.A>0, 1)<max(obj.options.min_pixel, 1));

Amask = (obj.A~=0);
ind_trim = false(size(ind));    % indicator of trimming neurons
ind_del = false(size(ind));     % indicator of deleting neurons
ctr = obj.estCenter();      %neuron's center
gSiz = obj.options.gSiz;        % maximum size of a neuron

% time
T = size(obj.C, 2);
t = 1:T;
if ~isnan(obj.Fs)
    t = t/obj.Fs;
    str_xlabel = 'Time (Sec.)';
else
    str_xlabel = 'Frame';
end

%% keep the log information
if ~save_img
    try
        log_file =  obj.P.log_file;
        flog = fopen(log_file, 'a');
        log_data = matfile(obj.P.log_data, 'Writable', true); %#ok<NASGU>
        manual_intervention.before = obj.obj2struct();
        
        fprintf(flog, '[%s]\b', get_minute());
        fprintf(flog, 'Start manual interventions:\n');
    end
end

% Exclude neurons from the analysis that are lying within close proximity
% of the border and only display the ones that would be taken
ind_f = ind;
ind = ind(~del);
ind_del(del) = 1;

g = figure;
imagesc(Cn); 
colormap gray;hold on;
for mm = 1:numel(ind)
    cont = obj.Coor{ind(mm)}; 
        if size(cont,2) > 1
            plot(cont(1,1:end),cont(2,1:end),'Color', 'b', 'linewidth', 1.5); 
        end

end
  
F = getframe ;
[X,~] = frame2im(F);
All_comp = imresize(X,[size(Cn)], 'bicubic');
close(g)


%% start viewing neurons
figure('position', [1300, 100, 1024, 1024]);
m=1;
while and(m>=1, m<=length(ind))
    
    %% full-frame view
    subplot(241); cla;
    obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    axis equal; axis off;

    %% zoomed-in view
    subplot(242); cla;
    obj.image(obj.A(:, ind(m)).*Amask(:, ind(m))); %
    %     imagesc(reshape(obj.A(:, ind(m)).*Amask(:,ind(m))), obj.options.d1, obj.options.d2));
     axis off;axis equal;

    x0 = ctr(ind(m), 2);
    y0 = ctr(ind(m), 1);
    
    if ~isnan(x0)
        xlim(x0+[-gSiz, gSiz]*2);
        ylim(y0+[-gSiz, gSiz]*2);
    end
    
    if ind_del(ind(m))
        title(sprintf('Neuron %d out of %d', m, length(ind)), 'color', 'r');
    else
        title(sprintf('Neuron %d out of %d', m, length(ind)));
    end
    
    %% Position in Image
    subplot(243); cla;
    cont = obj.Coor{ind(m)}; 
    %figure('position', [1800, 800, 512, 512]);
    imagesc(Cn); colormap gray;hold on;
    axis equal; axis off;
    if size(cont,2) > 1
        plot(cont(1,1:end),cont(2,1:end),'Color', 'r', 'linewidth', 1.5); 
    end
    
    %% Position in Image next to all others
    subplot(244); cla;
    %figure('position', [1800, 800, 512, 512]);
    imagesc(All_comp); hold on;
    
    cont = obj.Coor{ind(m)};
    axis equal; axis off;
    if size(cont,2) > 1
        plot(cont(1,1:end),cont(2,1:end),'Color', 'r', 'linewidth', 1.5); 
    end
    hold off
    
    %% temporal components
    subplot(2,4,5:8);cla;
    if ~isempty(C2)
        plot(t, C2(ind(m), :)*max(obj.A(:, ind(m))), 'k', 'linewidth', 2); hold on;
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))), 'r');
        plot(t, repelem(mean(C2(ind(m), :)) + 3*std(C2(ind(m), :)), size(t, 2))*max(obj.A(:, ind(m))), '--b', 'linewidth', 2);
    else
        
        plot(t, obj.C(ind(m), :)*max(obj.A(:, ind(m))));
    end
    xlim([t(1), t(end)]);
    xlabel(str_xlabel);
    
    %% save images
    if save_img
        drawnow();
        saveas(gcf, sprintf('neuron_%d.png', ind(m)));
        m = m+1;
    else
        fprintf('Neuron %d, keep(k, default)/delete(d)/delete all(da)/backward(b)/end(e):    ', m);
        
        temp = input('', 's');
        if temp=='d'
            ind_del(ind(m)) = true;
            m = m+1;
        elseif strcmpi(temp, 'b')
            m = m-1;
        elseif strcmpi(temp, 'da')
            ind_del(ind(m):end) = true;
            break;
        elseif strcmpi(temp, 'k')
            ind_del(ind(m)) = false;
            m= m+1;
        elseif strcmpi(temp, 'e')
            break;
        elseif ~isnan(str2double(temp))
            m = m + floor(str2double(temp));
            m = max(m, 1);
            m = min(m, length(ind));
            fprintf('jump to neuron %d / %d\n', m, length(ind));
        else
            m = m+1;
        end
    end
    
end



if save_img
    cd(cur_cd);
else
%     if ~isempty(ind(ind_trim))
%         obj.A(:, ind(ind_trim)) = obj.A(:,ind(ind_trim)).*Amask(:, ind(ind_trim));
%         try
%             fprintf(flog, '\n\tFollowing neurons were trimmed:\n');
%             ids_trimmed = ind(ind_trim);
%             for m=1:length(ids_trimmed)
%                 fprintf(flog, '%2d, ', ids_trimmed(m));
%             end
%             fprintf(flog, '\n');
%         end
%     end
    
    if ~isempty(ind_f(ind_del))
        try
            fprintf(flog, '\tDeleting manually selected neurons:\n');
        end
        obj.delete(ind_f(ind_del));
    end
    %     obj.Coor = obj.get_contours(0.9);
    
    
    return;
end
try
    fprintf(flog, '[%s]\b', get_minute());
    fprintf(flog, 'Finished the manual intervention.\n');
    fprintf(flog, '[%s]\b', get_minute());
    if obj.options.save_intermediate
        manual_intervention.after = obj.obj2struct(); %#ok<STRNU>
        tmp_str = get_date();
        tmp_str=strrep(tmp_str, '-', '_');
        eval(sprintf('log_data.manual_%s = manual_intervention;', tmp_str));
        
        fprintf(flog, '\tThe results were saved as intermediate_results.manual%s\n\n', tmp_str);
    end
    fclose(flog);
end
end

