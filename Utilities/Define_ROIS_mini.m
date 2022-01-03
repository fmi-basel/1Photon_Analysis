%% Define ROIs for Motion Registration
% Author Julian Hinz, Luthi lab 2018 (julian.hinz@fmi.ch)
% based on code of the Function "Crop_it"
function [Position_all] = Define_ROIS_mini(I, num)


gg = figure('position', [500, 200, 1000, 800]); 
imshow(I, [0 1.8]); box off; set(gca, 'Xtick',[], 'Ytick', [])

Position = struct;

if isempty(num)
    num = 0;
    p(1) = 1;
    h = imrect;
    while p(1) > 0
      p = wait(h);          
      if ~isempty(p) && p(1) > 0                                    
          rectangle('Position', p, 'LineWidth', 1, 'EdgeColor','r');
          num = num + 1;
          Position.(['G_' num2str(num)]) = p;
      end
    end    
else
    h = imrect;
    for i = 1:num
        p = wait(h);  
        rectangle('Position', p, 'LineWidth', 1, 'EdgeColor','r');
        Position.(['G_' num2str(i)]) = p;       
    end
end

close(gg)

%%% Get Positions

Position_all = cell(num, 2);
% Indices_include = cell(num, 1);
% Crop_String = cell(num, 1);

for i = 1:num
    Temp = round(Position.(['G_' num2str(i)]));
    Position_all{i, 1} = Temp(1):Temp(3)+Temp(1);
    Position_all{i, 2} = Temp(2):Temp(4)+Temp(2);
    
%     bw = poly2mask([Temp(1); Temp(3); Temp(3); Temp(1)], [Temp(2);  Temp(2); Temp(4); Temp(4)], size(I, 1), size(I, 2));
%     Indices_include{i} = find(bw == 1);
%     [Indices_x, Indices_y] = ind2sub([size(bw)], Indices_include{i});
%     try
%     Crop_String{i} = [Indices_y(1), Indices_x(1), Indices_y(end) - Indices_y(1), Indices_x(end)-Indices_x(1)];    
%     catch
%         disp('Stop')
%     end
end

end
