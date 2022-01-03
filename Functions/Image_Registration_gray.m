%% Image Registration

function [Registered_Images, Shifts] = Image_Registration_gray(Data, Template, Precision, psf, varargin)
    
    % Extract the indices for Motion correction 
    if size(varargin, 2) == 1
        X_Dim = varargin{1}; Y_Dim = varargin{1};
    elseif size(varargin, 2) == 2
        X_Dim = varargin{1}; Y_Dim = varargin{2};
    end
    
    % Filter Data for registration
    Data_y = imfilter(Data, psf, 'replicate');
    Template_y = imfilter(Template, psf, 'replicate');

    % Subselect region for Motion correction
    Moving_sm = Data_y(Y_Dim, X_Dim, :);
    fixed_smoothed = Template_y(Y_Dim, X_Dim);

    % Preallocate Space
    Shifts = nan(size(Moving_sm, 3), 2);
    Registered_Images = nan(size(Data, 1), size(Data, 2), size(Data,3)); %Preallocate space  

    parfor jj = 1:size(Moving_sm, 3) %run through the frames
        DFToutput = dftregistration(fft2(fixed_smoothed), fft2(Moving_sm(:,:,jj)), Precision); % Calculate shift on contrast enhanced images
        Registered_Images(:,:,jj) = imtranslate_old(Data(:, :, jj), [DFToutput(3) DFToutput(4)]); %Save
        Shifts(jj, :) = [DFToutput(3) DFToutput(4)];
    end
      
end
