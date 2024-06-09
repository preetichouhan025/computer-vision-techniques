% load the input image
Iorig = imread('medial1.png');

% one channel to work as a grayscale image
I  = im2double(Iorig(:,:,2));

% Number of octaves - Space levels
num_octaves = 6;
% Number of scales within each octave - Scale levels
m = 4;

% an array to store the output for each space and scale levels
gaussian_pyramid = cell(num_octaves, m);

input_image = I;

% for each octave level
for octave = 1:num_octaves

    % for ecah scale level
    for scale = 1:m
        
        % compute the Sigma for this scale level based on the scale and the
        % octave 
        sigma = 2^((octave - 1) + (scale-1)/m);

        % compute the 2D Gaussian of the image
        gaussian = imgaussfilt(input_image, sigma);

        % Store the results for this scale and ocatve in the gaussian
        % pyramid array
        gaussian_pyramid{octave, scale} = gaussian;
    end

    input_image = imresize(gaussian, 0.5, 'bilinear');
end

% Create an array to store the Difference of Gaussian images
DOG_gaussian_pyramid = cell(num_octaves, m);

% for each octave level
for octave = 1:num_octaves

    % for each scale at the given octave level
    for scale = 1:m

        % get the 2D Gaussian of the image for given octave and scale
        current_gaussian = gaussian_pyramid{octave, scale};


        % if the scale is less than 4 then its fine to get the next local
        % gradients from the same octave level and no resizing is required
        if scale < 4
            next_gaussian = gaussian_pyramid{octave, scale+1};

        % if the scale level is 4 then the next gaussian is from scale 1 
        % of the next octave. Also a upsampling is required here so taht
        % the differences can be calculated correctly
        elseif scale == 4
            if octave < 6
                next_gaussian = gaussian_pyramid{octave+1, 1};
                next_gaussian = imresize(next_gaussian, size(current_gaussian), 'bilinear');
            end
        end 

        % compute the DoG for this scale level 
        difference_of_gaussian = (scale^2) * (next_gaussian - current_gaussian);

        % Prepare the image for display in the composite figure
        % We are skipping  the top most box - therefore we end up with 23
        % DOG
        
        if (octave == 6 ) && (scale ==4 )
            continue;
        else
            DOG_gaussian_pyramid{octave, scale}  = difference_of_gaussian;   
        end
    end
end

%  minimum contrast threshold to keep or remove the keypoints
min_contrast_threshold = 1;
% initialize an empty array to store the SIFT keypoints (𝑥, 𝑦, 𝜎)
keypoints = []; 
 % for each octave
for octave = 1:num_octaves
    % for each scale
    for scale = 1:m

        % skip calculating for the top most gaussian and the first gaussian
        if (octave == 6 ) && (scale >= 3 )
            continue;
        elseif (octave ==1) && (scale ==1)
            continue;
        else
            % get the current DOG
            current_DoG = DOG_gaussian_pyramid{octave, scale};
            
            % if the scale is less then 4 then the next DOF is from the
            % same octave level and no reszing is required
            if scale < 4
                above_DoG = DOG_gaussian_pyramid{octave, scale + 1};

            % if the scale is 4 then teh next DOG is from next octave with scale 1
            % here perform the upsampling for better neighbourhood
            % comparsion
            elseif scale == 4
                above_DoG = DOG_gaussian_pyramid{octave+1,1};
                above_DoG = imresize(above_DoG, size(current_DoG), 'bilinear');
            end

            % is scale is greater than 1 then the previous DOG is from the
            % same octave level and no resizing is required
            if scale > 1
                below_DoG = DOG_gaussian_pyramid{octave, scale - 1};

            % if the scale is 1 then the below DOG is from the previous
            % octave with scale 4
            % we need to perform the downsampling to better compare the
            % neighbour pixels
            elseif scale ==1
                below_DoG = DOG_gaussian_pyramid{octave-1, 4};
                below_DoG = imresize(below_DoG, size(current_DoG), 'bilinear');
            end     
    
            [rows, cols] = size(current_DoG);

            % neighbourhood comparsions
            for i = 2:rows-2
                for j = 2:cols-2

                    % for ecah pixel value as the keypoint
                    pixel_value = current_DoG(i, j);
    
                    % get the pixel neighbors from above nad below levels
                    neighbors = [current_DoG(i-1:i+1, j-1:j+1), above_DoG(i-1:i+1, j-1:j+1), below_DoG(i-1:i+1, j-1:j+1)];
                    
                    % compute the pixel intensity differences
                    is_max = any(pixel_value == max(neighbors(:)));
                    is_min = any(pixel_value == min(neighbors(:)));
                    
                    % keypoint contrast
                    contrast = abs(max(neighbors(:)) - min(neighbors(:)));
                    
                    % check if the keypoint's contrast is above the
                    % threshold or not to decide whether to keep it or skip
                    if (is_max || is_min) && contrast > min_contrast_threshold
                       
                        x = j;
                        y = i;
                        % compute the sigma value from the scale and octave
                        sigma_value = 2^((octave - 1) + (scale-1)/m);
                        keypoints = [keypoints; [x, y, sigma_value]];
                    end
                end
            end
        end
    end
end


% Visualize filtered keypoints on the original image
figure('Name','Keypoints Detected using Difference of Gaussian Pyramid','NumberTitle','off');
imshow(Iorig);
title('Keypoints Detected usng Difference of Gaussian Pyramid - medial1.png');
hold on;

for i = 1:size(keypoints, 1)
    x = keypoints(i, 1);
    y = keypoints(i, 2);
    sigma = keypoints(i, 3);

    if sigma < 2
        color = 'r';
    elseif sigma < 3
        color = 'g';
    else
        color = 'b';
    end
    
    % draw a circle at each keypoint location
    viscircles([x, y], 1.5 * sigma, 'EdgeColor', color, 'LineWidth', 1);
end
