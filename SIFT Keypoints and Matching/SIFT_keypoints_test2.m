% load the input image 1
Iorig1 = imread('medial1.png');
% one channel to work as a grayscale image
Iorig1  = im2double(Iorig1(:,:,2));

% load the input image 2
Iorig2 = imread('medial2.png');
% one channel to work as a grayscale image
Iorig2  = im2double(Iorig2(:,:,2));

% detect the SHIFT keypoints (x, y, sigma, theta orientation, fetaure vector) 
first_image_KeyPoints = detectKeypoints(Iorig1);
second_image_KeyPoints = detectKeypoints(Iorig2);

% Initialize arrays to store the matched or the consistent keypoints
% between the two same images with different orientation
matched_keypoint_image1 = [];
matched_keypoint_image2 = [];

% a threshold for matching 
% this threshold decides how much teh keypoints feature vector needs to be
% matching to be considered consistent
match_threshold = 1;

% for all teh keypoints in teh image 1
for i = 1:size(first_image_KeyPoints, 1)
    feature_vector1 = first_image_KeyPoints(i, 4:end);
    
    % Initialize variables to keep track of the best match
    best_match_index = -1;
    best_match_value = -1;
     
    % for all the keypoints in the second image
    for j = 1:size(second_image_KeyPoints, 1)
        feature_vector2 = second_image_KeyPoints(j, 4:end);
        
        % use the cross-correlation to measure similarity between the
        % feature vectors of the two images
        correlation = xcorr(feature_vector1, feature_vector2);
        
        % the max correlation value
        max_correlation = max(correlation);
        % to check if we need to keep the keypoint or not
        if max_correlation > best_match_value
            best_match_index = j;
            best_match_value = max_correlation;
        end
    end
    
    % check if the best match is greater than the threshold
    % if yes then store the matched keypoint
    if best_match_value > match_threshold
        % Store the matched keypoints and their positions along with
        % the scale
        matched_keypoint_image1 = [matched_keypoint_image1; first_image_KeyPoints(i, 1:3)];
        matched_keypoint_image2 = [matched_keypoint_image2; second_image_KeyPoints(best_match_index, 1:3)];
    end
end

figure('Name','Keypoints Detected using Difference of Gaussian Pyramid','NumberTitle','off');
% display medial1.png
imshow(Iorig1);
title('Consistent Keypoints - medial1.png');
hold on;
plot(matched_keypoint_image1(:, 1), matched_keypoint_image1(:, 2), 'ro', 'MarkerSize', 10);

% display medial2.png
figure('Name','Keypoints Detected using Difference of Gaussian Pyramid','NumberTitle','off');
imshow(Iorig2);
title('Consistent Keypoints - medial2.png');
hold on;
plot(matched_keypoint_image2(:, 1), matched_keypoint_image2(:, 2), 'ro', 'MarkerSize', 10);


function [keypoints, DOG_gaussian_pyramid] = detectKeypoints(I)
    
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
    
    %  Create an array to store the Difference of Gaussian images
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
            difference_of_gaussian = next_gaussian - current_gaussian;
    
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

    num_bins = 16; % number of bins in the histogram
    
    keypoints = []; % Initialize an empty array to store the SHIFT keypoints
    
    for octave = 1:num_octaves
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
    
                        % Specify a minimum contrast threshold
                        min_contrast_threshold = 0.18; % Adjust this value as needed
                        
                        % check if the keypoint's contrast is above the
                        % threshold or not to decide whether to keep it or skip
                        if (is_max || is_min) && contrast > min_contrast_threshold
                            
                            x = j;
                            y = i;
                            % compute the sigma value from the scale and octave
                            scale_sigma = 2^((octave - 1) + (scale-1)/m);

                            % get a m x m window centered at the keypoint
                            window = current_DoG(i-7:i+8, j-7:j+8);
                         
                            % orientation histogram
                            histogram = histcounts(window(:), num_bins);
                            
                            
                            % fidn the peak 
                            [peaks, ~] = findpeaks(histogram);
                            % location of the peak
                            [~, max_peak_idx] = max(peaks);
                                   
                            % Shift the histogram so that the peak of the histogram aligns with the origin
                            feature_vector = circshift(histogram, [0, num_bins / 2 - max_peak_idx]);
                            
                            % final SHIFT KEYPOINT
                            feature_vector = [x, y, scale_sigma, feature_vector];
                            keypoints = [keypoints; feature_vector];
                        end
                    end
                end

            end
        end
    end
end    