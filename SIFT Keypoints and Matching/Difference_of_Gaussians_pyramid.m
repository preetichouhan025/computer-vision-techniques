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

%  composite image to display the DOG (Difference of Gaussian) pyramid
composite_image_dog = ones(1150, 512*4+100);

current_y = 1130;

% for each octave level
for octave = 1:num_octaves
    current_x = 20;
    % for ecah scale at the given octave level
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
        [h, w] = size(difference_of_gaussian);

        % Prepare the image for display in the composite figure
        % We are skipping  the top most box - therefore we end up with 23
        % DOG
        
        if (octave == 6 ) && (scale ==4 )
            continue;
        else
            composite_image_dog(current_y - h + 1:current_y, current_x:current_x + w - 1) = difference_of_gaussian;   
        end

        current_x = current_x + 512 + 20;
    end

    % Update the y-position for the next scale
    current_y = current_y - h - 10;
end

% Display the DOG pyramid
figure('Name','Difference of Gaussians (DOG) Pyramid with 6 ocatves and 4 scales on each octave','NumberTitle','off');
imshow(composite_image_dog);
title('Difference of Gaussians (DOG) Pyramid with 6 ocatves and 4 scales on each octave');
