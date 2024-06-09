% load the input image
Iorig = imread('medial1.png');

% get one channel to work as a grayscale image
I  = im2double(Iorig(:,:,2));

% number of levels in the Gaussian pyramid
num_levels = 6;

% a composite image to display the gaussian pyramid
composite_image_output = ones(1150, 512);

current_y = 1130;

input_image = I;

for level = 1:num_levels
    % sigma value for the given octave
    sigma = 2^(level - 1);
    
    % compute the 2D Gaussian of the image 
    output_image = imgaussfilt(input_image, sigma);

    [h, w] = size(output_image);
    
    % display image in the composite figure
    composite_image_output(current_y - h +1 :current_y, 1:w) = output_image;

    % update the x-position for the next level
    current_y = current_y - h - 10;

    % downsample the input image for the next octave level
    input_image = imresize(output_image, 1/2, 'bilinear');
end

% Display the composite figure
figure('Name','6 Octave Levels in the Gaussian Pyramid','NumberTitle','off');
imshow(composite_image_output);
title('6 Octave Levels in the Gaussian Pyramid');
