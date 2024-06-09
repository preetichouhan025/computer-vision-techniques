addpath("./");

folder_name = "Backyard"; % or Basketball
smoothing_factor = 0.2;
start_index = 7;
end_index = 14;
% window size for computing the local averaging
window_size = 25;

% smoothing the frame with respect to the 3rd dimension which is time
smooth_frames(folder_name,smoothing_factor,start_index,end_index);

% Gaussian 2D smoothing
for i = start_index: end_index
    input_smoothed_image = imread(fullfile(folder_name,strcat('image_smoothed_',num2str(i),'.png')));
    input_gaussian_smoothed_image = imgaussfilt(input_smoothed_image,2);
    imwrite(input_gaussian_smoothed_image,fullfile(folder_name,strcat('image_smoothed_',num2str(i),'.png')));
    
end    

% compute motion field estimates overlayed on original images
for i = start_index:end_index-1
    demo_optical_flow(folder_name,i,i+1, window_size);
end
