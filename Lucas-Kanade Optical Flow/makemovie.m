addpath("./");
folder_name = "Basketball";
smoothing_factor = 0.2;
start_index = 7;
end_index = 14;
% window size for computing the local averaging
window_size = 25;
vidObj = VideoWriter(strcat(folder_name, '_Video_Window_size_50','.avi'));

vidObj.FrameRate = 1.5;
open(vidObj);


% smoothing the frame with respect to the 3rd dimension which is time
smooth_frames(folder_name,smoothing_factor,start_index,end_index);

% Gaussian 2D smoothing
for i = start_index: end_index
    input_smoothed_image = imread(fullfile(folder_name,strcat('image_smoothed_',num2str(i),'.png')));
    input_gaussian_smoothed_image = imgaussfilt(input_smoothed_image,2);
    imwrite(input_gaussian_smoothed_image,fullfile(folder_name,strcat('image_smoothed_',num2str(i),'.png'))); 
end 

for i=start_index:end_index-1  
    demo_optical_flow(folder_name,i,i+1, window_size); %i+1 for subsequent frames
    fprintf('Frame #: %d\n',i)
    frame = getframe(gcf);
    writeVideo(vidObj,frame);
end

close(vidObj);
