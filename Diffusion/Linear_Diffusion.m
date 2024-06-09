% load the input Original image
Iorig = imread('wheatcypress.png');

% Convert channels to double
Ir = im2double(Iorig(:,:,1));
Ig = im2double(Iorig(:,:,2));
Ib = im2double(Iorig(:,:,3));

% sigma (defines the smoothing levels)
sigmas = [2, 4, 8, 12, 16];

% Δt - time step parameter
delta_t = .25;

% Initialize two array to store the results of blurred resulst and the
% smoothened results
% defining the space for 5 sigma values to store three results from three
% channels
blurred_output = cell(5, 3);
smoothed_output = cell(5,3);

% Initialize two arrays to store the final diffused images for 5 sigma
% values
diffused_images_blurred = cell(5, 1);
diffused_images_smoothed = cell(5, 1);

% Loop through each color channel and each value of sigma
for s = 1:length(sigmas)
    
    sigma = sigmas(s);
    
    % Blurring the image channels using imgaussfilt i.e. 2D Gaussian
    blurred_R = imgaussfilt(Ir, sigma);
    blurred_G = imgaussfilt(Ig, sigma);
    blurred_B = imgaussfilt(Ib, sigma);

    % Store the results
    blurred_output{s, 1} = blurred_R;
    blurred_output{s, 2} = blurred_G;
    blurred_output{s, 3} = blurred_B;

    % Store the diffused image for this sigma
    diffused_images_blurred{s} = cat(3, blurred_R, blurred_G, blurred_B);   
    
    % Smoothing with  heat equation
    smoothed_R = Ir; % Initialize the smoothed image
    smoothed_G = Ig;
    smoothed_B = Ib;

    % total difussion time
    n = round((0.5*(sigma^2))/ delta_t);

    for i = 1:n
        % Compute gradients using the gradient function
        [Ix_R, Iy_R] = gradient(smoothed_R);
        [Ix_G, Iy_G] = gradient(smoothed_G);
        [Ix_B, Iy_B] = gradient(smoothed_B);

        [Ir_dx, ~] = gradient(Ix_R);
        [~, Ir_dy] = gradient(Iy_R);

        [Ig_dx, ~] = gradient(Ix_G);
        [~, Ig_dy] = gradient(Iy_G);

        [Ib_dx, ~] = gradient(Ix_B);
        [~, Ib_dy] = gradient(Iy_B);
        
        % Update I using the discretized heat equation
        smoothed_R = (delta_t * (Ir_dx + Ir_dy)) + smoothed_R;
        smoothed_G = (delta_t * (Ig_dx + Ig_dy)) + smoothed_G;
        smoothed_B = (delta_t * (Ib_dx + Ib_dy)) + smoothed_B;
    end

    smoothed_output{s, 1} = smoothed_R;
    smoothed_output{s, 2} = smoothed_G;
    smoothed_output{s, 3} = smoothed_B;

    % Store the diffused image for this sigma
    diffused_images_smoothed{s} = cat(3, smoothed_R, smoothed_G, smoothed_B);
end

% Compare the results and generate surface plots
% for each sigma 
for s = 1:length(sigmas)
    sigma = sigmas(s);
    % for each channel
    for c = 1:3
        % Calcuate the absolute difference between the two results
        diff_image = abs(smoothed_output{s, c} - blurred_output{s, c});

        % Display the difference
        figure;
        surf(diff_image);
        colormap copper;
        shading interp;
        title(['Sigma = ', num2str(sigma), ', Channel ', num2str(c)]);
        xlabel('X-axis');
        ylabel('Y-axis');
        zlabel('Absolute Difference between 2D Gaussian and Heat Equation');
    end
end


for i = 1:length(sigmas)
    % Display the images generated after the 2D Gussian blurring and Heat
    % equtaion smoothing
    % Display the diffused color images side by side for each sigma
    subplot(2, 5, i);
    imshow(diffused_images_smoothed{i});
    title(['Heat Equation , Sigma = ', num2str(sigmas(i))]);
    
    subplot(2, 5, i + 5);
    imshow(diffused_images_blurred{i});
    title(['2D Gaussian, Sigma = ', num2str(sigmas(i))]);
end

