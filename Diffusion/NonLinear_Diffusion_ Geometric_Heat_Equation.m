% Load the original image
I = imread('wheatcypress.png');

% Convert channels to double
Ir = im2double(I(:,:,1));
Ig = im2double(I(:,:,2));
Ib = im2double(I(:,:,3));

% time step parameter
delta_t = .025; 

% Create arrays for the current and next iterations
% we only need two arrays i.e. current and previous for each channel
current_R = Ir;
current_G = Ig;
current_B = Ib;

for n = 1:50 % Number of iterations 
    % Calculate the gradients for each channel
    [Ix_R, Iy_R] = gradient(current_R);
    [Ix_G, Iy_G] = gradient(current_G);
    [Ix_B, Iy_B] = gradient(current_B);
    
    % Calculate the magnitude of the gradients with epsilon added for stability
    gradient_magnitude_R = sqrt(Ix_R.^2 + Iy_R.^2) + eps;
    gradient_magnitude_G = sqrt(Ix_G.^2 + Iy_G.^2) + eps;
    gradient_magnitude_B = sqrt(Ix_B.^2 + Iy_B.^2) + eps;
    
    
    % Calculate the divergence of the normalized gradient
    % curvature vector (divergence of normalized graident )
    div_gradient_R = divergence(Ix_R ./ gradient_magnitude_R, Iy_R ./ gradient_magnitude_R);
    div_gradient_G = divergence(Ix_G ./ gradient_magnitude_G, Iy_G ./ gradient_magnitude_G);
    div_gradient_B = divergence(Ix_B ./ gradient_magnitude_B, Iy_B ./ gradient_magnitude_B);
    
    % Update the next iteration's channels using the Geometric Heat Equation
    current_R = current_R + (delta_t * div_gradient_R);
    current_G = current_G + (delta_t * div_gradient_G);
    current_B = current_B + (delta_t * div_gradient_B);
    
end


% Combine the results to create a diffused color image
Idiffused = cat(3, current_R, current_G, current_B);

% Display the diffused color image
imshow(Idiffused);
title('Geometric heat equation Result n = 50');


% Step 4: Display the results and overlay a specific isophote
figure;
subplot(131);
imshow(current_R);
title('Red Channel');

subplot(132);
imshow(current_G);
title('Green Channel');

subplot(133);
imshow(current_B);
title('Blue Channel');


% Overlay a specific isophote (e.g., contour at intensity level 0.5)
figure;

% Subplot 1: Red Channel Contour
subplot(1, 3, 1);
contour(current_R, [0.5, 0.5], 'r');
title('Red Channel Isophote');

% Subplot 2: Green Channel Contour
subplot(1, 3, 2);
contour(current_G, [0.5, 0.5], 'g');
title('Green Channel Isophote');

% Subplot 3: Blue Channel Contour
subplot(1, 3, 3);
contour(current_B, [0.5, 0.5], 'b');
title('Blue Channel Isophote');

sgtitle('Isophote Contours for RGB Channels');


