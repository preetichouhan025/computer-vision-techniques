% load the original input image
I = imread('wheatcypress.png');
I = im2double(I); % Convert the image to a double for numerical computations

% set the parameters values
delta_t = 0.1; % time step parameter
K_values = [0.005, 0.01, 0.05, 0.1];
itertaions = 100; % number of iterations


for k = 1:length(K_values)
    K = K_values(k);

    % create arrays for the current and next iterations
    % we only need two arrays i.e. current and previous for each channel

    current_R = I(:, :, 1);
    current_G = I(:, :, 2);
    current_B = I(:, :, 3);

    next_R = zeros(size(current_R));
    next_G = zeros(size(current_G));
    next_B = zeros(size(current_B));

    % total number of iterations
    for n = 1:itertaions

        % compute the gradient of the current image
        [Ix_R, Iy_R] = gradient(current_R);
        [Ix_G, Iy_G] = gradient(current_G);
        [Ix_B, Iy_B] = gradient(current_B);

        gradient_magnitude_R = sqrt(Ix_R.^2 + Iy_R.^2);
        gradient_magnitude_G = sqrt(Ix_G.^2 + Iy_G.^2);
        gradient_magnitude_B = sqrt(Ix_B.^2 + Iy_B.^2);
        
        % compute the conduction coefficient c for each channel
        C_R = 1 ./ (1 + ((gradient_magnitude_R.^2) / K));
        C_G = 1 ./ (1 + ((gradient_magnitude_G.^2) / K));
        C_B = 1 ./ (1 + ((gradient_magnitude_B.^2) / K));

        % Calculate the Laplacian of the current image
        % Laplacian1 = del2(current_I);
        Laplacian_R = divergence(Ix_R .* C_R, Iy_R .* C_R);
        Laplacian_G = divergence(Ix_G .* C_G, Iy_G .* C_G);
        Laplacian_B = divergence(Ix_B .* C_B, Iy_B .* C_B);

        % compute the next iteration's image using the Perona-Malik equation
        next_R = current_R + (delta_t * Laplacian_R);
        next_G = current_G + (delta_t * Laplacian_G);
        next_B = current_B + (delta_t * Laplacian_B);

        % Swap the current and next iterations
        current_R = next_R;
        current_G = next_G;
        current_B = next_B;

    end

    % Display the result for the current K value
    output_image = cat(3, current_R, current_G, current_B);
    subplot(1, length(K_values), k);
    imshow(output_image);
    title(['K = ', num2str(K)]);

end

