% define/ manually selected 3D points on the basis of size of the cube given 
% i.e. each tile is 18 mm X 18 mm
points_3D = [
    0, 0, 0;
    36, 18, 0;
    36, 36, 0;
    0, 0, 36;
    0, 18, 18;
    0, 18, 36;
    18, 0, 36;
    18, 0,18;
    36, 0, 18;
    0, 36, 0;
    ];


% manually selected/ hand picked 2D points correspodning to above 3D points
points_2D = [
    3232, 1325;
    3624, 1303;
    3551, 1560;
    2635, 959;
    2866, 1420;
    2583, 1244;
    2888, 808;
    3173, 988;
    3408, 834;
    3093, 1893;
    ];

%  Calculate the P matrix
P = calculate_camera_matrix(points_3D, points_2D);

% display the calculated finite projective camera matrix
disp('The finite projective camera matrix P is:');
format long g;
disp(P);


% Testing the accuracy
% 3D points selected for testing the model
test_points_3D = [
   18, 36, 0;
   18, 0, 0;
   0, 36, 18;
   0, 36, 36;
   36, 0, 36;
];

% manually selected/ hand picked 2D points correspodning to above 3D points
true_points_2D = [
    3320, 1739;
    3470, 1171;
    2822, 1710;
    2547, 1512;
    3111, 665;
];


% Projecting the TEST 3D points onto 2D using the projection matrix P
% add 1 as the 4th Value for homogeneous matrix calculations
homogeneous_points_3D = [test_points_3D, ones(size(test_points_3D, 1), 1)];
projected_points_2D_homogeneous =  homogeneous_points_3D * P';
% get the coordinates x and y by dividing xw and yw with w
projected_points_2D = projected_points_2D_homogeneous(:, 1:2) ./ projected_points_2D_homogeneous(:, 3);

% rounding the numbers to 2 decimal places
projected_points_2D_rounded = round(projected_points_2D * 100) / 100;

disp("True 2D Points");
disp(true_points_2D);

% Display the original and the scaled matrices
disp('Predicted 2D points');
format long g;
disp(projected_points_2D_rounded);

%  calculate the positional error between the projected points and
%  true 2D points
positional_error = true_points_2D - projected_points_2D_rounded;

% mean of the positional error
mean_error = mean(positional_error);
% standard deviation of the positional error
std_error = std(positional_error);

% Display the mean and standard deviation of the error
disp('Mean positional error is:');
disp(mean_error);
disp('Standard deviation of the positional error is:');
disp(std_error);

function P = calculate_camera_matrix(points_3d, points_2d)
    % Data Normalization for the 3D points
    [mean_3D, std_3D] = calculate_3d_mean_std(points_3d);

    points_3d_normalized = (points_3d - mean_3D) ./ std_3D;

    M2 = [1/std_3D, 0, 0, -mean_3D(1)/std_3D;
          0, 1/std_3D, 0, -mean_3D(2)/std_3D;
          0, 0, 1/std_3D, -mean_3D(3)/std_3D;
          0, 0, 0, 1];

    
    % Data Normalization for the 2D points
    [mean_2D, std_2D] = calculate_2d_mean_std(points_2d);

    points_2d_normalized = (points_2d - mean_2D) ./ std_2D;

    M1 = [ 1/std_2D, 0, -mean_2D(1)/std_2D;
           0, 1/std_2D, -mean_2D(2)/std_2D;
           0, 0, 1];
    
    
    % matrix A of size 2N X 12 because each data points contributes 
    % into 2 rows
    N = size(points_2d, 1);
    A = zeros(2*N, 12);
    for i = 1:N
        X = points_3d_normalized(i, 1);
        Y = points_3d_normalized(i, 2);
        Z = points_3d_normalized(i, 3);
        x = points_2d_normalized(i, 1);
        y = points_2d_normalized(i, 2);

        A(2*i-1,:) = [X, Y, Z, 1, 0, 0, 0, 0, -x*X, -x*Y, -x*Z, -x];
        A(2*i,:)   = [0, 0, 0, 0, X, Y, Z, 1, -y*X, -y*Y, -y*Z, -y];

    end
    
    % using SVD to get eigne vector with smallest eigen value i.e. V column
    [~, ~, V] = svd(A);

    % V is the eigenvector corresponding to the smallest eigenvalue
    % P normalsied is eqivalent to V
    P_normalized = reshape(V(:, end), 4, 3)';
    
    % get the original P matrix
    P = inv(M1) * P_normalized * M2;
end


% calculate the mean and standard deviation for 2D points
function [mean_values, std_value] = calculate_2d_mean_std(points2D)

    % total number of 2D points
    n = size(points2D, 1);

    % mean calculation
    sumx = sum(points2D(:,1));
    sumy = sum(points2D(:,2));
    meanx = sumx/n;
    meany = sumy/n;
    mean_values = [meanx, meany];
    
    % standard deviation
    varx = sum((points2D(:, 1) - meanx).^2);
    vary = sum((points2D(:, 2) - meany).^2);
    std_value = sqrt((varx + vary)/(2*n));
end

% calculate the mean and standard deviation for 3D points
function [mean_values, std_value] = calculate_3d_mean_std(points3D)
    
    % total number of 3D points
    n = size(points3D, 1); 

    % mean
    sumX = sum(points3D(:,1));
    sumY = sum(points3D(:,2));
    sumZ = sum(points3D(:,3));
    
    meanX = sumX/n;
    meanY = sumY/n;
    meanZ = sumZ/n;
    mean_values = [meanX, meanY, meanZ];
    
    % standard deviation
    varX = sum((points3D(:, 1) - meanX).^2);
    varY = sum((points3D(:, 2) - meanY).^2);
    varZ = sum((points3D(:, 3) - meanZ).^2);
    std_value = sqrt((varX + varY + varZ) / (3*n));
end

