% Iorig = imshow("images\IMG_2216.jpg");
% imshow(Iorig);

% the camera calibration matrix K- intrinsic parameter already given
K = [3063.73307388903, 0, 2023.71937454517;
    0, 3050.01366910237, 1499.05930779402;
    0, 0, 1];


%==============================  IMAGE 1 ===============================
 
fprintf("\n--------------------- IMAGE 1 Data ----------------------------\n");
% define manually selected 3D points on the basis of size of the cube given 
% i.e. each tile is 18 mm X 18 mm
points_3D_image1 = [ 36, 18, 0; 36, 36, 0; 18, 0, 0; 18, 18, 0; 36, 0, 0; 36, 0, 18;];


% manually selected/ hand picked 2D points correspodning to above 3D points
points_2D_image1 = [3624, 1303; 3551, 1560; 3470, 1171; 3404, 1461; 3694, 1017; 3408, 834;];


%  Calculate the P matrix
P1 = calculate_camera_matrix(points_3D_image1, points_2D_image1);

% get the Rotation matrix R and Transaltion matrix T using  P & K matrix
[R1, T1] = decompose_camera_matrix(P1, K);

camera_position1 = T1;

% get the camera coordinates
camera_coordinates1 = get_camera_cooridnates(R1, points_3D_image1, T1);

fprintf("The camera coordinates are \n");
disp(camera_coordinates1);


%==============================  IMAGE 2 ===============================
 
fprintf("\n--------------------- IMAGE 2 Data ----------------------------\n");
% define manually selected 3D points on the basis of size of the cube given 
% i.e. each tile is 18 mm X 18 mm
points_3D_image2 = [ 18, 36, 0; 18, 0, 0; 0, 36,18; 0, 36,36; 36, 0, 36; 36, 18, 0;];


% manually selected/ hand picked 2D points correspodning to above 3D points
points_2D_image2 = [ 3177, 1820; 3276, 1230; 2646, 1684; 2459, 1461; 3140, 669; 3569, 1402; ];


%  Calculate the P matrix
P2 = calculate_camera_matrix(points_3D_image2, points_2D_image2);

% get the Rotation matrix R and Transaltion matrix T using  P & K matrix
[R2, T2] = decompose_camera_matrix(P2, K);

camera_position2 = T2;

% get the camera coordinates
camera_coordinates2 = get_camera_cooridnates(R2, points_3D_image2, T2);

fprintf("The camera coordinates are \n");
disp(camera_coordinates2);



%==============================  IMAGE 3 ===============================
 
fprintf("\n--------------------- IMAGE 3 Data ----------------------------\n");
% define manually selected 3D points on the basis of size of the cube given 
% i.e. each tile is 18 mm X 18 mm
points_3D_image3 = [0, 0,0; 18, 18, 0; 36, 36, 0; 0, 18, 0; 36, 0, 36; 
    18, 0, 18; 0, 0, 36;];


% manually selected/ hand picked 2D points correspodning to above 3D points
points_2D_image3 = [2517, 1406; 2873, 1651; 3221, 1875; 2495, 1714; 3019, 812;
    2767, 1098; 2290, 929;];


%  Calculate the P matrix
P3 = calculate_camera_matrix(points_3D_image3, points_2D_image3);

% get the Rotation matrix R and Transaltion matrix T using  P & K matrix
[R3, T3] = decompose_camera_matrix(P3, K);

camera_position3 = T3;

% get the camera coordinates
camera_coordinates3 = get_camera_cooridnates(R3, points_3D_image3, T3);

fprintf("The camera coordinates are \n");
disp(camera_coordinates3);

% ========================================================================

figure;
% % Plot the 3D world coordinates
% plot3(points_3D_image1(:,1), points_3D_image1(:,2), points_3D_image1(:,3), 'k.', MarkerSize=12);
% hold on;
% plot3(points_3D_image2(:,1), points_3D_image2(:,2), points_3D_image2(:,3), 'k.', MarkerSize=12);
% hold on;
% plot3(points_3D_image3(:,1), points_3D_image3(:,2), points_3D_image3(:,3), 'k.', MarkerSize=12);
% hold on;


% ============================================= Camera position 1 ========
plot3(camera_coordinates1(:,1), camera_coordinates1(:,2), camera_coordinates1(:,3), 's', MarkerSize=6, MarkerEdgeColor ='k', MarkerFaceColor='r');
hold on;
% Draw lines between consecutive points
for i = 1:size(camera_coordinates1,1)
    nextIndex = mod(i, size(camera_coordinates1,1)) + 1; % This will loop back to 1 after 4
    line([camera_coordinates1(i,1), camera_coordinates1(nextIndex,1)], [camera_coordinates1(i,2), camera_coordinates1(nextIndex,2)], [camera_coordinates1(i,3), camera_coordinates1(nextIndex,3)], 'Color', 'red');
end
hold on;

% Mark the camera position
plot3(camera_position1(1), camera_position1(2), camera_position1(3), 'ro', 'LineWidth', 2, 'MarkerSize', 10);
hold on;

% ============================================= Camera position 2 ========
plot3(camera_coordinates2(:,1), camera_coordinates2(:,2), camera_coordinates2(:,3), 's', MarkerSize=6, MarkerEdgeColor ='k', MarkerFaceColor='c');
hold on;

% Draw lines between consecutive points
for i = 1:size(camera_coordinates2,1)
    nextIndex = mod(i, size(camera_coordinates2,1)) + 1; % This will loop back to 1 after 4
    line([camera_coordinates2(i,1), camera_coordinates2(nextIndex,1)], [camera_coordinates2(i,2), camera_coordinates2(nextIndex,2)], [camera_coordinates2(i,3), camera_coordinates2(nextIndex,3)], 'Color', 'cyan');
end
hold on;

% Mark the camera position
plot3(camera_position2(1), camera_position2(2), camera_position2(3), 'co', 'LineWidth', 2, 'MarkerSize', 10);
hold on;

% ============================================= Camera position 3 ========
plot3(camera_coordinates3(:,1), camera_coordinates3(:,2), camera_coordinates3(:,3), 's', MarkerSize=6, MarkerEdgeColor ='k', MarkerFaceColor='b');
hold on;
% Draw lines between consecutive points
for i = 1:size(camera_coordinates3,1)
    nextIndex = mod(i, size(camera_coordinates3,1)) + 1; % This will loop back to 1 after 4
    line([camera_coordinates3(i,1), camera_coordinates3(nextIndex,1)], [camera_coordinates3(i,2), camera_coordinates3(nextIndex,2)], [camera_coordinates3(i,3), camera_coordinates3(nextIndex,3)], 'Color', 'blue');
end
hold on;

% Mark the camera position
plot3(camera_position3(1), camera_position3(2), camera_position3(3), 'bo', 'LineWidth', 2, 'MarkerSize', 10);
hold on;

% Label the axes
xlabel('X');
ylabel('Y');
zlabel('Z');

% axis equal;
grid on;

% Title and legend
title('World and Camera cooridnate Frames');

% Release the hold on the current plot
hold off;


% genertaing the camera coordiantes 
function Xc = get_camera_cooridnates(R, points_3D_image1, Cw)
    
    % subtract the camera world position from each world position point
    % it will translates all points relative to the camera position
    translated_points = points_3D_image1 - Cw';
    
    % Multiply by the rotation matrix to get points in camera coordinates
    Xc = (R * translated_points')';
end


function [R, T] = decompose_camera_matrix(P, K)

    % Normalize the matrix P using the third row magnitude
    norm_val = norm(P(3, 1:3));
    P_norm = P / norm_val;
    
    % get the 3x3 matrix from P
    P_tilda = P_norm(:, 1:3);
    
    % use P~ and P to get the ROtation matrix R
    R = rotation_matrix(P, P_tilda);   
    
    % get the Transaltion vector T using K, R and 4th column of P
    % R is a rotation matrix so R transpose is equal to R inverse
    T = -inv(R) * inv(K) * P(:, 4);

    fprintf("Rotation Matrix R is\n");
    disp(R);
    fprintf("Translation vector T is \n");
    disp(T);

end

function R = rotation_matrix(P, P_tilda)

    % calculate the rotation matrix Rz by angle theta
    theta = atan2(P(3,1), P(3,2));
    Rz_theta = [cos(theta), sin(theta), 0;
                -sin(theta), cos(theta), 0;
                 0,  0, 1];

    % multiply the rotation matrix Rz with P~
    % Z-rotation does not affect the 3rd column of P~
    P_tilda_Rz = P_tilda * Rz_theta;

    % calculate the rotation matrix Rx by angle beta
    beta = atan2(P_tilda_Rz(3, 2), P_tilda_Rz(3, 3));
    Rx_beta = [1, 0, 0; 
               0, cos(beta), sin(beta); 
               0, -sin(beta), cos(beta)];

    % multiply the rotation matrix Rx_beta to P_tilda_Rz
    P_tilda_Rz_Rx = P_tilda_Rz * Rx_beta;

    % calcluate the rotation matrix Ry by angle gamma
    gamma = atan2(P_tilda_Rz_Rx(2, 1), P_tilda_Rz_Rx(2, 2));
    Ry_gamma = [cos(gamma), 0, sin(gamma); 
                0, 1, 0; 
                -sin(gamma), 0, cos(gamma)];

    % get the final rotation matrix R by multiplying Ry_gamma, Rx_beta, and
    % Rz_theta
    R = Rz_theta * Rx_beta * Ry_gamma;
end

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

    % display the calculated finite projective camera matrix
    fprintf("\n The finite projective camera matrix P is: \n");
    format long g;
    disp(P);

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

