%% Feature Matchings
I_left = rgb2gray(imread('images\IMG_left.jpg'));
I_right = rgb2gray(imread('images\IMG_right.jpg'));

numOfFeatures = 0;  % <=0 means keep all features

points_left = detectSIFTFeatures(I_left);
points_right = detectSIFTFeatures(I_right);

if numOfFeatures > 0
    points_left = points_left.selectStrongest(numOfFeatures);
    points_right = points_right.selectStrongest(numOfFeatures);
end

[features_left, featurepoints_left] = extractFeatures(I_left, points_left);
[features_right, featurepoints_right] = extractFeatures(I_right, points_right);

indexPairs = matchFeatures(features_left, features_right);
matchedPoints_left = featurepoints_left(indexPairs(:, 1), :);
matchedPoints_right = featurepoints_right(indexPairs(:, 2), :);

figure; 
showMatchedFeatures(I_left,I_right,matchedPoints_left,matchedPoints_right);

matchedPoints_left = matchedPoints_left.Location;
matchedPoints_right = matchedPoints_right.Location;
N = size(matchedPoints_left, 1);

%% Estimation of F Matrix w/ RANSAC
TrialTimes = 1000;
DistanceThreshold = 5;
CountThreshold = 128;

BestFSoFar = [];
BestConsensusSet_Left = [];
BestConsensusSet_Right = [];
ShortestDistSqSoFar = intmax;

for j = 1 : TrialTimes
    % Randomly pick 8 points
    p = randperm(size(matchedPoints_left, 1), 8);
    SelectedMatchedPoints_Left = matchedPoints_left(p, :);
    SelectedMatchedPoints_Right = matchedPoints_right(p, :);
    
    F = MakeFMatrix_EightPoints(SelectedMatchedPoints_Left, SelectedMatchedPoints_Right);

    ConsensusSet_Left = [];
    ConsensusSet_Right = [];

    for i = 1 : size(matchedPoints_left, 1)
        DistSq = GetDistanceSq(matchedPoints_left(i, :), matchedPoints_right(i, :), F);
        if DistSq <= DistanceThreshold
            ConsensusSet_Left = [ConsensusSet_Left; matchedPoints_left(i, :)];
            ConsensusSet_Right = [ConsensusSet_Right; matchedPoints_right(i, :)];
        end
        
    end

    if size(ConsensusSet_Left, 1) >= CountThreshold
        % Refit the Model
        F = MakeFMatrix(ConsensusSet_Left, ConsensusSet_Right);
        TotalDistSq = GetTotalDistanceSq(ConsensusSet_Left, ConsensusSet_Right, F);
        if TotalDistSq < ShortestDistSqSoFar
            BestConsensusSet_Left = ConsensusSet_Left;
            BestConsensusSet_Right = ConsensusSet_Right;
            BestFSoFar = F;
            ShortestDistSqSoFar = TotalDistSq;
        end
    end
end

BestFSoFar = BestFSoFar / norm(BestFSoFar);
disp(ShortestDistSqSoFar);
disp(BestFSoFar);
disp(rank(BestFSoFar));

FMatrix = BestFSoFar;   % Our Final Result 
%% Rectification
[U, S, V] = svd(FMatrix);
e1 = V(:, 3);
e1 = e1 / e1(3);
e2 = U(:, 3);
e2 = e2 / e2(3);

[H1Matrix, H2Matrix] = MakeH1H2Matrix(e1, e2, FMatrix);

I_left = imread('images\IMG_left.jpg');
I_left = im2double(I_left(:, :, 1));
I_right = imread('images\IMG_right.jpg');
I_right = im2double(I_right(:, :, 1));

Permutation = randperm(size(BestConsensusSet_Left, 1), 10);
SelectedMatchedPoints_Left = BestConsensusSet_Left(Permutation, :);
SelectedMatchedPoints_Right = BestConsensusSet_Right(Permutation, :);

epiLines_left = epipolarLine(FMatrix', SelectedMatchedPoints_Right);
points_left = lineToBorderPoints(epiLines_left,size(I_left));

epiLines_right = epipolarLine(FMatrix, SelectedMatchedPoints_Left);
points_right = lineToBorderPoints(epiLines_right,size(I_right));

figure;
imshow(I_left);
viscircles(SelectedMatchedPoints_Left, 10);
line(points_left(:,[1,3])',points_left(:,[2,4])', 'Color', "red", 'LineWidth', 1);
frame = getframe(gca);
I_leftWithLines = frame2im(frame);
I_LeftRectified = RectifyImage(I_leftWithLines, H1Matrix);

figure;
imshow(I_right);
viscircles(SelectedMatchedPoints_Right, 10);
line(points_right(:,[1,3])',points_right(:,[2,4])', 'Color', "red", 'LineWidth', 1);
frame = getframe(gca);
I_rightWithLines = frame2im(frame);
I_RightRectified = RectifyImage(I_rightWithLines, H2Matrix);

RectifiedImages = zeros(size(I_LeftRectified, 1), size(I_LeftRectified, 2) * 2, size(I_LeftRectified, 3), 'uint8');
RectifiedImages(:, 1 : size(I_LeftRectified, 2), :) = I_LeftRectified;
RectifiedImages(:, size(I_LeftRectified, 2) + 1 : end, :) = I_RightRectified;
figure;
imshow(RectifiedImages);
%% FUNCTIONS

function[I_Rectified] = RectifyImage(I, H)
    I_Rectified = zeros(size(I, 1), size(I, 2), 3);
    I_Rectified = uint8(I_Rectified);

    for i = 1 : size(I, 1)
        for j = 1 : size(I, 2)
            Index = inv(H) * [j, i, 1].';
            Index = Index / Index(3);
            Index = floor(Index);

            if 1 <= Index(1) && Index(1) <= size(I, 2)
                if 1 <= Index(2) && Index(2) <= size(I, 1)
                    I_Rectified(i, j, :) = I(Index(2), Index(1), :);
                end
            end
        end
    end
end

function [H1, H2] = MakeH1H2Matrix(e1, e2, F)
    H1 = [1, 0, 0;
        -e1(2) / e1(1), 1, 0;
        -1 / e1(1), 0, 1];

    B = [-H1(3,1),0,0,H1(2,1),0,0,-F(1, 1);
        0,0,0,1,0,0,-F(1,2);
        -1,0,0,0,0,0,-F(1,3);
        0,-H1(3,1),0,0,H1(2,1),0,-F(2,1);
        0,0,0,0,1,0,-F(2,2);
        0,-1,0,0,0,0,-F(2,3);
        0,0,-H1(3,1),0,0,H1(2,1),-F(3,1);
        0,0,0,0,0,1,-F(3,2);
        0,0,-1,0,0,0,-F(3,3)];

    [~, S, V] = svd(B);
    p_hat = V(:, size(V, 2));

    H2 = [1, 0, 0;
          p_hat(1:3).';
          p_hat(4:6).'];
end

function [DistSq] = GetDistanceSq(matchedPoint_l, matchedPoint_r, F)
    DistSq = ([matchedPoint_r, 1] * F * transpose([matchedPoint_l, 1])) ^ 2;
end

function [TotalDistSq] = GetTotalDistanceSq(matchedPoints_l, matchedPoints_r, F)
    TotalDistSq = 0;
    for i = 1 : size(matchedPoints_l, 1)
        TotalDistSq = TotalDistSq + GetDistanceSq(matchedPoints_l(i, :), matchedPoints_r(i, :), F);
    end
end

function [A] = MakeAMatrix(matchedPoints_l, matchedPoints_r)
    A = zeros(size(matchedPoints_l, 1) + 1, 9);
    for i = 1:size(matchedPoints_l, 1)
        x1 = matchedPoints_l(i, 1);
        y1 = matchedPoints_l(i, 2);
        x2 = matchedPoints_r(i, 1);
        y2 = matchedPoints_r(i, 2);

        A(i, :) = [x1*x2, y1*x2, x2, x1*y2, y1*y2, y2, x1, y1, 1];
    end
end

function [A_normalized, M1, M2] = MakeAMatrix_Normalized(matchedPoints_l, matchedPoints_r)
    %Find their sigma1, sigma2, x2d_bar, y2d_bar, x3d_bar, y3d_bar, z3d_bar
    xl_bar = mean(matchedPoints_l(:, 1));
    yl_bar = mean(matchedPoints_l(:, 2));
    xr_bar = mean(matchedPoints_r(:, 1));
    yr_bar = mean(matchedPoints_r(:, 2));
    
    xl_power2sum = sum((matchedPoints_l(:, 1) - xl_bar).^2);
    yl_power2sum = sum((matchedPoints_l(:, 2) - yl_bar).^2);
    xr_power2sum = sum((matchedPoints_r(:, 1) - xr_bar).^2);
    yr_power2sum = sum((matchedPoints_r(:, 2) - yr_bar).^2);
    sigma1 = sqrt((xl_power2sum + yl_power2sum)/(2*size(matchedPoints_l, 1)));
    sigma2 = sqrt((xr_power2sum + yr_power2sum)/(2*size(matchedPoints_l, 1)));
    
    %Build matrix M1 and M2
    tmp1 = diag([1/sigma1, 1/sigma1, 1]);
    tmp2 = diag([1,1,1]);
    tmp2(1, 3) = -xl_bar;
    tmp2(2, 3) = -yl_bar;
    M1 = tmp1 * tmp2;
    
    tmp1 = diag([1/sigma2, 1/sigma2, 1]);
    tmp2 = diag([1,1,1]);
    tmp2(1, 3) = -xr_bar;
    tmp2(2, 3) = -yr_bar;
    M2 = tmp1 * tmp2;

    A_normalized = ones(size(matchedPoints_l, 1), 9);
    for i = 1:size(matchedPoints_l, 1)
        point_l_normed = transpose(M1 * transpose([matchedPoints_l(i, :), 1]));
        point_r_normed = transpose(M2 * transpose([matchedPoints_r(i, :), 1]));
        
        x1 = point_l_normed(:, 1);
        y1 = point_l_normed(:, 2);
        x2 = point_r_normed(:, 1);
        y2 = point_r_normed(:, 2);

        A_normalized(i, :) = [x1*x2, y1*x2, x2, x1*y2, y1*y2, y2, x1, y1, 1];
    end
end

function [F] = MakeFMatrix_EightPoints(matchedPoints_l, matchedPoints_r)
    A = MakeAMatrix(matchedPoints_l, matchedPoints_r);
        
    [~, ~, V] = svd(A);
    
    F = transpose(reshape(V(:, size(V, 2)), 3, 3));
    F = ForceRank2(F);
end

function [F] = MakeFMatrix(matchedPoints_l, matchedPoints_r)
    [A_normalized, M1, M2] = MakeAMatrix_Normalized(matchedPoints_l, matchedPoints_r);
    
    [~, ~, V] = svd(A_normalized);
    F_normalized = transpose(reshape(V(:, size(V, 2)), 3, 3));
    F_normalized = ForceRank2(F_normalized);

    F = transpose(M2) * F_normalized * M1;
end

function [F_rank2] = ForceRank2(F)
    [U, S, V] = svd(F);
    S(3, 3) = 0;
    F_rank2 = U * S * V';
end
