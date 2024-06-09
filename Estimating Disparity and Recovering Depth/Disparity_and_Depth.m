
%--------------  MiddleBury Images (takes approx 10 mins)  -------------%

% read the input rectified images 
leftImg = imread("images\middlebury_left.png");
rightImg = imread("images\middlebury_right.png");
% Set maximum disparity
max_disparity = 130; 
% focal length and baseline for depth calculation
focal_lenght = 1733.74;
baseline = 536.62;
%-----------------------------------------------------------------------%


% Uncomment below input lines to compute disaprity map for different input

%-----------  Elephant & Horse Images (takes approx 12 mins)  ----------%

% % read the input rectified images 
% leftImg = imread("images\I1Rect.jpg");
% rightImg = imread("images\I2Rect.jpg");
% % Set maximum disparity
% max_disparity = 72; 
% % focal length and baseline for depth calculation
% focal_lenght = 3060;
% baseline = 524;
% %-----------------------------------------------------------------------%


% Convert images to grayscale 
left_I1 = im2gray(leftImg);
right_I2 = im2gray(rightImg);

% Set penalty factors P1 and P2
P1 = 5;
P2 = 50;

% disparity map computation 
disparityMap = compute_disparity_map(left_I1, right_I2, max_disparity, P1, P2); 
% depth map computation
depthMap = disparityToDepthMap(disparityMap, baseline, focal_lenght);

% Visualization 
visualizeMaps(disparityMap, depthMap);



% computes the disaprty map for left and right rectified images 
function disparityMap = compute_disparity_map(leftImg, rightImg, maxDisparity, P1, P2)

    cost = compute_cost(leftImg, rightImg, maxDisparity);
    combined_cost = combine_cost(cost, P1, P2);
    disparityMap = get_disparity(combined_cost);
    disparityMap = post_process_disparity(disparityMap);
end


% calculate the bitwise XOR of 7X5 window for both images with current
% pixel as the center pixel
function cost = compute_cost(leftImg, rightImg, maxDisparity)

    leftIMG_count = block_count(leftImg);
    rightIMG_count = block_count(rightImg);

    %initialize the cost
    [rows, cols] = size(leftImg);
    cost = zeros(rows, cols, maxDisparity);

    % cost of each pixel pair for each disparity
    for d = 1:maxDisparity
        
        % get the XOR to get differing bits for each disparity d
        shiftedRight = [zeros(rows,d,'uint32'), rightIMG_count(:,1:end-d)];
        differingBits = bitxor(leftIMG_count, shiftedRight);
        
        % count the same bits ( number of 1)
        for bitIdx = 0:16  % 17-bit number
            cost(:,1:end-d,d)=cost(:,1:end-d,d)+double(bitget(differingBits(:,1:end-d),bitIdx+1));
        end
    end
end


function counts = block_count(img)
    [rows,cols] = size(img);
    counts = zeros(rows, cols, 'uint32');

    % Define the window size and center pixel position 
    windowSize = [7, 5];
    center = (windowSize(1) * windowSize(2)) / 2;

    % Pad the image to handle the borders
    paddedImg = padarray(img, [3 2], 'replicate', 'both');

    % Compute the pixel comparsion for each pixel as the center pixel 
    for i = 1:(center - 1)
        % Compute offsets for the center-symmetric pairs
        [rowOffset, colOffset] = ind2sub(windowSize, i);
        rowOffset = rowOffset - 4;
        colOffset = colOffset - 3;

        bit = uint32(paddedImg(4:(end-3), 3:(end-2)) >= ...
                     paddedImg(4 - rowOffset:(end-3) - rowOffset, ...
                               3 - colOffset:(end-2) - colOffset));

        counts = bitor(counts, bitshift(bit, i - 1));
    end
end



% Aggregates teh cost from different directiosn to smoothen the cost to
% depict real world continous intensity changes 
function combined_cost = combine_cost(cost, P1, P2)
    [height, width, ndisp] = size(cost);
     % 3 paths aggregation of costs
    directional_cost = zeros(height, width, ndisp, 5);

    % top-to-bottom path
    for x = 1:width
        for y = 2:height
            prevCosts = squeeze(directional_cost(y-1,x,:,1));
            minPrevCost = min(prevCosts);

            for d = 1:ndisp
                costSame = prevCosts(d);
                costPrevD = prevCosts(max(d-1,1))+P1;
                costNextD = prevCosts(min(d+1,ndisp))+P1;
                directional_cost(y,x,d,1) = cost(y,x,d) + min([costSame,costPrevD,costNextD,minPrevCost+P2])-minPrevCost;
            end
        end
    end

    % left-to-right path
    for y = 1:height
        for x = 2:width
            prevCosts = squeeze(directional_cost(y,x-1,:,2));
            minPrevCost = min(prevCosts);

            for d = 1:ndisp
                costSame = prevCosts(d);
                costPrevD = prevCosts(max(d-1, 1))+P1;
                costNextD = prevCosts(min(d+1, ndisp))+P1;
                directional_cost(y,x,d,2) = cost(y,x,d)+min([costSame,costPrevD,costNextD,minPrevCost+P2])-minPrevCost;
            end
        end
    end


    % top-left corner to the bottom-right corner of the image
    for y = 2:height
        for x = 2:width
            prevCosts = squeeze(directional_cost(y-1,x-1,:,3));
            minPrevCost = min(prevCosts);

            for d = 1:ndisp
                costSame = prevCosts(d);
                costPrevD = prevCosts(max(d-1,1))+P1;
                costNextD = prevCosts(min(d+1,ndisp))+P1;
                directional_cost(y,x,d,3) = cost(y, x, d)+min([costSame,costPrevD,costNextD,minPrevCost+P2])-minPrevCost;
            end
        end
    end

    %  top-right corner to the bottom-left corner of the image
    for y = 2:height
        for x = width-1:-1:1
            prevCosts = squeeze(directional_cost(y-1, x+1, :, 4));
            minPrevCost = min(prevCosts);

            for d = 1:ndisp
                costSame = prevCosts(d);
                costPrevD = prevCosts(max(d-1, 1)) + P1;
                costNextD = prevCosts(min(d+1, ndisp)) + P1;
                directional_cost(y, x, d, 4) = cost(y, x, d) + min([costSame, costPrevD, costNextD, minPrevCost + P2]) - minPrevCost;
            end
        end
    end

    % right-to-left path
    for y = 1:height
        for x = width-1:-1:1
            prevCosts = squeeze(directional_cost(y, x+1, :, 5));
            minPrevCost = min(prevCosts);

            for d = 1:ndisp
                costSame = prevCosts(d);
                costPrevD = prevCosts(max(d-1, 1)) + P1;
                costNextD = prevCosts(min(d+1, ndisp)) + P1;
                directional_cost(y, x, d, 5) = cost(y, x, d) + min([costSame, costPrevD, costNextD, minPrevCost + P2]) - minPrevCost;
            end
        end
    end

    % combined costs from different paths
    combined_cost = sum(directional_cost, 4);
end


% now given the disaprity values for each pixel pair, get the best disapriy
% and also make sure the pixel pairs are unique
function disparityMap = get_disparity(aggregated_cost)
    [height, width, num_disp] = size(aggregated_cost);
    %initialise the disparity map
    disparityMap = zeros(height, width);
    selectionMap = zeros(height, width, num_disp); 

    for y = 1:height
        for x = 1:width
            costs=squeeze(aggregated_cost(y, x, :));

            % penalize the disparities already selected
            % so that it won't be selected again to ensure uniqueness
            for  d=1:num_disp
                if selectionMap(y,x,d) == 1
                    costs(d) = costs(d)*20; 
                end
            end

            %select disparity with minimum cost
            [~,minIndex] = min(costs);
            disparityMap(y,x) = minIndex;

            %update the selection map
            selectionMap(y,x,minIndex) = 1;
        end
    end
end


function disparityMap = post_process_disparity(disparityMap)

    % Apply median filtering
    disparityMap = medfilt2(disparityMap, [5 5]);
    % handleing the invalid disparities values
    disparityMap(disparityMap < 0) = -1;
end


function depthMap = disparityToDepthMap(disparity_map, baseline, focal_length)

    depthMap = focal_length * baseline ./ double(disparity_map);
    % handle any invalid disparities
    depthMap(disparity_map == -1) = 0; 
end


function visualizeMaps(disparity_map, depth_map)

    % display the disparity map
    figure;
    imshow(disparity_map,[]);
    title("Disparity Map");
    colormap jet;
    colorbar;  

    % display the disparity map as a 3D surface

    %grid of coordinates
    [rows, cols] = size(depth_map);
    [X, Y] = meshgrid(1:cols, 1:rows);

    % disparity map as a surface
    figure;
    surf(X, Y, depth_map, disparity_map,'EdgeColor', 'none');
    zlabel('Depth');
    title('Disparity Map as a 3D Surface');
    colormap('jet');  
    colorbar;
    rotate3d on;

end
