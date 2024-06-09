
% standard deviation (sigma) and filter size
mean = 0;
sigma = 3;
gaussian = @(x, s) (1/(sqrt(2*pi)*s))*exp(-((x-mean).^2)/(2*(s.^2))); 

% 1D Gaussian kernel along the rows centered at X,Y=0,0
%multiply by 3 so final filter is square --> we want a 6sigma+1 sized filter
x = [-3*sigma:3*sigma];
y = gaussian(x, sigma);
gaussian_filter = y./sum(y); %normalize so they sum to 1

%second derivative filter
second_deriv = [1 -2 1];
filter_conv = conv(gaussian_filter, second_deriv, "same");

%create array
filter = repmat(filter_conv, (6*sigma)+1, 1);

%plots
figure
subplot(1, 2, 1)
imagesc(filter)
axis square
colorbar
subplot(1, 2, 2)
surf(filter)
colorbar
sgtitle('2D Gaussian filter convolved with second derivative filter (\sigma = 4)')


sigma_y = 5;
x2 = -2*sigma_y:(4*sigma_y)/(6*sigma):2*sigma_y;
y2 = gaussian(x2, sigma_y);
gaussian_filter_y = y2./sum(y2);

%multiply entry-wise
filter2 = filter .*(gaussian_filter_y');

%plot
figure
subplot(1, 2, 1)
imagesc(filter2)
axis square
colorbar
subplot(1, 2, 2)
surf(filter2)
colorbar
sgtitle('2D Gaussian filter convolved with second derivative filter multiplied by Gaussian in y (\sigma_y = 5)')


%put original filter in a larger matrix
filter2_larger = zeros(size(filter2, 1)+floor(size(filter2, 1)/2));
lower_bound = floor(size(filter2_larger, 1)/2) - floor(size(filter2, 1)/2)+1;
upper_bound = floor(size(filter2_larger, 1)/2) + floor(size(filter2, 1)/2)+1; 
filter2_larger(lower_bound:upper_bound, lower_bound:upper_bound) = filter2;

%make new templates
filter_45 = zeros(size(filter2, 1));
filter_90 = zeros(size(filter2, 1));
filter_135 = zeros(size(filter2, 1));

diff1 = ceil(size(filter_45, 1)/2); %to account for index shifts in new matrix
diff2 = ceil(size(filter2_larger, 1)/2); %to account for index shifts in old matrix

%90 degrees
angle = pi/2;
rotationMatrix = [cos(angle) sin(angle); -sin(angle) cos(angle)];
for i=1:size(filter_90, 1)
    for j=1:size(filter_90, 2)
        coords = round(rotationMatrix*[i-diff1;j-diff1]);
        coords = coords + diff2;
        filter_90(i, j) = filter2_larger(coords(1), coords(2));
    end;
end

%45 degrees
angle = pi/4;
rotationMatrix = [cos(angle) sin(angle); -sin(angle) cos(angle)];
for i=1:size(filter_45, 1)
    for j=1:size(filter_45, 2)
        coords = round(rotationMatrix*[i-diff1;j-diff1]);
        coords = coords + diff2;
        filter_45(i, j) = filter2_larger(coords(1), coords(2));
    end
end

%135 degrees
angle = 3*pi/4;
rotationMatrix = [cos(angle) sin(angle); -sin(angle) cos(angle)];
for i=1:size(filter_135, 1)
    for j=1:size(filter_135, 2)
        coords = round(rotationMatrix*[i-diff1;j-diff1]);
        coords = coords + diff2;
        filter_135(i, j) = filter2_larger(coords(1), coords(2));
    end
end


%plots
figure
subplot(1, 2, 1)
imagesc(filter_90)
axis square
colorbar
subplot(1, 2, 2)
surf(filter_90)
colorbar
sgtitle('Filter rotated at 90 degrees')

figure
subplot(1, 2, 1)
imagesc(filter_45)
axis square
colorbar
subplot(1, 2, 2)
surf(filter_45)
colorbar
sgtitle('Filter rotated at 45 degrees')

figure
subplot(1, 2, 1)
imagesc(filter_135)
axis square
colorbar
subplot(1, 2, 2)
surf(filter_135)
colorbar
sgtitle('Filter rotated at 135 degrees')



%read image
im_rgb = im2double(imread("Paolina.png"));
I = rgb2gray(im_rgb);

%filter images with filter
I_0 = conv2(I, filter2 , "same");
I_45 = conv2(I, filter_45, "same");
I_90 = conv2(I, filter_90, "same");
I_135 = conv2(I, filter_135, "same");

%create the binary image
threshold = 0.0005; %set threshold to remove weird artifacts

%find zero crossings in vertical direction --> compare x to x-1
shift_left_0 = circshift(I_0, 1, 2);
shift_down_0 = circshift(I_0, 1, 1);
shift_diagdown_0 = circshift(I_0, [-1 1]);
shift_diagup_0 = circshift(I_0, [1 1]);

I_0_bin = (I_0.*shift_left_0<0 & abs(I_0-shift_left_0)>threshold) | (I_0.*shift_down_0<0 & abs(I_0-shift_down_0)>threshold) ...
        | (I_0.*shift_diagdown_0<0 & abs(I_0-shift_diagdown_0)>threshold) | (I_0.*shift_diagup_0<0 & abs(I_0-shift_diagup_0)>threshold); 

%zero crossings in horizontal direction --> compare y to y-1
shift_left_90 = circshift(I_90, 1, 2);
shift_down_90 = circshift(I_90, 1, 1);
shift_diagdown_90 = circshift(I_90, [-1 1]);
shift_diagup_90 = circshift(I_90, [1 1]);

I_90_bin = (I_90.*shift_left_90<0 & abs(I_90-shift_left_90)>threshold) | (I_90.*shift_down_90<0 & abs(I_90-shift_down_90)>threshold) ...
        | (I_90.*shift_diagdown_90<0 & abs(I_90-shift_diagdown_90)>threshold) | (I_90.*shift_diagup_90<0 & abs(I_90-shift_diagup_90)>threshold); 

%zeros crossings in diagonal direction
shift_left_45 = circshift(I_45, 1, 2);
shift_down_45 = circshift(I_45, -1, 1);
shift_diagdown_45 = circshift(I_45, [-1 1]);
shift_diagup_45 = circshift(I_45, [1 1]);

I_45_bin = (I_45.*shift_left_45<0 & abs(I_45-shift_left_45)>threshold) | (I_45.*shift_down_45<0 & abs(I_45-shift_down_45)>threshold) ...
        | (I_45.*shift_diagdown_45<0 & abs(I_45-shift_diagdown_45)>threshold) | (I_45.*shift_diagup_45<0 & abs(I_45-shift_diagup_45)>threshold); 

%zero crossings in diagonal direction
shift_left_135 = circshift(I_135, 1, 2);
shift_down_135 = circshift(I_135, -1, 1);
shift_diagdown_135 = circshift(I_135, [-1 1]);
shift_diagup_135 = circshift(I_135, [1 1]);

I_135_bin = (I_135.*shift_left_135<0 & abs(I_135-shift_left_135)>threshold) | (I_135.*shift_down_135<0 & abs(I_135-shift_down_135)>threshold) ...
        | (I_135.*shift_diagdown_135<0 & abs(I_135-shift_diagdown_135)>threshold) | (I_135.*shift_diagup_135<0 & abs(I_135-shift_diagup_135)>threshold); 

%overlay
figure
tiledlayout(2, 2,'TileSpacing','tight')
nexttile
RGB_Image1 = cat(3, I_0_bin, I_0_bin, I_0_bin)*255; %creates the three color channels
RGB_Image1(:,:,1)*100;
RGB_Image1(:,:,2)*100;
RGB_Image1(:,:,3) = 0;
overlaid = imfuse(I, RGB_Image1,'blend');
imshow(overlaid)
title("Image with 0 degree rotated filter")

nexttile
RGB_Image3 = cat(3, I_45_bin, I_45_bin, I_45_bin)*255;
RGB_Image3(:,:,1) = 0;
RGB_Image3(:,:,2)*100;
overlaid = imfuse(I, RGB_Image3,'blend');
imshow(overlaid)
title("Image with 45 degree rotated filter")

nexttile
RGB_Image2 = cat(3, I_90_bin, I_90_bin, I_90_bin)*255;
RGB_Image2(:,:,1) = 0;
RGB_Image2(:,:,3) = 0;
overlaid = imfuse(I, RGB_Image2,'blend');
imshow(overlaid)
title("Image with 90 degree rotated filter")

nexttile
RGB_Image4 = cat(3, I_135_bin, I_135_bin, I_135_bin)*255;
RGB_Image4(:,:,1)*100;
RGB_Image4(:,:,2) = 0;
overlaid = imfuse(I, RGB_Image4,'blend');
imshow(overlaid)
title("Image with 135 degree rotated filter")