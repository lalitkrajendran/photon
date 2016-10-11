% this is the filepath where the images are stored
filepath = ['/home/barracuda/a/lrajendr/Projects/camera_simulation/results/images/bos/error-analysis/dot-size/processing/' case_name '/cropped-images/'];

files = dir([filepath '*.tif']);

img1 = imread([filepath files(1).name]);
img2 = imread([filepath files(2).name]);

T = zeros(size(img1,1), size(img1,2), 3);
T(:,:,1) = img1;
T(:,:,2) = img2;

figure
imshow(T)
