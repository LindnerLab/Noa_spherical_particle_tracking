% This code is based on: http://site.physics.georgetown.edu/matlab/tutorial.html
% and https://www.mathworks.com/help/images/ref/imfindcircles.html
% The codes find circles for each fram of an image sequence
% Presenting the circles that are found in the center plane!
% Saving the coordinates of each circle for each image in a .txt file

clc;
close all;
clear all;

% Read the whole image sequence
image_folder = 'E:\Lars\Github\Noa_spherical_particle_tracking\RawData\';    %  Enter name of folder from which you want to upload pictures 
folder = pwd;
filenames = [dir(fullfile(image_folder, '*.tif'))];                               % read all images with specified extention
total_images = numel(filenames);                                                % count total number of photos present in that folder
mkdir raw_images_crop
mkdir particles_marked
mkdir particles_coordinates

% crop images
full_name = fullfile(image_folder, filenames(5).name);           % specify images names with full path and extension
Img_temp = imread(full_name);      
[Img_temp, rect] = imcrop(Img_temp);
Img_size = size(Img_temp);

Img = zeros(Img_size(1),Img_size(2), total_images);
Img_inv = zeros([size(Img_temp), total_images]);
clear Img_temp Img_size

for n = 1:total_images
    full_name = fullfile(image_folder, filenames(n).name);           % specify images names with full path and extension
    Img(:,:,n) = double(imcrop(imread(full_name),rect));
end

% calculate average image for background removal      
average_fig = mean(Img,3);

% pass images through a filter, find circles and save x,y coordinates to *.txt files
for n = 1:total_images
    % Subtract background and invert image
    Img_inv(:,:,n) = 2^16-Img(:,:,n)+average_fig-1;
    
    % Find circles
    Img_bpass = bpass(Img_inv(:,:,n),0.5,30); % filter - smoothing and subtracting the background (picture, 1, need to adjust to fit with particle diameter)
    [allcenters1, radii1, ~] = imfindcircles(Img_bpass,[7,12]); %finds the particles with defined centers - the particles in the focus plane. [ radii range ]
    
    % Plot original image
    figure(n);
    subplot(2,1,1)
    colormap('gray'),imagesc(Img(:,:,n))
    title('Original Image')
    daspect([1 1 1])
    set(gcf, 'Visible', 'off') %dont show figure while code is running
    set(gcf, 'color', 'w');
    
    % Plot inverted and filtered image (with background removed)
    subplot(2,1,2)
    colormap('gray'); imagesc(Img_bpass); % display image
    set(gcf, 'Visible', 'off') %dont show figure while code is running
    title('Detected particles')
    
    % Indicate the found particles
    viscircles(allcenters1,radii1, 'EdgeColor','r');    
%             %%If there are 2 particle sizes in the video             
%             [centers2, radii2, metric2] = imfindcircles(b,[40,60]); %finds the particles with defined centers - the particles in the focus plane. [ radii range ]
%             allcenters2 = centers2(:,:);
%             allm2 = metric2(:);
%             viscircles(allcenters2,allm2, 'EdgeColor','r');
    daspect([1 1 1])
    
    % Save the combined figure
    baseFileName = sprintf('%03d.tif',n);
    fullFileName = fullfile(folder,'particles_marked', baseFileName);
    saveas(gcf, fullFileName,'tiff');
    
    % Save cropped images
    fullFileName = fullfile('raw_images_crop',  baseFileName);
    imwrite(uint16(Img(:,:,n)), fullFileName);
    
    % write the coordinates of the particles in each image to a seperate file
    baseFileName1 = sprintf('%03d.txt',n);
    fullFileName1 = fullfile(folder,'particles_coordinates', baseFileName1);
    dlmwrite(fullFileName1,allcenters1,'delimiter','\t','precision',3);
end
