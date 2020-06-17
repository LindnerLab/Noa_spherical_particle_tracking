% This code uses an convolutional method to determine the locations of
% circular particles in an image, given an expected particle radius, r.
% Particles can be separated from false positives by comparing any of the
% 'regionprops' properties to a reference value and selecting only the
% particles that have a property that is 'greater' or 'smaller' than the
% reference 'Value'. Examples of properties are 'Area' and 'Eccentricity'.
% More info can be found on https://www.mathworks.com/help/images/ref/regionprops.html
% The codes find circles for each frame of an image sequence
% Presenting the circles that are found in the center plane!
% Saving the coordinates of each circle for each image in a .txt file

clc;
close all;
clear all;

% Read the whole image sequence
image_folder = 'E:\Lars\Github\Noa_spherical_particle_tracking\RawData\';    %  Enter name of folder from which you want to upload pictures 
verbose = true;                                         % If verbose is true, the script will create and save a figure of the particle finding
r = 10;                                                 % Expected radius of the particles in pixels


SelectionCriteria = struct('Property',[],'Value',[],'Criteria',[]);
% The particle shouldn't be too small
SelectionCriteria(1).Property = 'Area';
SelectionCriteria(1).Value = pi*(0.8*r)^2;
SelectionCriteria(1).Criteria = 'Greater';

% The particle shouldn't be too large either
SelectionCriteria(2).Property = 'Area';
SelectionCriteria(2).Value = pi*(1.2*r)^2;
SelectionCriteria(2).Criteria = 'Smaller';

% The particle should be roughly round (0 is perfectly round, 1 is a line)
SelectionCriteria(3).Property = 'Eccentricity';
SelectionCriteria(3).Value = 0.4;
SelectionCriteria(3).Criteria = 'Smaller';


%% Below This No Input Are Required!
folder = pwd;
filenames = [dir(fullfile(image_folder, '*.tif'))];     % read all images with specified extention
total_images = numel(filenames);                        % count total number of photos present in that folder

mkdir raw_images_crop
mkdir particles_marked
mkdir particles_coordinates

% crop images
full_name = fullfile(image_folder, filenames(5).name);  % specify images names with full path and extension
Img_temp = imread(full_name);                           % Read img need to indicate crop
[Img_temp, rect] = imcrop(Img_temp);                    % Manually crop the example image and save the crop
Img_size = size(Img_temp);                              % Determine the size of the cropped image (needed for pre-allocating memory).

Img = zeros(Img_size(1),Img_size(2), total_images);     % Pre-allocate memory for the cropped images (is faster than allocating space when needed.
clear Img_temp Img_size                                 % Delete the temporary image, as it is no longer needed

for n = 1:total_images
    full_name = fullfile(image_folder, filenames(n).name);      % specify images names with full path and extension
    Img(:,:,n) = double(imcrop(imread(full_name),rect));        % Import, crop, convert the image to double and invert it
end

% calculate average image for background removal      
% Subtract background and invert image
Img_inv = 2^16 - Img + 1 - 5E4;
Img_inv(Img_inv < 0) = 0;
average_fig = mean(Img_inv,3);                          % Take the average over all the images (averages values of all pixels with the same z-index)
Img_bg = Img_inv - average_fig;                         % Subtract the background signal

% pass images through a filter, find circles and save x,y coordinates to *.txt files
for n = 1:total_images    
    % Find circles
    Img_bpass = bpass(Img_bg(:,:,n),0.5,30);            % filter - smoothing and subtracting the background (picture, 1, need to adjust to fit with particle diameter)
    [allcenters1] = FindParticlesConvolution(Img_bpass,r,SelectionCriteria);

    baseFileName = sprintf('%03d',n);                   % Create string with image number padded with zeros   
    if verbose
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
        colormap('gray'); imagesc(uint16(Img_bpass)); % display image
        set(gcf, 'Visible', 'off') %dont show figure while code is running
        title('Detected particles')

        % Indicate the found particles
        viscircles(allcenters1,ones(size(allcenters1,1),1)*r, 'EdgeColor','r','LineWidth',1);    
    %             %%If there are 2 particle sizes in the video             
    %             [centers2, radii2, metric2] = imfindcircles(b,[40,60]); %finds the particles with defined centers - the particles in the focus plane. [ radii range ]
    %             allcenters2 = centers2(:,:);
    %             allm2 = metric2(:);
    %             viscircles(allcenters2,allm2, 'EdgeColor','r');
        daspect([1 1 1])

        % Save the combined figure
        fullFileName = fullfile(folder,'particles_marked', [baseFileName,'.tif']);
        saveas(gcf, fullFileName,'tiff');
    end
    
    % Save cropped images
    fullFileName = fullfile(folder,'raw_images_crop',  [baseFileName,'.tif']);
    imwrite(uint16(Img(:,:,n)), fullFileName);
    
    % write the coordinates of the particles in each image to a seperate file
    fullFileName = fullfile(folder,'particles_coordinates', [baseFileName,'.txt']);
    dlmwrite(fullFileName,allcenters1,'delimiter','\t','precision',3);
end
