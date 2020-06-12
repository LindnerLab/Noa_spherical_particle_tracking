% This code is based on: http://site.physics.georgetown.edu/matlab/tutorial.html
% and https://www.mathworks.com/help/images/ref/imfindcircles.html
% The codes find circles for each fram of an image sequence
% Presenting the circles that are found in the center plane!
% Saving the coordinates of each circle for each image in a .txt file

clc;close all;clear all;

% Read the whole image sequence
image_folder = '\Users\bnoa5\OneDrive\Desktop\Particle_flow_analysis_Noa30032020\Results_analysis_June_2020\water-glycerol\80um_MPs\7ml_min\raw_images';    %  Enter name of folder from which you want to upload pictures 
folder = pwd;
filenames = dir(fullfile(image_folder, '*.tif'));                               % read all images with specified extention
total_images = numel(filenames);                                                % count total number of photos present in that folder
mkdir raw_images_crop
mkdir particles_marked
mkdir particles_coordinates



% % crop images
full_name = fullfile(image_folder, filenames(5).name);           % specify images names with full path and extension
b = (imread(full_name));      
[J, rect] = imcrop(b);


for n = 1:total_images
        full_name = fullfile(image_folder, filenames(n).name);           % specify images names with full path and extension
        b = (imread(full_name));
        b = imcrop(b,rect) ;           % crop image
        [filepath,name,ext] = fileparts(filenames(n).name) ;
        baseFileName = sprintf('%03d.tif',n);
        fullFileName = fullfile('raw_images_crop',  baseFileName);
        imwrite(b, fullFileName);
end



% calculate average image for background removal
image_folder2 = '\Users\bnoa5\OneDrive\Desktop\Particle_flow_analysis_Noa30032020\Results_analysis_June_2020\water-glycerol\80um_MPs\7ml_min\raw_images_crop';    %  Enter name of folder from which you want to upload pictures 
filenames = dir(fullfile(image_folder2, '*.tif'));                               % read all images with specified extention
total_images = numel(filenames);                                                % count total number of photos present in that folder

    for n = 1:total_images
            full_name = fullfile(image_folder2, filenames(n).name);           % specify images names with full path and extension
            a = double(imread(full_name));      
            if n ==1
                sumimage = a;
            else
                sumimage = sumimage + a;
            end
    end
    
    average_fig = sumimage / total_images;
 %   colormap('gray'),imagesc(average_fig)

% pass images through a filter, find circles and save x,y coordinates to *.txt files

for n = 1:total_images
h = figure; 
    
subplot(2,1,1)
            full_name= fullfile(image_folder2, filenames(n).name);           % it will specify images names with full path and extension
            a = double(imread(full_name));                                                   % Read images
            colormap('gray'),imagesc(a)
            title('Original Image')
            a = a - average_fig; %remove background
            a = 255-a;  % convert black to white
            daspect([1 1 1])
            set(gcf, 'Visible', 'off') %dont show figure while code is running
            set(gcf, 'color', 'w');

subplot(2,1,2)
            a = bpass(a,0.5,30); % filter - smoothing and subtracting the background (picture, 1, need to adjust to fit with particle diameter)
            colormap('gray'); imagesc(a); % display image
            title('Detected particles')

% % find circles
            [centers1, radii1, metric1] = imfindcircles(a,[7,12]); %finds the particles with defined centers - the particles in the focus plane. [ radii range ]
            allcenters1 = centers1(:,:);
            allm1 = metric1(:);
            viscircles(allcenters1,allm1, 'EdgeColor','r');    
%             %%If there are 2 particle sizes in the video             
%             [centers2, radii2, metric2] = imfindcircles(b,[40,60]); %finds the particles with defined centers - the particles in the focus plane. [ radii range ]
%             allcenters2 = centers2(:,:);
%             allm2 = metric2(:);
%             viscircles(allcenters2,allm2, 'EdgeColor','r');
            daspect([1 1 1])

            baseFileName = sprintf('%03d.tif',n);
            fullFileName = fullfile(folder,'particles_marked', baseFileName);
            figure(n);
            export_fig(fullFileName,n);
            set(gcf, 'Visible', 'off') %dont show figure while code is running

% write the coordinates of the particles in each image to a seperate file
            baseFileName1 = sprintf('%03d.txt',n);
            fullFileName1 = fullfile(folder,'particles_coordinates', baseFileName1);
            dlmwrite(fullFileName1,centers1,'delimiter','\t','precision',3);
  
end

