% http://htmlpreview.github.io/?https://raw.githubusercontent.com/tinevez/simpletracker/master/publish/html/TestSimpleTracker.html
% This code plots the data from all .txt files found by the code
% Tracking_Noa_18042020 - giving an indication of all the particles that
% were detected finding their trajectories and plotting them
%
%The code finds the velocity of each particle along the trajectory and
% plots it

clc;close all;clear all;

pix2micro = 3.48101E-06;    %distance in pixel=316.0396, distance=1.10E-03m - pixel/m 
%average_flow_velocity = 7;        %ml/minute
U = 9.64E-02; %m/sec % !!!!!!! Adjust according to experiment
fps = 3000; % !!!!!!! Adjust according to experiment
dt = 1/fps;
w = 1.10E-03; %channel width 
nx = 1024;
ny = 512;

% Define x,y coordinates and normalize
x = 0.55;
y = 0.5;


v_mag = zeros(nx,ny);
v_all = zeros(nx,ny);


% Read the whole image sequence
file_folder = 'E:\Lars\Github\particles_coordinates';    %  Enter name of folder from which you want to upload pictures with full path
folder = pwd;
filenames = dir(fullfile(file_folder, '*.txt'));                               % read all images with specified extention
total_files = numel(filenames);                                                % count total number of photos present in that folder
%mkdir particles_marked

points = cell(total_files,1); % allocating a cell array for saving all the data collected in the txt files

% %figure (1) % Plots all detected particles in black 
% % for n = 1:total_files
% % 
% %     full_name= fullfile(file_folder, filenames(n).name);           % specify images names with full path and extension
% %     fid = fopen(full_name, 'r');                   % Read images
% %     Mydata = textscan(fid,'%f %f');           % Skip header lines
% %     fclose(fid);
% %     Mydata=cell2mat(Mydata); % this can be plotted to see all the particles detected in an overlay
% %     points(n,:) = {Mydata}; % contains all coordinates of all points found in a single file to perform tracking from here.   
% %     %plot((Mydata(:,1)*pix2micro/w)-1.6,(Mydata(:,2)*pix2micro/w)-0.8,'o', 'MarkerSize', 3, 'LineWidth', 1); hold on 
% %     daspect([1 1 1])
% %     set(gcf, 'color', 'w');
% %     xlim([-1.2 1.2])
% %     ylim([-0.5 0.5])
% %     xlabel('\it y / w ');  ylabel('\it z / w'); 
% % end
% % title('All detected particles')


for n = 1:total_files % Plots all detected particles in black (background for trajectories)

    full_name= fullfile(file_folder, filenames(n).name);           % specify images names with full path and extension
    fid = fopen(full_name, 'r');                   % Read images
    Mydata = textscan(fid,'%f %f');           % Skip header lines
    fclose(fid);
    Mydata=cell2mat(Mydata); % this can be plotted to see all the particles detected in an overlay
    points(n,:) = {Mydata}; %contains all coordinates of all points found in a single file to perform tracking from here.   
    plot((Mydata(:,1)*pix2micro/w)-x,((Mydata(:,2)*pix2micro/w)-y),'o', 'MarkerSize', 2, 'LineWidth', 1, 'color', 'k'); hold on % Hold the plot normalize axis and align vortex center to x,y (0,0)
    daspect([1 1 1])
    set(gcf, 'color', 'w');
    
end


%linking between the coordinates with Nearestneighbor method.
max_linking_distance = 100; %Inf; defines the maximal possible length of a trajectory
max_gap_closing = 1; % defines maximal linking gap
debug = true;

[ tracks, adjacency_tracks ] = simpletracker(points,...
    'MaxLinkingDistance', max_linking_distance, ...
    'MaxGapClosing', max_gap_closing, ...
    'Debug', debug);

adjacency_tracks_filter = adjacency_tracks(cellfun('length',adjacency_tracks)>=5); % determine the minimum length of a trajectory (to plot)
tracks_filter = adjacency_tracks(cellfun('length',adjacency_tracks)>=5); % determine the minimum length of a trajectory (to plot)

% Plot tracks - plot each track in a given color. Normally we would have to retrieve the points coordinates in the given |points| initiall cell
% arrat, for each point in frame. To skip this, we simple use the adjacency_tracks, that can pick points directly in the concatenated
% points array |all_points|.

n_tracks = numel(adjacency_tracks_filter);
colors = jet(n_tracks);
all_points = vertcat(points{:});
all_points = round(all_points);
all_track_points = zeros(length(all_points),2);

figure (1)
for i_track = 2 : 1 : n_tracks   % Use the adjacency tracks to retrieve the points coordinates.
       
    track = adjacency_tracks_filter{i_track};
    track_points = all_points(track, :);
    all_track_points(track,:) = all_points(track, :); 
    
    plot((track_points(:,1)*pix2micro/w)-x, (track_points(:, 2)*pix2micro/w)-y, '-', 'Color', colors(i_track, :),'MarkerSize', 2, 'LineWidth', 1); hold on; %normalize axis and align vortex center to x,y (0,0)
    daspect([1 1 1])
    set(gcf, 'color', 'w');
    xlim([-0.5 0.5])
    ylim([-0.5 0.5])
    xlabel('\it y / w ');  ylabel('\it z / w'); 
    title('Calculated trajectories')

end


  figure (2)  %plot the trajectories with the velocity of each segment
  numColors = 255; % resolution of color map
  cmap = jet(numColors);
  colormap(jet(numColors));
  
  for i_track = 1 : 1 : n_tracks% find velocity of each particle (linked to a trajectory) as a function of time

        track = adjacency_tracks_filter{i_track}; % find the filtered tracks
        track_points = all_points(track, :); % find all points from each track
        % v_average = zeros(size(i_track,adjacency_tracks_filter(1)));
        
        
        distance = zeros(size(track_points(1))); % allocate a space to save the distance a particle traveled between frames
        distanceX = zeros(size(track_points(1))); % allocate a space to save the distance X a particle traveled between frames
        distanceY = zeros(size(track_points(1))); % allocate a space to save the distance Y a particle traveled between frames
        v = zeros(size(track_points(1))); % allocate a space to find the velocity of a particle
        
        for k = 1 : length(track_points)-1    % length of each trajectory
            %             t(k) = (k-1)*dt; % find the total time a particle was detected
            distanceX(k) = ((abs  (track_points ((k+1),1) )) - (abs (track_points (k,1)) ))^2;  %calculate the distance of a particle between frames (k-2) is used for smoothing
            distanceY(k) = ((abs  (track_points ((k+1),2) )) - (abs (track_points (k,2)) ))^2;  %calculate the distance of a particle between frames (k-2) is used for smoothing
            distance = sqrt(distanceX + distanceY);
            v = (distance*pix2micro) / dt;
            track_points((k+1),3) = v(k);
            v_mag(track_points((k+1),1),track_points((k+1),2)) =   track_points((k+1),3);
            v_all(k,i_track) = v(k); % collect all velocities in a matrix
        end
   
         vMax = max(v_all(:))/U; % Find the maximal velocity of all the trajectories
         %vMax = 1;
         
         for k = 1 : length(track_points)-2
             vAvg = ((track_points((k+1),3) + track_points((k+2),3))/2)/U;  % calculate average velocity for this segment
             if   vAvg ~= 0 % this is to avoid data that is 0
                 if vAvg < 3 %vMax/3 %to avoid extremely high velocity values determined by a single vMax value
                     idx = ceil((vAvg/(3)*numColors)); % look up color for current velocity normalize by maximum velocity  %    idx = ceil(vAvg/vMax*numColors); % normalize by maximum velocity
                     color = cmap(idx,:);
                     plot([(track_points(k+1,1)*pix2micro/w)-x,(track_points((k+2),1)*pix2micro/w)-x],[(track_points(k+1,2)*pix2micro/w)-y,(track_points((k+2),2)*pix2micro/w)-y],'o-','MarkerSize', 4, 'Color',color, 'LineWidth', 2);   hold on;
                 end
             end
         end
  end
 
  
            h = colorbar;
            ylabel(h, '\it v_{mag} / U', 'FontName', 'Times')
            caxis([0  vMax/3]) 
            daspect([1 1 1])
            set(gcf, 'color', 'w');
            set(gca,'FontSize', 16, 'FontName', 'Times', 'LineWidth',1.5)
            xlim([-0.5 0.5])
            ylim([-0.5 0.5])
            set(gca,'YTick',-0.5:0.5:0.5);
            xlabel('\it y / w ');  ylabel('\it z / w');
            %   legend ();
  %          title('Trajectories with velocity of each segment')

            
 % Collect points only from the tracks that were filtered 
 all_track_points_x  = nonzeros(all_track_points(:,1));%(all_track_points ~= 0);
 all_track_points_y  = nonzeros(all_track_points(:,2));%(all_track_points ~= 0);

 
%figure(3)
% plot the density of the particles
%[bins in x, bins in y] 1024 x 512 pixel (25 pixel per bin - approx size of a particle)
   %     densityplot((all_points(1:2:end,1)*pix2micro/w)-1.7,(all_points(1:2:end,2)*pix2micro/w)-0.8, [51,25]); %https://www.mathworks.com/matlabcentral/fileexchange/65166-densityplot-x-y-varargin 
   densityplot((all_track_points_x(1:1:end)*pix2micro/w)-x,(all_track_points_y(1:1:end)*pix2micro/w)-y, [25,20]); %https://www.mathworks.com/matlabcentral/fileexchange/65166-densityplot-x-y-varargin
   daspect([1 1 1])
   h = colorbar;
   caxis([0 5])
   ylabel(h, 'Number of particles', 'FontName', 'Times')
   set(gcf, 'color', 'w');
   set(gca,'Color',[0.0 0.0 0.6])
   set(gca,'FontSize', 16, 'FontName', 'Times', 'LineWidth',1.5)
   xlim([-0.4 0.35])
   ylim([-0.35 0.37])
   set(gca,'YTick',-0.35:0.35:0.35);
   set(gca,'XTick',-0.35:0.35:0.35);
   xlabel('\it y / w ');  ylabel('\it z / w');
%   title('Particle density')

