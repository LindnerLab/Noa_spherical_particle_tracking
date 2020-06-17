function [Pparticles] = FindParticlesConvolution(Img,r,SelectionCriteria)
%FINDPARTICLESCONVOLUTION finds circular particles in pre-treated images.
%   FindParticlesConvolution finds particles of radius r in image Img that
%   meet all the requirements set in SelectionCriteria
%
%   Outputs
%       Pparticles - 2D-array where the first column denotes the
%               x-coordinates of all the particles found and the second
%               column denotes the y-coordinates of the particles found.
%   
%   Inputs
%       Img   - 2D-array double, with dark background and light features.
%               Make sure that the features are well separated, as
%               regionprops will count overlapping/touching regions as one
%               particle.
%
%       r     - Double, expected radius of the particles in pixels.
%
%       SelectionCriteria - Struct with 3 fields: 'Property', 'Value'
%               and 'Criteria', where 'Property' should be
%               a string of one of the regionprops properties (see
%               https://www.mathworks.com/help/images/ref/regionprops.html
%               for more information on this), 'Value' should be a double
%               to which the regionproperty 'Property' of the regions will
%               be compared, and 'Criteria' should be a string equal to
%               either 'greater' or 'smaller', to indicate if the
%               regionproperty 'Property' should be smaller or greater than
%               the reference value 'Value'. This function supports
%               multiple selection criteria, where each selection criteria
%               is a new entry in the struct.

    Img_size = size(Img);                           % Determine the size of the image, needed to determine the size of the mask (as they need to be equal in size).
    
    xCenter = r;                                    % X-position of the center of the object to look for. It is centered at (r,r), meaning that the peaks in the convoluted image will be at (x+r,y+r), instead of (x,y)!
    yCenter = r;                                    % Y-position of the center of the object to look for.
    xCircle = r*cos(0:pi/15:2*pi) + xCenter;        % X-coordinates of a circle centered at (xCenter, yCenter) with radius r.
    yCircle = r*sin(0:pi/15:2*pi) + yCenter;        % Y-coordinates of a circle centered at (xCenter, yCenter) with radius r.
    
    [X, Y] = meshgrid(1:Img_size(2),1:Img_size(1)); % Create a mesh of X and Y coordinates needed to create the mask (mesh should be equal in size to the image!)
    mask = inpolygon(X,Y,xCircle,yCircle);          % Check for each x and y in the mesh if it is inside or outside the 'object' on the mask.
    
    Img_fft = fft2(Img);                            % Perform a 2D-fast fourier transform on the image
    mask_fft = fft2(mask);                          % Perform a 2D-fast fourier transform on the mask
    
    Img_conv_fft = Img_fft.*mask_fft;               % Calculate the FT convolution of the image by multiplying the FT image and FT mask
    Img_conv = ifft2(Img_conv_fft);                 % Obtain the convoluted image by performing an 2D-inverse FFT
    Img_conv_norm = Img_conv./max(max(Img_conv));   % Normalize the image such that the maximum value will always be 1 (background removal in pretreatment makes sure that the lowest value is always 0)
    Img_conv_bin = imbinarize(Img_conv_norm,0.5);   % Binarize the image, such that it can be used in regionprops
      
    Pparticles = ParticleSelection(Img_conv_bin, r, SelectionCriteria); % Determine if the found patches are particles according to the SelectionCriteria
end
