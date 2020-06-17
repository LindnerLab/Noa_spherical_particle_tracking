function [Pparticles] = FindParticlesConvolution(Img,r,threshold)
    Img_size = size(Img);
%     Img = Img - threshold;
%     Img(Img < 0) = 0;
    
    xCenter = r;
    yCenter = r;
    xCircle = r*cos(0:pi/15:2*pi) + xCenter;
    yCircle = r*sin(0:pi/15:2*pi) + yCenter;
    
    [X, Y] = meshgrid(1:Img_size(2),1:Img_size(1));
    mask = inpolygon(X,Y,xCircle,yCircle);
    
    Img_fft = fft2(Img);
    mask_fft = fft2(mask);
    
    Img_conv_fft = Img_fft.*mask_fft;
    Img_conv = ifft2(Img_conv_fft);
    Img_conv_norm = Img_conv./max(max(Img_conv));
    Img_conv_bin = imbinarize(Img_conv_norm,0.5);
    
    Pcentroid = regionprops(Img_conv_bin,'Centroid');
    Parea = regionprops(Img_conv_bin,'Area');
    Peccentricity = regionprops(Img_conv_bin,'Eccentricity');
    particle = Pcentroid(logical([[Parea.Area] > pi*(r-2)^2] .* [[Parea.Area] < pi*(r+2)^2] .* [[Peccentricity.Eccentricity] < 0.4]));
    
    if ~isempty(particle)
        for i = 1:length(particle)
            Pparticles(i,:) = [particle(i).Centroid];
        end
        
        Pparticles = Pparticles-r;
    else
        Pparticles = [];
    end
end

