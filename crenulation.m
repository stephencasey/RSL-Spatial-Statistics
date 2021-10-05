% Analysis to evaluate oversampling, mapping resolution, and patch
% perimeter scaling
%
% Reducing the resolution of a 2-D simple object (like a square) will
% result in a proportional decrease in the relative perimeter of that
% object. In other words, downsampling a square of width = 2 to a square
% of width = 1, changes the relative perimeter by 1/2. We expect this
% linear relationship with oversampled objects too: the reduction in scale
% is roughly equal to the perimeter reduction (across the scale of
% oversampling). When data is not oversampled (or we are in the range
% beyond the oversampling), the reduction in perimeter may be greater as
% one loses fine scale crenulation due to downsampling.
%
% This script downsamples maps and plots the loss of perimeter as a
% function of that downsampling. If a map is oversampled, we expect an
% equal reduction in perimeter across the oversampled range (a horizontal
% line). As one reaches the actual resolution of the map, the relationship
% will diverge from this horizontal line as detail is lost (assuming the
% map contains detail and is not simply rectangles or simple shapes). This
% can be helpful in determining if a map is oversampled.
%
% Likewise, by evaluating the magnitude of perimeter loss at different
% scales beyond the apparent oversampling range, one may be able to infer
% the scales at which crenulation is most (or least) present and thereby 
% the landscape perimeter scaling.
%
% The perimeter difference is calulated using two methods. Method 1 uses
% the actual difference in perimeter between each downsampling step, while
% Method 2 uses the relative difference at each step to the original
% unmodified perimeter.
%
% The input is a series of binary maps with a common base filename followed
% by a digit saved as a .mat file (e.g. examplemap3.mat). The output is a
% matrix containing the reduction in perimeter values using Method 2, as
% well as plots showing the average reduction in perimeter for both
% methods.


clear
clc

% Directory containing maps with beginning of filename
filedirectory = 'D:\Backup\RSL Backup\RSL Patterning Analysis\Matlab Analysis\ENPmap_analysis\Maps\ENPmapall';
scalechange=1/1.3;      % Reduction in scale for each step
scalesteps=20;          % Total number of downsampling steps
nmaps = 16;              % Number of maps in directory

scalechange1=scalechange.^(1:scalesteps)';
for mapn=1:nmaps
filename=[filedirectory int2str(mapn)];
load(filename);
z=1-(data==0);
z=logical(z);

disp(mapn)
% Uncomment to display an image of each map as it loops
% imshow(1-z)
for n=1:scalesteps
    z1=imresize(z,scalechange1(n));
    c=bwlabel(z1,4);
    rstats=regionprops(c,'Perimeter');
    perims=cat(1,rstats.Perimeter);
    perimall(n,mapn)=sum(perims);
end

% Method 1
perimchange1(:,mapn)=(perimall(2:end,mapn)./perimall(1:end-1,mapn));
% Method 2
perimchange2(:,mapn)=(perimall(2:end,mapn)./perimall(1,mapn))./(scalechange1(2:end)./(scalechange1(1)^2));    

% % Uncomment to plot individual map perim change
% figure(mapn)
% plot((perimchange2(:,mapn)),'x')
% semilogx(1./scalechange1(1:end-1),(perimchange2(:,mapn)),'x')
% xlabel('Pixel Size')
% ylabel('Reduction in Perimeter')
% hold on
% plot([1 1000],[scalechange scalechange])
% hold off

end
disp(perimchange2)

% Plot mean perim change
figure(mapn+1)
plot(mean(perimchange2,2),'x')
semilogx(1./scalechange1(1:end-1),mean(perimchange2,2),'x')
title('Mean Values (Method 2)')
xlabel('Pixel Size')
ylabel('Reduction in Perimeter')
hold on
plot([1 1000],[scalechange scalechange])
hold off

figure(mapn+2)
plot(mean(perimchange1,2),'x')
semilogx(1./scalechange1(1:end-1),mean(perimchange1,2),'x')
title('Mean Values (Method 1)')
xlabel('Pixel Size')
ylabel('Reduction in Perimeter')
hold on
plot([1 1000],[scalechange scalechange])
hold off