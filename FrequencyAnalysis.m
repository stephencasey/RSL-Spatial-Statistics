% Periodogram (Radial & Angular Spectrum) Analysis for Images

clear
clc

map_n=33;
analysisname='1995';
filedirectory=['C:\Users\thorn\OneDrive\Desktop\RSL_Github\RSL\RSL\' analysisname];

allangles=1;    % Omnidirectional R-Spectrum = 1; Horizontal R-Spectrum = 0
                % A-Specrum is only computed for omnidirectional case

orientation=zeros(map_n,1);
c=zeros(map_n,1);

excelfilename= [filedirectory '\FourierTransform\rspectrum_test.xlsx'];
delete(excelfilename)

for mapn=2:map_n
clearvars -except excelfilename c mapn filedirectory map_n orientation allangles
filename=[filedirectory '\Maps\1995all' int2str(mapn) '.mat'];
load(filename)

data=ismember(data,[1 2 3]);

data(1,:)=[];
data(:,1)=[];

smalldim=min(size(data));
data=data(1:smalldim,1:smalldim);

b_length=length(data);
z=b_length/2;
in=(0:1:z);
[x,y]=meshgrid(in,in);
distance=sqrt(x.^2+y.^2);
sort_distance=sort(distance(:));
distance(distance>z-1)=z-1;
distance=[fliplr(distance) distance(:,2:end)];
distance=[flipud(distance);distance(2:end,:)];

% Make extent circular for A-Spectrum
if allangles==1
    distance2=distance-(z-2);
    distance2(distance2<0)=0;
    data(data==1)=data(data==1)-.5*distance2(data==1)./max(distance2(:));
    data(data==0)=data(data==0)+.5*distance2(data==0)./max(distance2(:));
end

cell_angle=rot90(atan(y./x).*180/pi);
cell_angle=[cell_angle;90+rot90(cell_angle(:,2:end),3)];
cell_angle=[rot90(cell_angle(:,2:end),2) cell_angle];

cell_angle(distance>=b_length/2)=nan;

b=fftshift(fft2(data));
b_log=log(abs(b));

n1=1;
n=1;
r_spectrum=zeros((b_length+1)/2,1);
dis=zeros((b_length+1)/2,1);
if allangles == 0
    for n2=2:1:((b_length+1)/2+1)
        r_spectrum(n)=mean((b_log(distance>=n1 & distance < n2 & cell_angle>=80 & cell_angle<100)));
        dis(n)=mean([n1,n2]);
        n=n+1;
        n1=n2;
    end;
else
    for n2=2:1:((b_length+1)/2+1)
        r_spectrum(n)=mean((b_log(distance>=n1 & distance < n2)));
        dis(n)=mean([n1,n2]);
        n=n+1;
        n1=n2;
    end;
end

c(mapn)=corr(dis,r_spectrum,'type','Spearman');

p=figure(1);
xaxis=1:1:(b_length+1)/2;
plot(xaxis,r_spectrum)
title('R-Spectrum')
xlabel('Wavenumber')
ylabel('Power')
filename= [filedirectory '\FourierTransform\' 'rspectrum' int2str(mapn) '.png'];
saveas(p,filename,'png')

im1=figure(2);
imshow(b_log,[]), colormap(jet(64)) 
axis on

filename= [filedirectory '\FourierTransform\' 'periodogram' int2str(mapn) '.png'];
saveas(im1,filename,'png')

excel_index = [idx2A1(mapn) '2'];
xlswrite(excelfilename,r_spectrum,'r_spectrum',excel_index);
disp(mapn)

% A-Spectrum
if allangles == 1
    a_spectrum=zeros(180,1);
    a_spectrum(1)=mean(log(abs(b(cell_angle>=179.5 | cell_angle < .5))));
    n1=.5;
    for n=2:180
        a_spectrum(n)=mean(log(abs(b(cell_angle>=n1 & cell_angle < n1+1))));
        n1=n1+1;
    end
    figure(3)
    plot(0:1:359,[a_spectrum;a_spectrum])
    title('S-Spectrum')
    xlabel('Angle (degrees)')
    ylabel('Amplitude')
    % Max A-Spectrum indicates orientation of landscape
    orientation(mapn)=find(max(a_spectrum)==a_spectrum)-1;
end

figure(4)
imshow(1-data)
drawnow
end

