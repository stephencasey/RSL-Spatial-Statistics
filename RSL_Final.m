% RSL_Final calculates patterning metrics from input maps & prints them to an excel file
%
% This script takes rasterized maps of Everglades ridge-slough landscapes, 
% converts them to a binary representation and performs numerous
% statistical analyses on the spatial characteristics.
%
% If running the map making section of this script, the inputs are a series
% of raw rasterized maps with multiple vegetation classes as a .txt file,
% or if not running the map making section, the inputs are binary maps
% saved as .mat files
%
% Maps provided with this script were produced by the South Florida Water 
% Management District (Rutchey et al. 2005)
%
% This is the primary script used for the analysis presented in Casey et.
% al. 2015, however distribution testing and frequency analysis are both
% accomplished in separate scripts.
%
% References:
% Casey, Stephen T., et al. "Hydrologic controls on aperiodic spatial 
% organization of the ridge–slough patterned landscape." Hydrology and 
% Earth System Sciences 20.11 (2016): 4457-4467.
%
% Rutchey, K., Vilchek, L., and Love, M.: Development of a vegetation map for Water Conservation Area 3, Technical Publication
% ERA Number 421, South Florida Water Management District,
% West Palm Beach, FL, USA, 2005.


%%  Map Making
% 
clear
clc

% Input Values

mapy=6000;          % vertical extent
mapx=6000;          % horizontal extent
pixelsize=1;       % size of pixels (m^2)
xmin=1;             % delete lengths less than this (m)
regressionmin=1;    % minimum number of pixels for patch regression & PAFRAC
TI_threshold=10000; % lower threshold for tree island classification (m^2)
mapn=33;            % total number of maps
row_n=33;           % number of map rows (set equal to number of maps)
paretocutoff=10;    % lower cutoff for pareto fitting (m^2)
scale=1;            % Scale to resample map to
rotatemap90=0;      % Rotate map 90 degrees?
makemap=0;          % Run entire map making section or use saved maps (1=yes 0=no)

pixelsize=pixelsize/scale;
% Sets the file directory
analysisname='1995';
filedirectory=['C:\Users\thorn\OneDrive\Desktop\RSL_Github\RSL\RSL\' analysisname '\'];
% Set the base filename for the new maps
mapname='1995all';   

mapy=mapy/pixelsize;
mapx=mapx/pixelsize;

% Creates output folders
mkdir([filedirectory 'MapImages'])
mkdir([filedirectory 'Elongation&Isotropy'])
mkdir([filedirectory 'OutputData'])
mkdir([filedirectory 'CDFs'])
mkdir([filedirectory 'FourierTransform'])
% Converts from vegetation classes to 0=slough, 1=ridge, 2=tree island, 3=other 
if makemap == 1
    map_n=1;
    for row=1:row_n
        
        nodatafix=1;    % turn on/off nodata fix 
         fprintf('Map %d \n', row);
        clearvars -except mapx mapy subodh4000_index runs reps rotatemap90 scale cattailall nodatafix paretocutoff mapname filedirectory row mapy mapx mapn others TI_threshold pixelsize map_n row_n PAmin mapsize xmin regressionmin
        filename=[filedirectory 'Maps\' int2str(row) '.txt'];
        
        if exist(filename,'file')==2
            
        data_all=load(filename);
        
        if rotatemap90 == 1
            data_all=rot90(data_all);
        end
        ranges=size(data_all);

        extra(1)=rem(ranges(1),mapy);
        extra(2)=rem(ranges(2),mapx);

        if extra(1)>0
            n=extra(1);
            while n>0
                data_all(end,:)=[];
                n=n-1;
                if n<1
                    break
                end
                data_all(1,:)=[];
                n=n-1;
            end
        end

        if extra(2)>0
            n=extra(2);
            while n>0
                data_all(:,end)=[];
                n=n-1;
                if n<1
                    break
                end
                data_all(:,1)=[];
                n=n-1;
            end
        end

        ranges=size(data_all);
        nmaps=ranges(2)/mapx;

        %fill in blank NODATA cells using the adjacent neighborhood
        fprintf('%d NODATA cells \n', sum(sum(data_all==-9999)));
        while nodatafix==1
            data_all(data_all==-9999)=9999;
            data_all=padarray(data_all,[1,1],-1);
            [nodata_i,nodata_j]=find(data_all==9999);
            for n=1:length(nodata_j)
                nodata_window=data_all(nodata_i(n)-1:nodata_i(n)+1,nodata_j(n)-1:nodata_j(n)+1);
                nodata_window(nodata_window==-1)=[];
                mode_nodata=histc(nodata_window(:),unique(nodata_window));
                data_all(nodata_i(n),nodata_j(n))=mode(nodata_window(:));
            end

            data_all(1,:)=[];
            data_all(:,1)=[];
            data_all(:,mapx+1)=[];
            data_all(mapy+1,:)=[];
            if sum(sum(data_all==9999))==0
                nodatafix=0;
            end
        end

    % biomeconvert is a list of all the keys that determine what vegetation
    % class is converted to ridge, slough, tree island and other
    
    %     % No Conversion
    %     data_all=data_all+1;
    %     biomeconvert=[0 1 2 3];
    
    %     1995 converted maps
    %     biomeconvert=[3 3 3 1 3 0 1 1 2 2 2 2 1 2 2 1 2 3 3 1 0 3 0 1 2 2 1 1 3 2 0 3 3 2 2 1 2 3 0 2 2 2 2 1 2 3 2 3 3 3 3 3 3];

        % 1995final maps
          biomeconvert=[2 1 2 2 2 3 0 3 2 4 0 3 0 2 2 2 2 0 0 0 0 2 2 0 1 2 1 1 0 1 2 2 2 3];

        % 2004 maps
     %    biomeconvert=[3 3 3 1 3 0 1 1 2 2 2 2 1 2 2 1 2 3 3 1 0 3 0 1 2 2 1 1 3 2 0 3 3 2 2 1 2 3 0 2 2 2 2 1 2 3 2 3 3 3 3 3 3];

    %   % Nungesser maps
    %     biomeconvert=[1 0];

    %   % ENP maps
    %     biomeconvert=[0 1 2 2 1 2 3 0 2 0 3 1 0 3 0 3 2 3 3 3 3 3 3 3 3 3 2 2 3 3 3 3 3 3 3 3 3 3 3 2 2 3 2 3 3 3 3 3];

        rowdata=zeros(size(data_all));
        for k=1:numel(data_all)
            rowdata(k)=biomeconvert(data_all(k));
        end
        cattailall=zeros(1,nmaps);
        others=zeros(1,nmaps);
        for n=1:nmaps
            data=rowdata(:,(n*mapx-mapx)+1:n*mapx);
            cattailall(map_n)=sum(sum(data==4))/(mapx*mapy);
            data(data==4)=1;
            
            if scale ~= 1
            data1=data==1;
            data2=data==2;
            data3=data==3;

            data1=imresize(data1,scale);
            data2=imresize(data2,scale);
            data3=imresize(data3,scale);

            datanew=zeros(size(data1));
            datanew(data1)=1;
            datanew(data2)=2;
            datanew(data3)=3;

            data=datanew;
            end

            filename=[filedirectory 'Maps\' mapname int2str(map_n) '.mat'];
            save(filename,'data');

            filename=[filedirectory mapname int2str(map_n) '.mat'];
            save(filename,'data');

            others(map_n)=sum(sum(data==3));
            colormap_data=[0 0 0; 0 0 0; 0 0 0; 0 0 0];
            map=label2rgb(data,colormap_data);

            filename=[filedirectory '\MapImages\BW' int2str(map_n) '.png'];
            delete(filename)
            imwrite(map,filename,'png');

            map_n=map_n+1;
        end
        else
        end   
    end
    mapn=map_n-1;
    row_n=row-1;

    % Create maps for analysis with Fragstats

    for site=1:mapn
        clearvars data
    filename= [filedirectory 'Maps\' mapname int2str(site)];

    load(filename)

    data=data+1;

    data(data==3)=2;
    data(data==4)=2;

    filename= [filedirectory 'Maps\ascmap' int2str(site) '.asc'];
    dlmwrite(filename, data, 'delimiter', ' ')
    end



    %% Converts to Tree Island Classification (0=slough, 1=ridge, 2=TIaff, 3=TI, 4=Other)
    %
    % This is used to visualize tree islands and their patch connectivity
    % (see TI maps in MapImages folder; red=TreeIsland, green=DirectConnections)
    
    for map_n=1:mapn
        clearvars -except cattailall paretocutoff mapname filedirectory map_n mapn mapy mapx TI_threshold pixelsize mapsize xmin PAmin
        filename=[filedirectory 'Maps\' mapname int2str(map_n)];
        load(filename)
        
        neighborh=4;
        
        data_treeisland=data==2;
        data_r=ismember(data,[1 2]);
        
        c_ti=bwlabel(data_treeisland,neighborh);
        ti_areas=regionprops(c_ti,'area');
        ti_areas=cat(1,ti_areas.Area);
        
        clearvars data_treeisland
        
        threshold_index=find(ti_areas>=TI_threshold/(pixelsize^2));
        data_threshold=ismember(c_ti,threshold_index);
        
        clearvars c_ti
        
        data_ti=zeros(mapy,mapx);
        data_ti=data_ti~=0;
        for i=1:mapy
            zeta=bwlabel(data_r(i,:));
            [B first1]=unique(zeta,'first');
            first1(1)=[];
            [B last1]=unique(zeta,'last');
            last1(1)=[];
            diff1=last1-first1;
            for n=1:length(first1)
                if sum(data_threshold(i,first1(n):last1(n))==1)>0
                    data_ti(i,first1(n):last1(n))=1;
                end;
            end;
        end;
        
        for j=1:mapx
            zeta=bwlabel(data_r(:,j))';
            [B first1]=unique(zeta,'first');
            first1(1)=[];
            [B last1]=unique(zeta,'last');
            last1(1)=[];
            diff1=last1-first1;
            for n=1:length(first1)
                if sum(data_threshold(first1(n):last1(n),j)==1)>0
                    data_ti(first1(n):last1(n),j)=1;
                end;
            end;
        end;
        
        data_threshold=data_ti;
        data_ti=zeros(mapy,mapx);
        data_ti=data_ti~=0;
        for i=1:mapy
            zeta=bwlabel(data_r(i,:));
            [B first1]=unique(zeta,'first');
            first1(1)=[];
            [B last1]=unique(zeta,'last');
            last1(1)=[];
            diff1=last1-first1;
            for n=1:length(first1)
                if sum(data_threshold(i,first1(n):last1(n))==1)>0
                    data_ti(i,first1(n):last1(n))=1;
                end;
            end;
        end;
        
        for j=1:mapx
            zeta=bwlabel(data_r(:,j))';
            [B first1]=unique(zeta,'first');
            first1(1)=[];
            [B last1]=unique(zeta,'last');
            last1(1)=[];
            diff1=last1-first1;
            for n=1:length(first1)
                if sum(data_threshold(first1(n):last1(n),j)==1)>0
                    data_ti(first1(n):last1(n),j)=1;
                end;
            end;
        end;
        
        clearvars -except cattailall paretocutoff mapname filedirectory data_ti map_n PAmin tstart mapn mapy mapx TI_threshold pixelsize map_size xmin regressionmin
        
        filename=[filedirectory 'Maps\' mapname int2str(map_n)];
        load(filename)
        
        new_data=data;
        new_data(data_ti==1)=2;
        new_data(data==2)=3;
        new_data(data==3)=4;
        data=new_data;
        
        colormap_data=[0 0 0; 0 1 0; 1 0 0; 0 0 1];
        map=label2rgb(data,colormap_data);
        
        filename=[filedirectory 'MapImages\TI' int2str(map_n) '.png'];
        delete(filename)
        imwrite(map,filename,'png');
        
        filename=[filedirectory 'Maps\' mapname int2str(map_n)];
        save(filename,'data')
        
    end
    
end

%% RSL Lengths and Widths Analysis


for site=1:mapn
       
    cattail=0;     
    

    clearvars -except regressionmin mapy mapx cattail cattailall mapname paretocutoff filedirectory site tstart mapn xmin pixelsize mapsize PAmin PAexp2 regressionmin
 
    filename=[filedirectory 'Maps\' mapname int2str(site)];
    load(filename)
    
    % set 0 = slough; 1,2,3 = ridge
    data=ismember(data,[1 2 3]);
    
        colormap_data=[0 0 0; 0 0 0; 0 0 0; 0 0 0];
        map=label2rgb(data,colormap_data);
        
        filename=[filedirectory '\MapImages\BW' int2str(site) '.png'];
        
        imwrite(map,filename,'png');
    
    % Downsampler
    % data=imresize(data,1/10,'box');
    
    % Calculate proportion of landscape cover for each class
    density=zeros(1,5);                             
    density(1) = sum(sum(data==0))/numel(data);
    density(2) = sum(sum(data==1))/numel(data);
    density(3) = sum(sum(data==2))/numel(data);
    density(4) = sum(sum(data==3))/numel(data);
    density(5) = sum(sum(data==4))/numel(data);
    
    neighborh=4;    
    
    ranges = size(data);
    rangex=ranges(2);
    rangey=ranges(1);
    
    %% calculates slough widths & lengths
    
    data_s=data==0;
    
    z=1;
    first=zeros(round(rangex*rangey/2),1)';
    last=zeros(round(rangex*rangey/2),1)';
    for i=1:rangey
        zeta=bwlabel(data_s(i,:));
        [~, first1]=unique(zeta,'first');
        first1(1)=[];
        [~, last1]=unique(zeta,'last');
        last1(1)=[];
        
        first1l=length(first1);
        first(z:(z+first1l-1))=first1;
        last(z:(z+first1l-1))=last1;
        z=z+first1l;
    end;
    
    first(first==0)=[];
    last(last==0)=[];
    
    % wslough and lslough are a list of the continuous length & width transects
    wslough=last-first+1;
    
    z=1;
    first=zeros(round(rangex*rangey/2),1)';
    last=zeros(round(rangex*rangey/2),1)';
    for j=1:rangex
        zeta=bwlabel(data_s(:,j))';
        [~, first1]=unique(zeta,'first');
        first1(1)=[];
        [~, last1]=unique(zeta,'last');
        last1(1)=[];
        
        first1l=length(first1);
        first(z:(z+first1l-1))=first1;
        last(z:(z+first1l-1))=last1;
        
        z=z+first1l;
    end;
    
    first(first==0)=[];
    last(last==0)=[];
    
    lslough=last-first+1;
    
    fprintf('slough widths and lengths \n')

    clearvars data_s
    
    %% Calculates ridge widths and lengths (ungrouped)
    
    data_ti=ismember(data,[2 3]);
    data_r=ismember(data,[1 2 3]);
    c=bwlabel(data_r,neighborh);
    data_slough=logical(1-data_r);
    cslough=bwlabel(data_slough,neighborh);
   
    p = label2rgb(c, 'hsv', 'k','shuffle');
    filename= [filedirectory 'MapImages\Grouped' int2str(site) '.png'];
    imwrite(p,filename,'png')   
    
    ridges=max(max(c));
    
    ridge_map=zeros(rangey,rangex);
    z=1;
    first=0;
    last=0;
    wridge_ti=0;
    for i=1:rangey
        zeta=bwlabel(data_r(i,:));
        zeta_c=c(i,:);
        [~, first1]=unique(zeta,'first');
        first1(1)=[];
        [~, last1]=unique(zeta,'last');
        last1(1)=[];
        
        wridgecenters=(first1+last1)/2;
        wnn1=wridgecenters(2:end)-wridgecenters(1:end-1);   
        
        wnn1l=length(wnn1);
        wnn(z:(z+wnn1l-1))=wnn1;
        
        first1l=length(first1);
        first(z:(z+first1l-1))=first1;
        last(z:(z+first1l-1))=last1;
        wridge_n(z:(z+first1l-1))=zeta_c(first1);
        
        for n=1:first1l
            if min(data_ti(i,first1(n):last1(n)))==1
                wridge_ti(z)=1;
            end;
            z=z+1;
        end;
    end;
    
    first((first==0))=[];
    last((last==0))=[];
    wridge_n((wridge_n==0))=[];
    zeta=zeros(size(wridge_n));
    zeta(1:length(wridge_ti))=wridge_ti;
    wridge_ti=zeta;
    wridge=last-first+1;
    
    if sum(wnn<0)>1
        error('wnn')
    end
    
    z=1;
    first=0;
    last=0;
    lridge_ti=0;
    for j=1:rangex
        zeta=bwlabel(data_r(:,j))';
        error
        zeta_c=c(:,j);
        [~, first1]=unique(zeta,'first');
        first1(1)=[];
        [~, last1]=unique(zeta,'last');
        last1(1)=[];
        
        %     for g=1:length(first1)
        %     if last1(g)-first1(g) > 2000
        %         ridge_map(first1(g):last1(g),j)=2;
        %     end
        %     end
        lridgecenters=(first1+last1)/2;
        lnn1=lridgecenters(2:end)-lridgecenters(1:end-1);   
        
        lnn1l=length(lnn1);
        lnn(z:(z+lnn1l-1))=lnn1;
        
        first1l=length(first1);
        first(z:(z+first1l-1))=first1;
        last(z:(z+first1l-1))=last1;
        lridge_n(z:(z+first1l-1))=zeta_c(first1);
        
        for n=1:first1l
            if min(data_ti(first1(n):last1(n),j))==1
                lridge_ti(z)=1;
            end;
            z=z+1;
        end;
    end;
    
    first((first==0))=[];
    last((last==0))=[];
    lridge_n((lridge_n==0))=[];
    zeta=zeros(size(lridge_n));
    zeta(1:length(lridge_ti))=lridge_ti;
    lridge_ti=zeta;
    
    lridge=last-first+1;
    
    
    wnn(wnn==0)=[];
    lnn(lnn==0)=[];
    
    clearvars data_ti
    
    fprintf('ridge ungrouped widths and lengths \n')
    
    
    %% Calculates grouped stats
    
    ridgestats=zeros(ridges,9);
    rstats=regionprops(c,'all');
    rstatso=regionprops(c','Orientation');
    
    rstatsslough=regionprops(cslough,'All');
    sloughareas=cat(1,rstatsslough.Area);
    sloughperims=cat(1,rstatsslough.Perimeter);
    sloughperims(sloughperims==0)=1;
    
    ridgestats(:,1)=cat(1,rstats.Area);
    ridgestats(:,2)=cat(1,rstats.MajorAxisLength);
    ridgestats(:,3)=cat(1,rstats.MinorAxisLength);
    ridgestats(:,4)=cat(1,rstats.Eccentricity);
    ridgestats(:,5)=cat(1,rstatso.Orientation);
    ridgestats(:,6)=cat(1,rstats.EulerNumber);
    ridgestats(:,7)=cat(1,rstats.Solidity);
    ridgestats(:,8)=cat(1,rstats.Extent);
%     ridgestats(:,9)=cat(1,rstats.Perimeter);
     
%     ridgeperims=ridgestats(:,9);
%     ridgeperims(ridgeperims==0)=1;


        data_pad=padarray(data_r,[1 1]);    
        cellperim=zeros(size(data_pad));
        for j=2:rangey+1 
            for k=2:rangex+1
                if data_pad(j,k)==1
                    cellperim(j,k)=sum([data_pad(j-1,k) data_pad(j+1,k) data_pad(j,k+1) data_pad(j,k-1)]==0);
                end
            end
        end
        cellperim=cellperim(2:rangey+1,2:rangex+1);
        perims=sum(sum(cellperim));
        
        perimindex_p=cellperim(cellperim~=0);
        perimindex_c= c(cellperim~=0);
        
        ridgeperims=zeros(ridges,1);
        for n1=1:ridges
            ridgeperims(n1)=sum(perimindex_p(perimindex_c==n1));
        end
        
        ridgeperims=ridgeperims.*pixelsize;

    grouped=zeros(ridges,23);
    grouped(:,1)=ridgestats(:,1);
    grouped(:,2)=ridgeperims;
    grouped(:,21)=ridgestats(:,4);
    grouped(:,22)=ridgestats(:,5);
    grouped(:,23)=ridgestats(:,7);
    
    grouped_ti=zeros(ridges,1);
    ti_ridges=unique(lridge_ti.*lridge_n);
    ti_ridges(1)=[];
    grouped_ti(ti_ridges)=1;
    
    
    fprintf('calculated orientation metrics\n')
    
    
    %% Top half vs bottom half
    
    bottomhalf=zeros(ridges,1);
    tophalf=zeros(ridges,1);
    for nr=1:ridges
        patchimage=rstats(nr,1).Image;
        pixellist=rstats(nr,1).PixelList;
        ulcorner=[min(pixellist(:,1)) min(pixellist(:,2))];
        centroid=rstats(nr,1).Centroid-(ulcorner-1);
        
        cy=centroid(2)-.5;
        
        
        psize=size(patchimage);
        
        
        patch_pad=padarray(patchimage,[1 1]);    
        cellperim=zeros(size(patchimage));
        for j=2:psize(1)+1 
            for k=2:psize(2)+1
                if patch_pad(j,k)==1
                    cellperim(j,k)=sum([patch_pad(j-1,k) patch_pad(j+1,k) patch_pad(j,k+1) patch_pad(j,k-1)]==0);
                end
            end
        end
        
        cellperim=cellperim(2:psize(1)+1,2:psize(2)+1);
        
        if rem(round(psize(1)),2)~=0
            tophalf(nr)=sum(sum(cellperim(1:(round(psize(1))+1)/2,:)));
            bottomhalf(nr)=sum(sum(cellperim(end-(round(psize(1))+1)/2+1:end,:)));
        else
            tophalf(nr)=sum(sum(cellperim(1:round(psize(1))/2,:)));
            bottomhalf(nr)=sum(sum(cellperim(end-round(psize(1))/2+1:end,:)));
        end
        
        
        
    end
    
    toverb=sum(tophalf)/sum(bottomhalf);
    toverb50000=sum(tophalf(grouped(:,1)>50000))/sum(bottomhalf(grouped(:,1)>50000));
    toverb10000=sum(tophalf(grouped(:,1)>10000))/sum(bottomhalf(grouped(:,1)>10000));
   
    tobPatch=mean(log(tophalf./bottomhalf));
    tobPatch10000=mean(log(tophalf(grouped(:,1)>10000)./bottomhalf(grouped(:,1)>10000)));
    tobPatch50000=mean(log(tophalf(grouped(:,1)>50000)./bottomhalf(grouped(:,1)>50000)));
    tobPatch500000=mean(log(tophalf(grouped(:,1)>500000)./bottomhalf(grouped(:,1)>500000)));
    
   
    
    [h,p10000]=ttest(log(tophalf(grouped(:,1)>10000)./bottomhalf(grouped(:,1)>10000)));
    
    excelfilename=[filedirectory 'topvsbottom.xlsx'];
    if site==1
        delete(excelfilename)
    end
    
    header={'Site' 'toverb' '>50000' '>10000' 'Patch' 'Patch10000' 'Patch50000' 'Patch500000' 'P10000'};
    xlswrite(excelfilename,header,'toverb','A1');
    excel_index = ['A' int2str(1+site)];
    xlswrite(excelfilename,[site toverb toverb50000 toverb10000 tobPatch tobPatch10000 tobPatch50000 tobPatch500000 p10000],'toverb',excel_index);  
     
   
    
    
    %% Grouped Lengths and widthds
    
    for n=1:ridges
        
        
        pixely=rstats(n).PixelList(:,2);
        pixelx=rstats(n).PixelList(:,1);
        ridgepic=rstats(n).Image;
        
        size_n=size(ridgepic);
        
        
        %Sum of Cells Method
        if min(size(ridgepic)) > 1
            l1=sum(ridgepic);
            w1=sum(ridgepic');
        else
            l1=size_n(1);
            w1=size_n(2);
        end;
        
        
        %Continuous Method
        l3=lridge((lridge_n == n));
        w3=wridge((wridge_n == n));
        
        
        grouped(n,3)=length(l1)';   %poly W
        grouped(n,4)=length(w1)';   %poly L
        grouped(n,5)=mean(w1)';   %Average sum based W
        grouped(n,6)=mean(l1)';   %Average sum based L
        grouped(n,7)=median(w1)';
        grouped(n,8)=median(l1)';
        grouped(n,9)=max(w1)';
        grouped(n,10)=max(l1)';
        grouped(n,11)=var(w1)';
        grouped(n,12)=var(l1)';
        grouped(n,13)=mean(w3)';   %Average continuous transect based W
        grouped(n,14)=mean(l3)';   %Average continuous transect based L
        grouped(n,15)=median(w3)';
        grouped(n,16)=median(l3)';
        grouped(n,17)=max(w3)';
        grouped(n,18)=max(l3)';
        grouped(n,19)=var(w3)';
        grouped(n,20)=var(l3)';
        
        
        
    end;
    
    fprintf('calculate grouped widths and lengths \n')
    
    
    %% Convert pixels to proper units
    
    wslough=wslough*pixelsize;
    lslough=lslough*pixelsize;
    wridge=wridge*pixelsize;
    lridge=lridge*pixelsize;
    grouped(:,1)=grouped(:,1)*pixelsize^2;
    grouped(:,2:20)=grouped(:,2:20)*pixelsize;
    wnn=wnn*pixelsize;
    lnn=lnn*pixelsize;
    sloughareas=sloughareas*pixelsize^2;
    sloughperims=sloughperims*pixelsize;
    
    
    
    %% Delete under map res
    
    wslough=wslough(wslough >= xmin);
    lslough=lslough(lslough >= xmin);
    
    
    wridge_ti=wridge_ti(wridge >= xmin);
    lridge_ti=lridge_ti(lridge >= xmin);
    wridge=wridge(wridge >= xmin);
    lridge=lridge(lridge >= xmin);
    
    ridgestats(grouped(:,1)<xmin^2,:)=[];
    grouped_ti(grouped(:,1)<xmin^2,:)=[];
    grouped(grouped(:,1)<xmin^2,:)=[];
    
    ridges=length(grouped(:,1));
    ti_ridges=length(grouped_ti(:,1));
    
    sloughperims(sloughareas<xmin^2)=[];
    sloughareas(sloughareas<xmin^2)=[];
    
    %% Split up values belonging to tree islands
    
    wridge_ok=wridge((wridge_ti==0));
    wridge_tree=wridge((wridge_ti==1));
    lridge_ok=lridge((lridge_ti==0));
    lridge_tree=lridge((lridge_ti==1));
    
    %% Construct output arrays
    
    a_wslough(:,1)=wslough';
    a_lslough(:,1)=lslough';
    a_wridge(:,1)=wridge';
    a_lridge(:,1)=lridge';
    a_wridge_ok(:,1)=wridge_ok';
    a_lridge_ok(:,1)=lridge_ok';
    a_wridge_tree(:,1)=wridge_tree'; 
    a_lridge_tree(:,1)=lridge_tree';
    
    wnn=wnn';
    lnn=lnn';
    
    fprintf('construct output arrays \n')
    
    
    
    %% Create Plots
    
    orientation = mean(ridgestats(:,5));
    
    areas100=grouped(:,1);
    areas100=areas100(areas100>=100);
 
    if length(areas100) >= 100
        
        areas=grouped(:,1);
        log_areas=log10(areas);
        
        fit_func=@(a,x) a(1)+a(2)*x;
        [ecc_fit,ecc_r,ecc_j,ecc_covb,ecc_mse]=nlinfit(log_areas,ridgestats(:,4),fit_func,[0 1]);
        fit_func=@(ecc_fit,x) ecc_fit(1)+ecc_fit(2)*x;
        fit_func_line(1)=fit_func(ecc_fit,min(log_areas));
        fit_func_line(2)=fit_func(ecc_fit,(max(log_areas)));
        
        p=figure(1);
        plot(log_areas,ridgestats(:,4),'.');
        xlim([min(log_areas) ceil(max(log_areas))]);
        ylim([min(ridgestats(:,4)) 1]);
        title(['Patch Eccentricity (Site ' num2str(site) ')'])
        xlabel('log(Patch Area)')
        ylabel('Eccentricity')
        hold on
        plot([min(log_areas) max(log_areas)],fit_func_line,'r')
        hold off
        filename= [filedirectory 'Elongation&Isotropy\Eccentricity' int2str(site) '.png'];
        saveas(p,filename,'png')
        
        
        plot(log_areas,ridgestats(:,5),'.');
        hold on
        plot([min(log_areas) max(log_areas)],[orientation orientation],'r');
        hold off
        xlim([min(log_areas) ceil(max(log_areas))]);
        title(['Patch Orientation (Site ' num2str(site) ')'])
        xlabel('log(Patch Area)')
        ylabel('Orientation (degrees)')
        filename= [filedirectory 'Elongation&Isotropy\Orientation' int2str(site) '.png'];
        saveas(p,filename,'png')
        
        
        plot(log_areas,ridgestats(:,7)*100,'.');
        xlim([min(log_areas) ceil(max(log_areas))]);
        title(['Solidity (Site ' num2str(site) ')'])
        xlabel('log(Patch Area)')
        ylabel('Solidity (percent)')
        filename= [filedirectory 'Elongation&Isotropy\Solidity' int2str(site) '.png'];
        saveas(p,filename,'png')
        
        
        lw_grouped = log10(grouped(:,6)./grouped(:,5));
        fit_func=@(a,x) a(1)+a(2)*x;
        [lw_fit,lw_r,lw_j,lw_covb,lw_mse]=nlinfit(log_areas,lw_grouped,fit_func,[0 1]);
        fit_func=@(lw_fit,x) lw_fit(1)+lw_fit(2)*x;
        fit_func_line(1)=fit_func(lw_fit,min(log_areas));
        fit_func_line(2)=fit_func(lw_fit,(max(log_areas)));
        
        LWr=1-(sum(lw_r.^2)/sum((lw_grouped-mean(lw_grouped)).^2));
        
        plot(log_areas,lw_grouped,'.');
        
        
        title(['log(L/W)1 (Site ' num2str(site) ')'])
        xlabel('log(Patch Area)')
        ylabel('log(L/W)')
        xlim([min(log_areas) max(log_areas)]);
        hold on
        plot([min(log_areas) max(log_areas)],fit_func_line,'r')
        text(5.5,0,[num2str(lw_fit(1)),'+',num2str(lw_fit(2)),'*X'])
        text(5.5,-.25,['r^2 = ',num2str(LWr)])
        hold off
        filename= [filedirectory 'Elongation&Isotropy\LW(1)' int2str(site) '.png'];
        saveas(p,filename,'png')
       
        
        
        
        lwscaling=lw_fit(2);
        
        lw_grouped2 = log10(grouped(:,14)./grouped(:,13));
        fit_func=@(a,x) a(1)+a(2)*x;
        [lw_fit2,lw_r,lw_j,lw_covb,lw_mse2]=nlinfit(log_areas,lw_grouped2,fit_func,[0 1]);
        fit_func=@(lw_fit,x) lw_fit(1)+lw_fit(2)*x;
        fit_func_line(1)=fit_func(lw_fit2,min(log_areas));
        fit_func_line(2)=fit_func(lw_fit2,(max(log_areas)));
        
        
        
        
        
        LWr2=1-(sum(lw_r.^2)/sum((lw_grouped2-mean(lw_grouped2)).^2));
        
        plot(log_areas,lw_grouped2,'.');
        
      
        excel_index = [idx2A1(site*2-1) '2'];
        filename= [filedirectory 'Elongation&Isotropy\LWgrouped.xlsx'];
        xlswrite(filename,areas,'LWgrouped',excel_index);
        excel_index = [idx2A1(site*2) '2'];
        xlswrite(filename,lw_grouped2,'LWgrouped',excel_index); 
        
        excel_index = [idx2A1(site*2-1) '2'];
        filename= [filedirectory 'Elongation&Isotropy\Orientation.xlsx'];
        xlswrite(filename,areas,'Orientation',excel_index);
        excel_index = [idx2A1(site*2) '2'];
        xlswrite(filename,ridgestats(:,5),'Orientation',excel_index);
       
        
        
        
        title(['log(L/W)2 (Site ' num2str(site) ')'])
        xlabel('log(Patch Area)')
        ylabel('log(L/W)')
        xlim([min(log_areas) max(log_areas)]);
        hold on
        plot([min(log_areas) max(log_areas)],fit_func_line,'r')
        text(5.5,0,[num2str(lw_fit2(1)),'+',num2str(lw_fit2(2)),'*X'])
        text(5.5,-.25,['r^2 = ',num2str(LWr2)])
        hold off
        filename= [filedirectory 'Elongation&Isotropy\LW(2)' int2str(site) '.png'];
        saveas(p,filename,'png')
        
        
        lw_max = log(grouped(:,18)./grouped(:,17));
        fit_func=@(a,x) a(1)+a(2)*x;
        [lw_fit,lw_r,lw_j,lw_covb,lw_mse]=nlinfit(log_areas,lw_max,fit_func,[0 1]);
        fit_func=@(lw_fit,x) lw_fit(1)+lw_fit(2)*x;
        fit_func_line(1)=fit_func(lw_fit,min(log_areas));
        fit_func_line(2)=fit_func(lw_fit,(max(log_areas)));
        
        LWr=1-(sum(lw_r.^2)/sum((lw_max-mean(lw_max)).^2));
        
        plot(log_areas,lw_max,'.');
        title(['Max log(L/W) (Site ' num2str(site) ')'])
        xlabel('log(Patch Area)')
        ylabel('log(MaxL/MaxW)')
        xlim([min(log_areas) max(log_areas)]);
        hold on
        plot([min(log_areas) max(log_areas)],fit_func_line,'r')
        text(5.5,0,[num2str(lw_fit(1)),'+',num2str(lw_fit(2)),'*X'])
        text(5.5,-.25,['r^2 = ',num2str(LWr)])
        hold off
        filename= [filedirectory 'Elongation&Isotropy\LWmax' int2str(site) '.png'];
        saveas(p,filename,'png')
       
        lwscalingmax=lw_fit(2);
        
        
        
        
    PAmin1=regressionmin*pixelsize^2;
    PAareas=grouped(:,1);
    PAperims=grouped(:,2);
    logPAperims=log10(PAperims(PAareas>=2500));
    logPAareas=log10(PAareas(PAareas>=2500));
    
    fit_func=@(a,x) a(1)+a(2)*x.^a(3);
    [PA_fit,PA_r,PA_j,PA_covb,PA_mse]=nlinfit(logPAperims,logPAareas,fit_func,[0 1 1]);
    fit_func=@(x) PA_fit(1)+PA_fit(2)*x^PA_fit(3);
    PAr=1-(sum(PA_r.^2)/sum((logPAareas-mean(logPAareas)).^2));
    
    fit_func1=@(a,x) a(1)+a(2)*x;
    [PA_fit1,PA_r1,PA_j,PA_covb,PA_mse]=nlinfit(logPAperims,logPAareas,fit_func1,[0 1]);
    fit_func1=@(x) PA_fit1(1)+PA_fit1(2)*x;
    PAr1=1-(sum(PA_r1.^2)/sum((logPAareas-mean(logPAareas)).^2));
    
    PAFRAC1=2/PA_fit1(2);
    PAexp=PA_fit(3);
    PAexp2(site)=PA_fit(3);
    
    logPAperims5=log10(PAperims(PAareas>=10000));
    logPAareas5=log10(PAareas(PAareas>=10000));
    
    fit_func5=@(a,x) a(1)+a(2)*x;
    [PA_fit5,PA_r5,PA_j,PA_covb,PA_mse]=nlinfit(logPAperims5,logPAareas5,fit_func5,[0 1]);
    fit_func5=@(x) PA_fit5(1)+PA_fit5(2)*x;    
    PAFRAC5=2/PA_fit5(2);
    PAr5=1-(sum(PA_r5.^2)/sum((logPAareas5-mean(logPAareas5)).^2));
    
    logPAperims10=log10(PAperims(PAareas>=25000));
    logPAareas10=log10(PAareas(PAareas>=25000));
    
    fit_func10=@(a,x) a(1)+a(2)*x;
    [PA_fit10,PA_r10,PA_j,PA_covb,PA_mse]=nlinfit(logPAperims10,logPAareas10,fit_func10,[0 1]);
    fit_func10=@(x) PA_fit10(1)+PA_fit10(2)*x;
    PAFRAC10=2/PA_fit10(2);
    PAr10=1-(sum(PA_r10.^2)/sum((logPAareas10-mean(logPAareas10)).^2));
    
    p=figure(1);
    plot(logPAperims,logPAareas,'.');
    xlim([min(logPAperims) ceil(max(logPAperims))]);
    ylim([min(logPAareas) ceil(max(logPAareas))]);
    title(['Patch Perimeter VS Areas (Site ' num2str(site) ')'])
    xlabel('log(Perimeter)')
    ylabel('log(Patch Area)')
    hold on
    fplot(fit_func,[min(logPAperims) ceil(max(logPAperims))])
    exponenttext = sprintf('{%d}',PA_fit(3));
    fplot(fit_func5,[min(logPAperims) ceil(max(logPAperims))],'r')
    text(2.25,5.7,['r^2 = ',num2str(PAr5)])
    text(2.25,5.9,['PAFRAC(15) = ',num2str(PAFRAC5)])
    text(2.25,6.2,['r^2 = ',num2str(PAr10)])
    text(2.25,6.4,['PAFRAC(25) = ',num2str(PAFRAC10)])
    text(2.25,6.7,['r^2 = ',num2str(PAr1)])
    text(2.25,6.9,['PAFRAC = ',num2str(PAFRAC1)])
    hold off
    filename= [filedirectory 'Elongation&Isotropy\PAFRAC' int2str(site) '.png'];
    saveas(p,filename,'png')
    
    
    if length(sloughperims)>1
    
    logPAperimsslough=log10(sloughperims(sloughareas>=PAmin1));
    logPAareasslough=log10(sloughareas(sloughareas>=PAmin1));
    
    
    fit_func2=@(a,x) a(1)+a(2)*x.^a(3);
    [PA_fit2,PA_r2,PA_j,PA_covb,PA_mse]=nlinfit(logPAperimsslough,logPAareasslough,fit_func2,[1 1 1]);
    fit_func2=@(x) PA_fit2(1)+PA_fit2(2).*x.^PA_fit2(3);
    PAr2=1-(sum(PA_r2.^2)/sum((logPAareasslough-mean(logPAareasslough)).^2));
    
   
    
    fit_funcslough=@(a,x) a(1)+a(2)*x;
    [PA_fitslough,PA_rslough,PA_j,PA_covb,PA_mse]=nlinfit(logPAperimsslough,logPAareasslough,fit_funcslough,[0 1]);
    fit_funcslough1=@(x) PA_fitslough(1)+PA_fitslough(2).*x;
    PAFRACslough=2/PA_fitslough(2);
    PArslough=1-(sum(PA_rslough.^2)/sum((logPAareasslough-mean(logPAareasslough)).^2));
    
    logPAperimsslough=log10(sloughperims(sloughareas>=1000));
    logPAareasslough=log10(sloughareas(sloughareas>=1000));
    
    [PA_fitslough,PA_rslough,PA_j,PA_covb,PA_mse]=nlinfit(logPAperimsslough,logPAareasslough,fit_funcslough,[0 1]);
    
    PAFRACslough5=2/PA_fitslough(2);
    PArslough=1-(sum(PA_rslough.^2)/sum((logPAareasslough-mean(logPAareasslough)).^2));
    
    logPAperimsslough=log10(sloughperims(sloughareas>=10000));
    logPAareasslough=log10(sloughareas(sloughareas>=10000));
    
    [PA_fitslough,PA_rslough,PA_j,PA_covb,PA_mse]=nlinfit(logPAperimsslough,logPAareasslough,fit_funcslough,[0 1]);

    PAFRACslough10=2/PA_fitslough(2);
    PArslough=1-(sum(PA_rslough.^2)/sum((logPAareasslough-mean(logPAareasslough)).^2));
    
    
    
    p=figure(1);
    plot(logPAperimsslough,logPAareasslough,'.');
    xlim([min(logPAperimsslough) ceil(max(logPAperimsslough))]);
    ylim([min(logPAareasslough) ceil(max(logPAareasslough))]);
    title(['Slough Perimeter VS Areas (Site ' num2str(site) ')'])
    xlabel('log(Perimeter)')
    ylabel('log(Patch Area)')
    hold on
    fplot(fit_func2,[min(logPAperimsslough) max(logPAperimsslough)])
    exponenttext2 = sprintf('{%d}',PA_fit2(3));
    PAexpslough=PA_fit2(3);
    text(3,6.5,[num2str(PA_fit2(1)),'+',num2str(PA_fit2(2)),'*X^',exponenttext2])
    text(3,6.25,['r^2 = ',num2str(PAr2)])
    fplot(fit_funcslough1,[min(logPAperimsslough) max(logPAperimsslough)],'r')
    text(3.5,5,[num2str(PA_fitslough(1)),'+',num2str(PA_fitslough(2)),'*X'])
    text(3.5,4.75,['r^2 = ',num2str(PArslough)])
    text(3.5,4.5,['PAFRAC(5) = ',num2str(PAFRACslough)])
    hold off
    filename= [filedirectory 'Elongation&Isotropy\PAFRAC_slough' int2str(site) '.png'];
    saveas(p,filename,'png')
    
    else
        PAexpslough=nan; 
        PAFRACslough=nan;
    end
      
%% Patch-Based Regression

    patchcutoff=PAmin1;
    
    lw1=lw_grouped(grouped(:,1)>=patchcutoff);
    lw2=lw_grouped2(grouped(:,1)>=patchcutoff);
    lwmax=lw_max(grouped(:,1)>=patchcutoff);
    patchareas=PAareas(grouped(:,1)>=patchcutoff);
    logperims=log10(PAperims(PAareas>=patchcutoff));
    logareas=log10(PAareas(PAareas>=patchcutoff));

    input=patchareas;
    
    output=lw1;
    
    
    fit_func1=@(a,x) a(1)+a(2)*log10(x);
    [fit,resid,jacob,covb,mserror]=nlinfit(input,output,fit_func1,[0 1]);
    rsquared=1-(sum(resid.^2)/sum((output-mean(output)).^2));
    lw1_r=rsquared;
    lw1_c=fit;
    lw1_p=linhyptest(fit(2),covb(2,2));
    
    
    
    output=lw2;
    [fit,resid,jacob,covb,mserror]=nlinfit(input,output,fit_func1,[0 1]);
    rsquared=1-(sum(resid.^2)/sum((output-mean(output)).^2));
    lw2_r=rsquared;
    lw2_c=fit;
    lw2_p=linhyptest(fit(2),covb(2,2));
    
    
    
    output=lwmax;
    stats=regstats(output,input,'linear',{'rsquare','adjrsquare','fstat','beta'});
    lwmax_r=stats.rsquare;
    lwmax_c=stats.beta;
    lwmax_p=stats.fstat.pval;
    
    
    
    
    fit_func1=@(a,x) a(1)+a(2)*x;
    fit_func2=@(a,x) a(1)+a(2)*x+a(3)*x.^2;
    fit_func3=@(a,x) a(1)+a(2)*x.^a(3);
    
    input=logareas;
    output=logperims;
    
    stats=regstats(output,input,'linear',{'rsquare','adjrsquare','fstat','beta'});
    PA_p=stats.fstat.pval;
    
    [fit,resid,jacob,covb,mserror]=nlinfit(input,output,fit_func1,[0 1]);
    rsquared=1-(sum(resid.^2)/sum((output-mean(output)).^2));
    PA_r=rsquared;
    PA_ar1=1-(1-rsquared)*((length(input)-1)/(length(input)-length(fit)));
    PA_c1=fit;
    [fit,resid,jacob,covb,mserror]=nlinfit(input,output,fit_func2,[0 1 1]);
    rsquared=1-(sum(resid.^2)/sum((output-mean(output)).^2));
    PA_r2=rsquared;
    PA_ar2=1-(1-rsquared)*((length(input)-1)/(length(input)-length(fit)));
    PA_c2=fit;
    [fit,resid,jacob,covb,mserror]=nlinfit(input,output,fit_func3,[0 1 1]);
    rsquared=1-(sum(resid.^2)/sum((output-mean(output)).^2));
    PA_r3=rsquared;
    PA_ar3=1-(1-rsquared)*((length(input)-1)/(length(input)-length(fit)));
    PA_c3=fit;
   
    
    
    
    dfj=length(input)-3;
    linquad_f=dfj*(PA_r2-PA_r)/(1-PA_r2);
    linquad_p=1-fcdf(linquad_f,2,dfj);
    linquad_p2=1-fcdf(linquad_f,2,3);
    
    linpow_f=dfj*(PA_r3-PA_r)/(1-PA_r3);
    linpow_p=1-fcdf(linpow_f,2,dfj);
    linpow_p2=1-fcdf(linpow_f,2,3);
    
        logx=logareas;
        logx_x=min(logx):(max(logx)-min(logx))/1000:max(logx);
        x2=10.^logx_x';
        liny=10.^fit_func1(PA_c1,logx_x);
        quady=10.^fit_func2(PA_c2,logx_x);
    
   
        
        
        
    output=log10(sloughperims(sloughareas>=PAmin1));
    input=log10(sloughareas(sloughareas>=PAmin1));
    
    [PAslough_fit1,PAslough_r1,PAslough_j1,PAslough_covb1,PAslough_mse1]=nlinfit(input,output,fit_func1,[0 1]);
    rsquared=1-(sum(PAslough_r1.^2)/sum((output-mean(output)).^2));
    PAslough_r1=rsquared;
    PAslough_ar1=1-(1-rsquared)*((length(input)-1)/(length(input)-length(PAslough_fit1)));
    PAslough_c1=PAslough_fit1;
    
    [PAslough_fit2,PAslough_r2,PAslough_j2,PAslough_covb2,PAslough_mse2]=nlinfit(input,output,fit_func2,[1 1 1]);
    rsquared=1-(sum(PAslough_r2.^2)/sum((output-mean(output)).^2));
    PAslough_r2=rsquared;
    PAslough_ar2=1-(1-rsquared)*((length(input)-1)/(length(input)-length(PAslough_fit2)));
    PAslough_c2=PAslough_fit2;
   
    dfj=length(input)-3;
    linquadslough_f=dfj*(PAslough_r2-PAslough_r1)/(1-PAslough_r2);
    linquadslough_p=1-fcdf(linquadslough_f,2,dfj);
    linquadslough_p2=1-fcdf(linquadslough_f,2,3);
        
        logx=input;
        logx_x=min(logx):(max(logx)-min(logx))/1000:max(logx);
        x2slough=10.^logx_x';
        linyslough=10.^fit_func1(PAslough_c1,logx_x);
        quadyslough=10.^fit_func2(PAslough_c2,logx_x);
    
        
    
    
    excelfilename=[filedirectory 'PatchRegression.xlsx'];
    if site==1
        delete(excelfilename)
    end
    
    header={'Site' 'LW1' 'LW2' 'LWmax' 'PA' 'PA1500'};
    xlswrite(excelfilename,header,'p-value','A1');
    xlswrite(excelfilename,header,'rsquared','A1');

    excel_index = ['A' int2str(1+site)];
    xlswrite(excelfilename,[site lw1_p lw2_p lwmax_p PA_p],'p-value',excel_index);  
    xlswrite(excelfilename,[site lw1_r lw2_r lwmax_r PA_r],'rsquared',excel_index);  
    xlswrite(excelfilename,[site linquad_p linquad_p2 linquadslough_p linquadslough_p2],'L vs Q p-value',excel_index);  
    xlswrite(excelfilename,[site linpow_p linpow_p2],'L vs P p-value',excel_index);  
    
    header={'Site' 'LW1_1' 'LW1_2' 'LW2_1' 'LW2_2' 'LWmax_1' 'LWmax_2' 'PA1' 'PA1' 'PA2' 'PA2' 'PA2' 'PA3'  'PA3'  'PA3' 'PAslough1' 'PAslough1' 'PAslough2' 'PAslough2' 'PAslough2'};
    xlswrite(excelfilename,header,'coefficients','A1');
    xlswrite(excelfilename,[site lw1_c lw2_c lwmax_c' PA_c1 PA_c2 PA_c3 PAslough_c1 PAslough_c2],'coefficients',excel_index);  
    
    header={'Site' 'PA1' 'PA2' 'PA3' 'Slough1' 'Slough2'};
    xlswrite(excelfilename,header,'adj rsquared','A1');
    xlswrite(excelfilename,[site PA_ar1 PA_ar2 PA_ar3 PAslough_ar1 PAslough_ar2],'adj rsquared',excel_index);  
    
    header={['Site ' int2str(site)]};
    excel_index = [idx2A1(site) '1'];
    xlswrite(excelfilename,header,'areas',excel_index);
    xlswrite(excelfilename,header,'perims',excel_index);
    xlswrite(excelfilename,header,'LW1',excel_index);
    xlswrite(excelfilename,header,'LW2',excel_index);
    xlswrite(excelfilename,header,'LWmax',excel_index);
    xlswrite(excelfilename,header,'orientation',excel_index);
    
    excel_index = [idx2A1(site) '2'];  
    xlswrite(excelfilename,PAareas,'areas',excel_index);
    xlswrite(excelfilename,PAperims,'perims',excel_index);
    xlswrite(excelfilename,lw_grouped,'LW1',excel_index);  
    xlswrite(excelfilename,lw_grouped2,'LW2',excel_index);  
    xlswrite(excelfilename,lw_max,'LWmax',excel_index);  
    xlswrite(excelfilename,ridgestats(:,5),'orientation',excel_index); 
    xlswrite(excelfilename,x2,'x 2',excel_index); 
    xlswrite(excelfilename,liny','lin y',excel_index);  
    xlswrite(excelfilename,quady','quad y',excel_index); 
    xlswrite(excelfilename,x2slough,'x 2slough',excel_index); 
    xlswrite(excelfilename,linyslough','lin yslough',excel_index);  
    xlswrite(excelfilename,quadyslough','quad yslough',excel_index); 
    
    
%% distribution plots


    
    areas1=grouped(:,1);
    areas1=areas1(areas1>=paretocutoff);
    pareas=ccdfcalc(areas1);

if length(areas1)>10
    
    fit_func=@(a,x) (x/min(areas1)).^(-a+1);
    [regress_fit,areas_r,areas_j,~,areas_mse]=nlinfit(unique(areas1),pareas,fit_func,[2]);
   
    fit_func2=@(x) (x/min(areas1)).^(-regress_fit+1);
    areasr=1-(sum(areas_r.^2)/sum((pareas-mean(pareas)).^2));
    
    p=figure(1);
    loglog(unique(areas1),pareas,'.')
    xlim([min(areas1) ceil(max(areas1))]);
    ylim([min(pareas) ceil(max(pareas))]);
    title(['Ridge Area CDF (regression) (Site ' num2str(site) ')'])
    xlabel('Patch size')
    ylabel('P(X > x)')
    hold on
    fplot(fit_func2,[min(areas1) max(areas1)])
    text(100,.1,['r^2 = ',num2str(areasr)])
    text(100,.2,['a = ',num2str(regress_fit)])
    filename= [filedirectory 'CDFs\Areas_CDF2_' int2str(site) '.png'];
    saveas(p,filename,'png')   
    hold off
   
    
    nx=length(areas1);
    xbound=paretocutoff-1;
    fit1=1+nx*(1/(sum(log(areas1./xbound))));
    pareto_pdf=@(x,alpha) ((alpha-1)*xbound^(alpha-1))*x.^-alpha;
   
    flike=@(fit) -sum(log(pareto_pdf(areas1,fit)));
    pareto_fit=fminsearch(flike, fit1)
    
    
    d_pareto_fit=1+nx*(1/(sum(log(areas1./(xbound-.5)))));
%     discrete_like=@(alpha) -(-nx*log(zeta(alpha))-alpha*sum(log(areas1)));
%     pareto_fit=fminsearch(discrete_like,fit1)
    
    
    fit_func2=@(x) (x/xbound).^(-pareto_fit+1);
    areasr=1-(sum((pareas-fit_func2(unique(areas1))).^2)/sum((pareas-mean(pareas)).^2));

   
    
    p=figure(1);
    loglog(unique(areas1),pareas,'.');
    title(['Ridge Areas CDF (mle) (Site ' num2str(site) ')'])
    xlabel('Patch Size')
    ylabel('P(X > x)')
    hold on
    fplot(fit_func2,[min(areas1) max(areas1)])
    text(100,.1,['r^2 = ',num2str(areasr)])
    text(100,.2,['a = ',num2str(pareto_fit)])
    hold off
    filename= [filedirectory 'CDFs\Areas_CDF1_' int2str(site) '.png'];
    saveas(p,filename,'png')   
    
    areas1=sloughareas;
    areas1=areas1(areas1>=paretocutoff);
    pareas=ccdfcalc(areas1);
    
    
    fit_func=@(a,x) (x/min(areas1)).^(-a+1);
    [sloughregress_fit,areas_r,areas_j,areas_covb,areas_mse]=nlinfit(unique(areas1),pareas,fit_func,[2]);
   
    fit_func2=@(x) (x/min(areas1)).^(-sloughregress_fit+1);
    areasr=1-(sum(areas_r.^2)/sum((pareas-mean(pareas)).^2));
    
    p=figure(1);
    loglog(unique(areas1),pareas,'.')
    xlim([min(areas1) ceil(max(areas1))]);
    ylim([min(pareas) ceil(max(pareas))]);
    title(['Slough Area CDF (regression) (Site ' num2str(site) ')'])
    xlabel('Patch size')
    ylabel('P(X > x)')
    hold on
    fplot(fit_func2,[min(areas1) max(areas1)])
    text(10000,.04,['r^2 = ',num2str(areasr)])
    text(10000,.03,['a = ',num2str(sloughregress_fit)])
    filename= [filedirectory 'CDFs\SloughAreas_CDF2_' int2str(site) '.png'];
    saveas(p,filename,'png')   
    hold off
    

    nx=length(areas1);
    xbound=paretocutoff;
    fit1=1+nx*(1/(sum(log(areas1./xbound))));
    sloughpareto_pdf=@(x,alpha) ((alpha-1)*xbound^(alpha-1))*x.^-alpha;
    flike=@(fit) -sum(log(sloughpareto_pdf(areas1,fit)));
    sloughpareto_fit=fminsearch(flike, fit1);
    
    fit_func2=@(x) (x/min(areas1)).^(-sloughpareto_fit+1);
    areasr=1-(sum((pareas-fit_func2(unique(areas1))).^2)/sum((pareas-mean(pareas)).^2));

    p=figure(1);
    loglog(unique(areas1),pareas,'.');
    title(['Slough Areas CDF (mle) (Site ' num2str(site) ')'])
    xlabel('Patch Size')
    ylabel('P(X > x)')
    hold on
    fplot(fit_func2,[min(areas1) max(areas1)])
    text(10000,.04,['r^2 = ',num2str(areasr)])
    text(10000,.03,['a = ',num2str(sloughpareto_fit)])
    hold off
    filename= [filedirectory 'CDFs\SloughAreas_CDF1_' int2str(site) '.png'];
    saveas(p,filename,'png')   
    
    p=figure(1);
    semilogy(unique(wridge'),ccdfcalc(wridge'),'.');
    title(['Ridge Widths CDF (Site ' num2str(site) ')'])
    xlabel('Width')
    ylabel('P(X > x)')  
    filename= [filedirectory 'CDFs\wridge_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')     
    p=figure(1);
    semilogy(unique(lridge'),ccdfcalc(lridge'),'.');
    title(['Ridge Lengths CDF (Site ' num2str(site) ')'])
    xlabel('Length')
    ylabel('P(X > x)')    
    filename= [filedirectory 'CDFs\lridge_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')   
    p=figure(1);
    semilogy(unique(wridge_ok'),ccdfcalc(wridge_ok'),'.');
    title(['Ridge WidthsOK CDF (Site ' num2str(site) ')'])
    xlabel('Width')
    ylabel('P(X > x)')    
    filename= [filedirectory 'CDFs\wridge_ok_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')   
    p=figure(1);
    semilogy(unique(lridge_ok'),ccdfcalc(lridge_ok'),'.');
    title(['Ridge LengthsOK CDF (Site ' num2str(site) ')'])
    xlabel('Length')
    ylabel('P(X > x)')
    filename= [filedirectory 'CDFs\lridge_ok_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')   
    p=figure(1); 
    semilogy(unique(wslough'),ccdfcalc(wslough'),'.');
    title(['Slough Widths CDF (Site ' num2str(site) ')'])
    xlabel('Width')
    ylabel('P(X > x)')    
    filename= [filedirectory 'CDFs\wslough_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')   
    p=figure(1);
    semilogy(unique(lslough'),ccdfcalc(lslough'),'.');
    title(['Slough Lengths CDF (Site ' num2str(site) ')'])
    xlabel('Length')
    ylabel('P(X > x)')
    filename= [filedirectory 'CDFs\lslough_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')  
    p=figure(1); 
    semilogy(unique(wnn),ccdfcalc(wnn),'.');
    title(['Ridge ENN Widths CDF (Site ' num2str(site) ')'])
    xlabel('Width ENN')
    ylabel('P(X > x)')    
    filename= [filedirectory 'CDFs\wnn_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')   
    p=figure(1);
    semilogy(unique(lnn),ccdfcalc(lnn),'.');
    title(['Ridge ENN Lengths CDF (Site ' num2str(site) ')'])
    xlabel('Length')
    ylabel('P(X > x)')
    filename= [filedirectory 'CDFs\lnn_CDF' int2str(site) '.png'];
    saveas(p,filename,'png')   
    
    
    

end
    
    
    else
        PAFRAC1=nan;
        PAFRAC5=nan;
        PAFRAC10=nan;
        lwscaling=nan;
        lwscalingmax=nan;
        PAexp=nan;
        PAFRACslough=nan;
        PAFRACslough5=nan;
        PAFRACslough10=nan;
        PAexpslough=nan;
        pareto_fit=nan;
        regress_fit=nan;
        d_pareto_fit=nan;
        
    end
    
    %
    %
    % figure(5)
    % subplot(2,1,1), semilogy(wslough,wslough_p,'.');
    % title('Slough Widths (semi-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(wslough,wslough_p,'.');
    % title('Slough Widths (log-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    %
    % figure(6)
    % subplot(2,1,1), semilogy(lslough,lslough_p,'.');
    % title('Slough Lengths (semi-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(lslough,lslough_p,'.');
    % title('Slough Lengths (log-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    %
    % figure(7)
    % subplot(2,1,1), semilogy(lridge,lridge_p,'.');
    % title('Ridge Length (semi-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(lridge,lridge_p,'.');
    % title('Ridge Lengths (log-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    %
    % figure(8)
    % subplot(2,1,1), semilogy(wridge,wridge_p,'.');
    % title('Ridge Widths (semi-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(wridge,wridge_p,'.');
    % title('Ridge Widths (log-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    %
    % figure(9)
    % subplot(2,1,1), semilogy(grouped(:,1),grouped_p(:,1),'.');
    % title('Areas (semi-log)')
    % xlabel('Area')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(grouped(:,1),grouped_p(:,1),'.');
    % title('Areas (log-log)')
    % xlabel('Area')
    % ylabel('P(X > x)')
    %
    % figure(10)
    % cshow(c);
    
    
    % figure(1)
    % subplot(2,1,1), semilogy(lridge_tree,lridge_tree_p,'.');
    % title('Tree Islands Length (semi-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(lridge_tree,lridge_tree_p,'.');
    % title('Tree Islands Lengths (log-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    %
    % figure(2)
    % subplot(2,1,1), semilogy(wridge_tree,wridge_tree_p,'.');
    % title('Tree Islands Widths (semi-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(wridge_tree,wridge_tree_p,'.');
    % title('Tree Islands Widths (log-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    %
    % figure(3)
    % subplot(2,1,1), semilogy(lridge_ok,lridge_ok_p,'.');
    % title('Ridge Length without Tree Islands (semi-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(lridge_ok,lridge_ok_p,'.');
    % title('Ridge Lengths without Tree Islands (log-log)')
    % xlabel('Length')
    % ylabel('P(X > x)')
    %
    % figure(4)
    % subplot(2,1,1), semilogy(wridge_ok,wridge_ok_p,'.');
    % title('Ridge Widths without Tree Islands (semi-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    % subplot(2,1,2), loglog(wridge_ok,wridge_ok_p,'.');
    % title('Ridge Widths without Tree Islands (log-log)')
    % xlabel('Width')
    % ylabel('P(X > x)')
    %
    % figure(5)
    
    
    
    
    %% Construct Distribution matrix
    
    l1=length(grouped(:,1));
    l2=length(a_wridge(:,1));
    l3=length(a_lridge(:,1));
    l4=length(a_wslough(:,1));
    l5=length(a_lslough(:,1));
    l6=length(wridge_tree);
    l7=length(lridge_tree);
    l8=length(wridge_ok);
    l9=length(lridge_ok);
    l10=length(wnn);
    l11=length(lnn);
    
    l_all = [l1 l2 l3 l4 l5 l6 l7 l8 l9 l10 l11];
    dist_array=zeros(max(l_all),11);
    
    dist_array(1:l1,1)=grouped(:,1);
    dist_array(1:l2,2)=a_wridge(:,1);
    dist_array(1:l3,3)=a_lridge(:,1);
    dist_array(1:l4,4)=a_wslough(:,1);
    dist_array(1:l5,5)=a_lslough(:,1);
    dist_array(1:l6,6)=wridge_tree;
    dist_array(1:l7,7)=lridge_tree;
    dist_array(1:l8,8)=wridge_ok;
    dist_array(1:l9,9)=lridge_ok;
    dist_array(1:l10,10)=wnn;
    dist_array(1:l11,11)=lnn;
    

    
    %% Write Area Distributions
    
    % areas=(grouped(:,1));
    % cdf_areas=ecdfcalc(areas);
    % cdf_areas_all=zeros(length(areas),1);
    % for n=1:length(areas)
    %     cdf_areas_all(n)=cdf_areas(areas(n)==unique(areas));
    % end;
    %
    % area_dist=[areas 1-cdf_areas_all grouped_ti];
    filename= [filedirectory 'OutputData\dist_array' int2str(site) '.mat'];
    delete(filename)
    sloughpareto_fit=0;
    sloughregress_fit=0;
    save(filename,'dist_array','ridgeperims', 'cattail', 'density','orientation','PAFRAC1','PAFRAC5','PAFRAC10', 'lwscaling', 'PAexp', 'PAexpslough', 'PAFRACslough', 'PAFRACslough5', 'PAFRACslough10', 'sloughperims','sloughareas', 'lwscalingmax')
    
    


    
end








%% RSL Metrics
%Calculates & writes the final metrics in an excel file

clearvars -except paretocutoff mapsize mapn xmin filedirectory mapname


excelfilename=[filedirectory 'RegressionInputs.xlsx'];


delete(excelfilename)
header={'Site' 'Easting' 'Northing' 'PTI' 'AvgDepth' 'Slough %' 'Ridge %' 'TI %' 'TIaff %' 'RidgeAll %' 'Slough L' 'Slough W' 'Ridge L' 'Ridge W' 'RidgeOK L' 'RidgeOK W' 'Ridge MLE' 'RidgeRegress' 'Slough MLE' 'SloughRegress' 'Avg Ridge Size' 'Avg Slough Size' 'Powexp Cutoff' 'Hypexp Cutoff' 'L/W Scaling' 'Avg L/W' 'Avg RidgeOK L/W' 'PAFRAC100' 'PAFRAC1000' 'PAFRAC10000' 'PAFRAC100sloughs' 'PAFRAC1000sloughs' 'PAFRAC10000sloughs' 'TI/TIall' 'Cattail' 'Discrete Plaw Estimate' 'Total Edge'};
xlswrite(excelfilename,header,'RegressionInputs','A1');

xbound=paretocutoff;

for site=1:mapn
    disp(site);
    filename= [filedirectory 'OutputData\dist_array' int2str(site) '.mat'];
    load(filename)
    
    
    excel_index = ['F' int2str(1+site)];
    xlswrite(excelfilename,density(1),'RegressionInputs',excel_index);  
    pause(1)
    excel_index = ['G' int2str(1+site)];
    xlswrite(excelfilename,density(2),'RegressionInputs',excel_index); 
    pause(1)
    excel_index = ['H' int2str(1+site)];
    xlswrite(excelfilename,density(4),'RegressionInputs',excel_index);    
    pause(1)
    excel_index = ['I' int2str(1+site)];
    xlswrite(excelfilename,density(3)+density(4),'RegressionInputs',excel_index); 
    pause(1)
    excel_index = ['J' int2str(1+site)];
    xlswrite(excelfilename,density(2)+density(3)+density(4),'RegressionInputs',excel_index);
    pause(1)
    excel_index = ['K' int2str(1+site)];
    x=dist_array(:,5);
    xlswrite(excelfilename,mean(x(x>0)),'RegressionInputs',excel_index);
    pause(1)
    excel_index = ['L' int2str(1+site)];
    x=dist_array(:,4);
    xlswrite(excelfilename,mean(x(x>0)),'RegressionInputs',excel_index);  
    pause(1)
    excel_index = ['M' int2str(1+site)];
    x=dist_array(:,3);
    ridgel=mean(x(x>0));
    xlswrite(excelfilename,mean(x(x>0)),'RegressionInputs',excel_index);  
    pause(1)
    excel_index = ['N' int2str(1+site)];
    x=dist_array(:,2);
    ridgew=mean(x(x>0));
    xlswrite(excelfilename,mean(x(x>0)),'RegressionInputs',excel_index);
    pause(1)
    excel_index = ['O' int2str(1+site)];
    x=dist_array(:,9);
    ridgeOKl=mean(x(x>0));
    xlswrite(excelfilename,mean(x(x>0)),'RegressionInputs',excel_index);  
    excel_index = ['P' int2str(1+site)];
    x=dist_array(:,8);
    ridgeOKw=mean(x(x>0));
    xlswrite(excelfilename,mean(x(x>0)),'RegressionInputs',excel_index); 
    
    

    
    
    excel_index = ['U' int2str(1+site)];
    xlswrite(excelfilename,mean(dist_array(:,1)),'RegressionInputs',excel_index);
    excel_index = ['V' int2str(1+site)];
    xlswrite(excelfilename,mean(sloughareas),'RegressionInputs',excel_index);
    


    excel_index = ['Y' int2str(1+site)];
    xlswrite(excelfilename,lwscaling,'RegressionInputs',excel_index);    
    excel_index = ['Z' int2str(1+site)];
    xlswrite(excelfilename,(ridgel/ridgew),'RegressionInputs',excel_index);    
    excel_index = ['AA' int2str(1+site)];
    xlswrite(excelfilename,(ridgeOKl/ridgeOKw),'RegressionInputs',excel_index);
    
    
    excel_index = ['AB' int2str(1+site)];
    xlswrite(excelfilename,PAFRAC1,'RegressionInputs',excel_index);
    excel_index = ['AC' int2str(1+site)];
    xlswrite(excelfilename,PAFRAC5,'RegressionInputs',excel_index);
    excel_index = ['AD' int2str(1+site)];
    xlswrite(excelfilename,PAFRAC10,'RegressionInputs',excel_index);
    excel_index = ['AE' int2str(1+site)];
    xlswrite(excelfilename,PAFRACslough,'RegressionInputs',excel_index);
    
    
    excel_index = ['AF' int2str(1+site)];
    xlswrite(excelfilename,PAFRACslough5,'RegressionInputs',excel_index);
    excel_index = ['AG' int2str(1+site)];
    xlswrite(excelfilename,PAFRACslough10,'RegressionInputs',excel_index);
    
     
    excel_index = ['AH' int2str(1+site)];
    xlswrite(excelfilename,density(4)/(density(4)+density(3)),'RegressionInputs',excel_index);  
    
    excel_index = ['AI' int2str(1+site)];
    xlswrite(excelfilename,cattail,'RegressionInputs',excel_index);    
    
    
    excel_index = ['AK' int2str(1+site)];
    xlswrite(excelfilename,sum(ridgeperims),'RegressionInputs',excel_index); 
end

excelFileName = 'RegressionInputs';
excelFilePath = filedirectory ; 
sheetName = 'Sheet'; 


objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(fullfile(excelFilePath, excelFileName)); 


try
objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
catch

end

objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;


