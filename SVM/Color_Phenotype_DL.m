%write a script to visualize the trajectories. To do this, we want to use
%much of the same information that is already inside of Overlap

%okay. So...
%you'll have to point to all the usual directories so you can load the data
%from each frame. You'll need Centroid or Centers, as well as PixelIdxList.

%from there LOAD the data in the text file.
%determine the classification for each cell for each frame. 

%next, for each frame
% for every object in that frame (found in text file) 
% find deltar to all other centroids. 
%whichever is closest, use that PixelIdxList to color the image.


colors = {[1,1,0],[1,0,1],[0,1,0],[0,0,1],[0, 1, 1]};
%{[1,0,1],[1,1,0],[0,0,1],[1,0,0],[0, 1, 1]};
%{[1,1,0],[1,0,1],[0,1,0],[0,0,1],[0, 1, 1]};
%{[1,1,0],[1,0,1],[0,1,0],[0,0,1],[1,0,0], [0, 1, 1]};

%magenta yellow blue red  cyan
% AE       FP   BB   LA   IM

%yellow magenta green blue red cyan
% AE     FP      HB    LA   SB  IM

%yellow magenta green blue  cyan
% AE     FP      HB    LA    IM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Check that external functions exist within the directory.
DefaultDir=pwd;
if ~exist([DefaultDir, '/max_inscribed_circle.m'])
    error('The Matlab script "max_inscribed_circle.m" is not in the same directory')
end
if ~exist([DefaultDir, '/inpoly.m'])
    error('The Matlab script "inpoly.m" is not in the same directory')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We want to ask the user to identify the max projected stack of images
%which we will use to identify the cells. 
disp('Select where the .mat files are that were saved from cell_v2.py  ')
mat_path=uigetdir(pwd,'Select where the .mat files are that were saved from cell_v2.py');

d_ims=dir([mat_path,filesep,'*.mat']);
d_ims=d_ims(~ismember({d_ims.name},{'.','..'})); %delete hidden files
Matnames=d_ims;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load saved data
disp('Please guide to the final saved checkpoint .mat file after running CellTrack_DL.m ')
[file,path]=uigetfile('*.mat');
load(fullfile(path,file))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select where you want to save the printed images to.
disp('Select directory you wish to save printed images to ')
save_path=uigetdir(pwd,'Select directory you wish to save printed images to');
if strcmp(save_path,mat_path)
    error('You do not want gray path and save path to be the same')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


idcs=strfind(save_path,'/');

disp('Please select text file from SVM.py:   ')
[f,pth]=uigetfile('*.txt');
%fid = fopen([save_path(1:idcs(end)),save_path(idcs(end-1)+1:idcs(end)-1),'.txt'],'rt');
fid = fopen(fullfile(pth,f),'rt');
C = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ', 'HeaderLines',1);
fclose(fid);
%keep in mind, the columns here are: 
%frame cellnumber x y ae fp hs la sb

%run classification
Classifications=zeros(1,length(C{1}));
for i=1:length(C{1})
    M=[C{5}(i),C{6}(i),C{7}(i),C{8}(i),C{9}(i)];
    [val,in]=max(M);
    if val<0.7
        in=5;
    end
    Classifications(i)=in;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellmax=0;
maxframe=length(Matnames);

AllFrames = cat(1,AllFrameStats.Frame);

for im_num=1:maxframe

    load([Matnames(im_num).folder,filesep,Matnames(im_num).name])
    
    OneFrameStats = AllFrameStats(AllFrames==im_num);
    OneFramePos = cat(1,OneFrameStats.Centroid);
    OneFrameCenters = cat(1,OneFrameStats.MICCenter);%AllFrameCenters(AllFrames==im_num,:);
    
    Igray=image;
    
    ColorImagech1 = zeros(size(Igray));
    ColorImagech2 = zeros(size(Igray));
    ColorImagech3 = zeros(size(Igray));
    TransparencyData = zeros(size(Igray));
    
    %find indeces where current frame is.
    inds=find(C{1}==im_num);
    output_posx=C{3}(inds);
    output_posy=C{4}(inds);
    classes=Classifications(inds);
    
    for obj=1:length(inds)
        deltar=sqrt(((OneFrameCenters(:,1)-output_posx(obj)).^2)+...
            ((OneFrameCenters(:,2)-output_posy(obj)).^2));
        [~,i]=min(deltar);
        change_pixels=OneFrameStats(i).PixelIdxList;
        %the problem is right now this is linear indexed. 
        cs=colors{classes(obj)};
        ColorImagech1(change_pixels)=cs(1);
        ColorImagech2(change_pixels)=cs(2);
        ColorImagech3(change_pixels)=cs(3);
    end
    
    ColorImage=cat(3, ColorImagech1, ColorImagech2, ColorImagech3);
    cs=sum(ColorImage,3);
    TransparencyData(cs>0)=0.3;
    
    %%%CREATE AND WRITE FIGURE.
    figure(2000)
    E=imshow(Igray);
    hold on
    L=imshow(ColorImage);
    set(L,'AlphaData',TransparencyData)
    
    h=getframe;
    [framerows,framecols,~]=size(h.cdata);

    %%%FIXING OUTPUT SIZE ISSUE.
    if framerows==721
        im=h.cdata(2:end,:,:);
    else
        im=h.cdata;
    end
    %%%

    close all
    imwrite(im,[save_path,filesep,strcat(Matnames(im_num).name(1:end-4),'.tif')]) 
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define functions necessary for this algorithm

function Centers = Max_Fit_Circles(BW)
    B=bwboundaries(BW,4,'noholes'); %cell array
    [objs,~]=size(B);
    Centers = zeros(objs,2);
    for bb=1:objs
        grain=false(size(BW));
        for uu=1:length(B{bb})
            grain(B{bb}(uu,1),B{bb}(uu,2))=true;
        end
        grain=uint8(grain)*255;
        [~,cx,cy] = max_inscribed_circle(grain,[]);
        Centers(bb,1)=cx;
        Centers(bb,2)=cy;
    end
end
