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
disp('Select max-projected stack of images  ')
gray_path=uigetdir(pwd,'select max-projected stack of images');
gray_ims=dir([gray_path,filesep,'*tif']);
gray_ims=gray_ims(~ismember({gray_ims.name},{'.','..'})); %delete hidden files

convert=0;
if isempty(gray_ims)
    convert=1;
    gray_ims=dir([gray_path,filesep,'*png']);
    gray_ims=gray_ims(~ismember({gray_ims.name},{'.','..'})); %delete hidden files
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load saved data
disp('Please guide to the final saved .mat file')
[file,path]=uigetfile('*.mat');
load(fullfile(path,file))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Select where you want to save the printed images to.
disp('Select directory you wish to save printed images to ')
save_path=uigetdir(pwd,'Select directory you wish to save printed images to');
if strcmp(save_path,gray_path)
    error('You do not want gray path and save path to be the same')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %recursively ask user to identify z stack in consecutive order, until they
% %hit cancel. 
% last_check=0;
% while last_check==0
%     p=0;
%     Fpaths={};
%     while p==0
%         disp("Select z-stacks in order. If complete or made an error, press cancel")
%         test=uigetdir(pwd,'select z-stacks in order');
%         if isstr(test)~=0
%             Fpaths{end+1}=test;
%         else
%             p=1;
%         end
%     end
%     Fpaths=Fpaths';
%     disp(Fpaths)
%     input_hf = str2double(input('Is the printed order of directories correct? (Enter an integer) 1: yes, 2: no    ', 's'));
%     while isnan(input_hf) || fix(input_hf) ~= input_hf || input_hf<1 || input_hf>2
%       input_hf = str2double(input('Is the printed order of directories correct? (Enter an integer) 1: yes, 2: no    ', 's'));
%     end
%     if input_hf==1
%         last_check=1;
%     else
%         disp("Having user reselect the z-stacks...")
%     end
% end

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

%Code to process images and correctly identify cell numbers. 

%for each path, we want to analyze the ith frame, measuring each object and
%assigning a z along with it. These objects have cell numbers (indices)
%assigned to them depending only on where they appear on the image (l->r) by default.
%We want to correct these numbers so that we can accurately track the
%trajectory of any individual cell. 

% Fnames=cell(1,length(Fpaths));
% for fpath_num=1:length(Fpaths)
%     d_ims=dir([Fpaths{fpath_num},filesep,'*tif']);
%     d_ims=d_ims(~ismember({d_ims.name},{'.','..'})); %delete hidden files
%     if isempty(d_ims)
%         d_ims=dir([Fpaths{fpath_num},filesep,'*png']);
%         d_ims=d_ims(~ismember({d_ims.name},{'.','..'})); %delete hidden files
%     end
%     Fnames{fpath_num}=d_ims;
% end

cellmax=0;
%maxframe=length(Fnames{1}); %should double check all these filepaths have the same number of files in them
maxframe=length(gray_ims);

for im_num=1:maxframe
    %     for fpath_num=1:length(Fpaths)
    %         %load BW, make logical, run stats, positions, circles, frame number,
    %         %get cell numbers (add cellmax), assign z
    %         try
    %             fullpathname=[Fnames{fpath_num}(im_num).folder, filesep, Fnames{fpath_num}(im_num).name];
    %             fprintf('Reading %s\n',fullpathname)
    %             BW=imread(fullpathname);
    %         catch
    %             fullpathname=[Fpaths{fpath_num}, filesep, Fnames{fpath_num}(im_num).name];
    %             fprintf('Reading %s\n',fullpathname)
    %             BW=imread([Fpaths{fpath_num}, filesep, Fnames{fpath_num}(im_num).name]); %2016 Matlab doesn't have folder in Dir command.
    %         end
    %         BW=im2bw(BW,0.5);
    %         bpixels=[BW(:,1)',BW(:,end)',BW(1,:),BW(end,:)];
    %         if mean(bpixels)>0.5
    %             %means pixels were flipped.
    %             try
    %                 BW=imread(fullpathname);
    %             catch
    %                 BW=imread([Fpaths{fpath_num}, filesep, Fnames{fpath_num}(im_num).name]); %2016 Matlab doesn't have folder in Dir command.
    %             end
    %             BW=im2bw(BW,0.5);
    %             BW=imcomplement(BW);
    %         end
    %         BW=imfill(BW,'holes');
    %         BW=bwareaopen(BW,300,4); %sometimes in images two objects will be connected by 8 connectivity.
    %         BW=imclearborder(BW);
    %         
    %         if fpath_num==1
    %             CC = bwconncomp(BW>0,4);
    %             stats = regionprops(CC, {'Centroid','PixelIdxList'});
    %             pos = cat(1,stats.Centroid);
    %             
    %             centers = Max_Fit_Circles(BW); %see function
    %             
    %             OneFrameStats=stats;
    %             OneFramePos=pos;
    %             OneFrameCenters=centers;
    %         else
    %             CC = bwconncomp(BW>0,4);
    %             stats = regionprops(CC, {'Centroid','PixelIdxList'});
    %             
    %             if ~isempty(stats)
    %                 pos = cat(1,stats.Centroid);
    %                 
    %                 OneFrameStats=[OneFrameStats;stats];
    %                 OneFramePos=[OneFramePos;pos];
    %                 OneFrameCenters = [OneFrameCenters; Max_Fit_Circles(BW)];
    %             end
    %         end
    %     end
    
    OneFrameStats = AllFrameStats(AllFrames==im_num);
    OneFramePos = cat(1,OneFrameStats.Centroid);
    OneFrameCenters = AllFrameCenters(AllFrames==im_num,:);
    
    Igray=imread([gray_path,filesep,gray_ims(im_num).name]);
    
    if convert==1
        Igray=im2uint8(Igray);
    end
    
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
    imwrite(im,[save_path,filesep,gray_ims(im_num).name])
    
    
        
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%rgbImage = cat(3, grayImage, grayImage, grayImage);

%
%image(TheGrayscaleImage);
%colormap(gray(256));
%hold on
%image(TheColorImage, 'AlphaData', YourTransparencyData)
%YourTransparencyData should be 1 in the places you want TheColorImage to completely show up, 
%0 in the places where you want the grayscale image to be completely visible 
%and the color image not visible, and anything inbetween for the places you 
%want to "blend" the two images.

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
