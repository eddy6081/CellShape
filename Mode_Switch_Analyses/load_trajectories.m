function [lx, ly, framescheck, numcellcheck, classes] = load_trajectories(use_softmax, delete_IM)
    %{
    Author: Chris Eddy
    Date: 03/15/21
    Loading the trajectories from saved .txt file

    Input Arguments:
        use_softmax:    boolean option to use softmax classification - i.e. no
                        intermediate states

        delete_IM:      boolean option to 'delete' uncertain states, analysis proceeds
                        as if these were missing frames.

    Outputs:
        lx:     cell array of shape 1 x Number of Cells. Contains the x
                location of cell i in lx{i} over its trajectory.

        ly:     cell array of shape 1 x Number of Cells. Contains the y
                location of cell i in lx{i} over its trajectory.

        framescheck: cell array of shape 1 x Number of Cells. frames{i} contains the
                     frame numbers which cell i was observed in (may not be
                     continuous)
        
        numcellcheck: cell array of shape 1 x Number of Cells.
                      numcellcheck{i} should all be the same.

        p_vecs: cell array of shape 1 x Number of Cells. Contains the persistent
                and non-persistent vectors (in form [p_x, p_y, np_x, np_y]) for
                each frame. Non persistent vector points left from persistent as
                positive direction.

        classes: cell array of shape 1 x Number of Cells. Contains the
                 classification for each cell i and frame f in classes{i}(f). 
    %}
    %% LOAD DATA.
    %%%%COMPUTES both location and mean square displacement. 
    %the SVM file has the classification, the frame, and position. This should
    %be all we need.
    disp('Please select text file from SVM.py:   ')
    [fil,pth]=uigetfile('*.txt');
    %fid = fopen([save_path(1:idcs(end)),save_path(idcs(end-1)+1:idcs(end)-1),'.txt'],'rt');
    fid = fopen(fullfile(pth,fil),'rt');
    C = textscan(fid, '%f%f%f%f%f%f%f%f%f', 'MultipleDelimsAsOne',true, 'Delimiter',' ', 'HeaderLines',1);
    fclose(fid);
    %keep in mind, the columns here are: 
    %frame cellnumber x y ae fp hs la sb

    pixel2dist = 0.53763672; %microns per pixel

    IM_num = 5; %specify the intermediate class number. 
    %delete_IM=false;
    softmax = use_softmax;
    if softmax==false
        req_thresh=0.6;
    end

    %delete the last column if it has nan.
    if isnan(C{end}(1))
        C(end)=[];
    end
    %run classification
    Classifications=zeros(1,length(C{1}));
    for i=1:length(C{1})
        M=[];
        for el=5:length(C)
            M=[M, C{el}(i)];
        end
        [val,in]=max(M);
        if softmax==false
            if val<req_thresh
               in=IM_num;
            end
        end
        Classifications(i) = in;
    end

    %convert pixels to microns.
    locsx = C{3}(:)*pixel2dist;
    locsy = C{4}(:)*pixel2dist;
    %To see cell numbers.
    NumCells=C{2}(:);

    %Record frame numbers, since that is how we will determine continuity.
    frames=C{1}(:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %easier to use string finding patterns than anything else
    cs=int2str(Classifications);
    cs=strrep(cs,' ',''); %delete spaces
    % pausechar = strfind(classes,'23');
    %%%%%%%%%%%%%%
    if delete_IM == true  
        indsIM=strfind(cs,num2str(IM_num));
        cs(indsIM)=[];
        frames(indsIM)=[];
        NumCells(indsIM)=[];
        locsx(indsIM)=[];
        locsy(indsIM)=[];
    end

    %delete these items.
    %%%%%%%%%%%%%%
    
    %%
    %compute diff of frames to use only continuous lines.
    %turn classes into a string. 
    %cs=int2str(classes');
    cs=strrep(cs,' ',''); %delete spaces
    %fdiff=diff(frames);
    %inds=find(fdiff~=1);
    %inds=inds+1; %at these indeces we want to put them on a new line.
    %%%%%%%%%%%%%%%%%%%%%%BUG FIX%%%%%%%%%%%%%%%%%%%%%%%%
    %HUGE BUG: In the worst circumstances, sometimes you can have cell A and
    %cell A+1 where Cell A disappears at frame F, and cell A+1 appears at frame
    %F+1. Ex in 25Aligned data. It makes it look like the data is continuous,
    %but a huge leap is taken. 
    Ndiff = diff(NumCells);
    inds = find(Ndiff~=0); %was indsN
    inds = inds+1; %was indsN

    %inds = [inds;indsN];
    inds = unique(inds);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lx = cell(1,length(inds)+1);
    ly = cell(1,length(inds)+1);
    %p_vecs_before = cell(1,length(inds)+1);
    %p_vecs_after = cell(1,length(inds)+1);
    framescheck=cell(1,length(inds)+1);
    numcellcheck = cell(1,length(inds)+1);
    charcellarray=cell(1,length(inds)+1);
    classes = cell(1,length(inds)+1);

    for i=1:length(inds)+1
        if i==1
            charcellarray{i}=cs(1:inds(i)-1);
            framescheck{i}=frames(1:inds(i)-1);
            numcellcheck{i}=NumCells(1:inds(i)-1);
            lx{i} = locsx(1:inds(i)-1);
            ly{i} = locsy(1:inds(i)-1);
            %p_vecs_before{i} = vecs_before(1:inds(i)-1,:);
            %p_vecs_after{i} = vecs_after(1:inds(i)-1,:);
            classes{i} = Classifications(1:inds(i)-1)';
        elseif i==length(inds)+1
            charcellarray{i}=cs(inds(i-1):end);
            framescheck{i}=frames(inds(i-1):end);
            numcellcheck{i}=NumCells(inds(i-1):end);
            lx{i} = locsx(inds(i-1):end);
            ly{i} = locsy(inds(i-1):end);
            %p_vecs_before{i} = vecs_before(inds(i-1):end,:);
            %p_vecs_after{i} = vecs_after(inds(i-1):end,:);
            classes{i} = Classifications(inds(i-1):end)';
        else
            charcellarray{i}=cs(inds(i-1):inds(i)-1);
            framescheck{i}=frames(inds(i-1):inds(i)-1);
            numcellcheck{i}=NumCells(inds(i-1):inds(i)-1);
            lx{i} = locsx(inds(i-1):inds(i)-1);
            ly{i} = locsy(inds(i-1):inds(i)-1);
            %p_vecs_before{i} = vecs_before(inds(i-1):inds(i)-1,:);
            %p_vecs_after{i} = vecs_after(inds(i-1):inds(i)-1,:);
            classes{i} = Classifications(inds(i-1):inds(i)-1)';
        end
    end


    %% ERROR SEARCH
    for i=1:length(numcellcheck)
        if length(unique(numcellcheck{i}))>1 %the cell numbers for the trajectory should all be the same.
            disp('Cell number error')
        end
        %if sum(diff(framescheck{i})>1)>0 %the trajectories should all be continuous
        %    disp('Frame error')
        %end
    end

end
