function [Persist_Steps, Deflect_Steps, Persist_AMME_Steps, Deflect_AMME_Steps]=Compute_Persist_Steps(dt, locsx, locsy, frames, classes, p_vec)
    %{
    AUTHOR: Chris Eddy
    Date: 03/16/21
    Input Arguments:

        locsx:   cell array, 1 x NumCells. Contains the x positions of cell
                 i in locsx{i}

        locsy:   cell array, 1 x NumCells. Contains the y positions of cell
                 i in locsy{i}

        frames:  cell array, 1 x NumCells. Contains the frame numbers of cell
                 i in frames{i}
        
        classes: cell array, 1 x NumCells. Contains the classifications of
                 cell i in classes{i}

    Outputs:

        Persist_Steps: step sizes for each migration mode in the persistent
                       direction. How are we defining the persistent
                       direction? Do we just use previous steps? Or do we
                       use persistent vector?
        Deflect_Steps: perpendicular step direction from persistent.
                       Positive is to the right. 
    %}

    %%
    Persist_Steps=cell(5,5);
    Deflect_Steps=cell(5,5);

    Persist_AMME_Steps=cell(3,3);
    Deflect_AMME_Steps = cell(3,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%VERY IMPORTANT PARAMETER%%%%%%%%%%%%%
    %dt=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i=1:length(locsx)
        fdiff = diff(frames{i});
        
        %for amoeboid mesenchymal migration, change classes 
        classes_AMME = classes{i};
        classes_AMME(classes_AMME==4)=2;
        classes_AMME(classes_AMME==3)=1;
        classes_AMME(classes_AMME==5)=3; %change intermediate.
        
        for c=1:length(fdiff)-dt
            if sum(isnan(p_vec{i}(c,:)))==0 %make sure a persistent vector was defined
                frame_diff = cumsum(fdiff(c:end)); %compute the frame differences
                ind = find(frame_diff==dt); %should only be one, if it is not empty. 
                %not very computationally efficient. at most, you need to only
                %ever look locally around frame f to f - dt.
                %However, this is fine for now. 
                if ~isempty(ind)
                    diffx = locsx{i}(c+ind) - locsx{i}(c);
                    diffy = locsy{i}(c+ind) - locsy{i}(c);
                    dp = dot([diffx,diffy],[p_vec{i}(c,1), p_vec{i}(c,2)]); %/ norm([diffx,diffy]); %this is normalized....
                    if ~isinf(dp)
                        Persist_Steps{classes{i}(c),classes{i}(c+ind)}(end+1) = dp;
                        Persist_AMME_Steps{classes_AMME(c), classes_AMME(c+ind)}(end+1) = dp;
                    end
                    cp = cross([diffx,diffy,0],[p_vec{i}(c,1), p_vec{i}(c,2),0]); %clockwise is positive
                    %this should be the same as (but opposite sign) 
                    np = dot([diffx,diffy],[p_vec{i}(c,3), p_vec{i}(c,4)]);% / norm([diffx,diffy]); %counter clockwise is positive
                    if ~isinf(cp)
                        Deflect_Steps{classes{i}(c),classes{i}(c+ind)}(end+1) = cp(3);
                        Deflect_AMME_Steps{classes_AMME(c), classes_AMME(c+ind)}(end+1) = cp(3);
                    end
                end
            end
        end
    end
end