function [p_vecs]=Compute_Persist_Vectors2(locsx, locsy, frames)
    %[p_vecs_before, p_vecs_after]
    %{
    Author: Chris Eddy
    Date: 03/15/21
    Purpose: Compute_Persist_Vectors computes the persistent vector from
             vectors within a local time range (t-DT-past, ... t-1, t+1,
             ... t+DT+future). Here, we use the entire trajectory from t=1
             to the end to compute just a single persistent vector which is
             then copied for all time points. That is, the returned p_vecs
             here is the same for all frames in the trajectory. 
    Results:
             This leads to rather diffusive processes, meaning that the
             same persistent vector is typically not maintained during the
             entirety of the cell trajectory. This might be expected, as a
             cell doesn't need to remember information about the entirety
             of its trajectory. 

    Input Args:
        locsx: cell array, x locations of cell i in locsx{i}
        locsy: cell array, y locations of cell i in locsy{i}
        frames: cell array, time frame stamps of cell i in frames{i}
        DT: time distance from time t to compute vectors for SVD. Will compute 
            V(t) - V(t-DT), V(t) - V(t-(DT-1)),... V(t) - V(t-1) for past vectors
           V(t+DT) - V(t), V(t+(DT-1)) - V(t),... V(t-1) - V(t) for future vectors
        past: Additional time for past vectors (that is, DT+past)
        future: Additional time for future vectors (that is, DT+future)
        NumCells: cell array, cell numbers. cell i should all be the same cell
                  number in NumCells{i}

    Output Args:
        p_vecs: cell array, persistent and non-persistent direction vectors
                from positions V as 
                V(1+DT+future) - V(1), V(1+(DT-1)+future) - V(1),... V(2) - V(1)
                and computes orthogonal axes using SVD. 
    %}

    %% Calculate Persistent and NonPersistent direction for each trajectory.

    %%for each trajectory, calculate the persistent and non persistent
    %%direction using SVD
    %pull out positions for each trajectory
    %ucells = unique(NumCells);
    ucells = length(locsx);
    %vecs = zeros(length(ucells),4); %px, py, npx, npy
    p_vecs = cell(1, ucells);
    %p_vecs_before = cell(1, ucells);
    %p_vecs_after = cell(1, ucells);
    for ncell = 1:ucells
        p_vecs{ncell}=zeros(length(locsx{ncell}),4); %px, py, npx, npy.
    	%p_vecs_before{ncell}=zeros(length(locsx{ncell}),4); %px, py, npx, npy.
        %p_vecs_after{ncell} = zeros(length(locsx{ncell}),4); %px, py, npx, npy.
    end
    %a persistent vector for every point.

    %This method uses all vectors from time t within time frame +DT
    for cell_i =1:length(locsx)
        cell_frames = frames{cell_i};
        measure = [locsx{cell_i}, locsy{cell_i}];
        inds_t = [2:length(cell_frames)];
        if ~isempty(inds_t)
            df = cell_frames(inds_t) - cell_frames(1);
            M = measure(inds_t,:) - measure(1,:);
            M = M./[sign(df),sign(df)];
            M = M ./ sqrt(sum(M.^2,2));
            inf_inds = isinf(M(:));
            M(inf_inds)=0.0; %the vector was zero initially.
            nan_inds = isnan(M(:));
            M(nan_inds)=0.0; %the vector was zero initially.
            if size(M,1)>1
                V = M' * M;
                [v, d]= eig(V');
                [~,idx]=sort(diag(d),'descend');
                d=d(idx,idx);
                v = v(:,idx); %left of persistence is positive NP direction.
                v(:,2)=-v(:,2); %right is now positive NP direction.

                %sometimes, this gets the sign wrong too. SVD is
                %usually opposite. I think the thing to do is to take a
                %majority vote for directions.

                %simple test. Dot product. 
                sn = sign(sum(sign(M * v(:,1))));
                if sn==0
                    sn=1;
                end
                v(:,1) = v(:,1)*sn;
                v(:,2) = v(:,2)*sn;

            else
                %easiest, just use the vector as the direction. 
                v = [M(1,:), M(1,2), -M(1,1)];
                %sometimes, the above code places P in the wrong direction.
                %SVD does not. 
                %[~,~,v]=svd(M);
                %it seems we get some weird behaviors if, for example, all
                %the vectors are positive or negative. try:
                %T = [1,0 ; 1,0]; %both vectors point right. 
                %[~,~,V]=svd(T); %V has opposite signs. 
                %Even weirder if you try T=[1,0;0.7,0;-.1,1]
                %I would expect the svd vector to point right and up. 
                %svd(T) points left and up.
                %svd(-T) points right and down
            end
            %p_vecs_before{cell_i}(t,:) = reshape(v,1,[]);
            %v(:,1) = v(:,1) / norm(v(:,1));
            %v(:,2) = v(:,2) / norm(v(:,2));
            p_vecs{cell_i}(:,:) = repmat(reshape(v,1,[]),length(cell_frames),1);
        else
            p_vecs{cell_i}(:,:) = NaN(1,4);
        end
    end
end