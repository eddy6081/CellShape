function [p_vecs]=Compute_Persist_Vectors(locsx, locsy, frames, DT, past, future)
    %{
    Author: Chris Eddy
    Date: 03/15/21
    Purpose: Computes the persistent vector at time t from
             vectors within a local time range (t-DT-past, ... t-1, t+1,
             ... t+DT+future).

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
        p_vecs_before: cell array, persistent and non-persistent direction vectors from 
                       V(t) - V(t-DT-past), V(t) - V(t-(DT-1)-past),... V(t) - V(t-1)
                       calculated from SVD.
        p_vecs_after:  cell array, persistent and non-persistent direction vectors from 
                       V(t+DT+future) - V(t), V(t+(DT-1)+future) - V(t),... V(t+1) - V(t)
                       calculated from SVD.
                    ***Contains the persistent
                    ***and non-persistent vectors (in form [p_x, p_y, np_x, np_y]) for
                    ***each frame. Non persistent vector points left from persistent as
                    ***positive direction.
        p_vecs: cell array, persistent and non-persistent direction vectors
                from positions V as 
                V(t) - V(t-DT-past), V(t) - V(t-(DT-1)-past),... V(t) - V(t-1), 
                V(t+DT+future) - V(t), V(t+(DT-1)+future) - V(t),... V(t+1) - V(t)
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

    %DT = 3;% time points away from t
    X = past; %more memory for past.
    Y = future; %more memory for future;
    %This method uses all vectors from time t within time frame +DT
    for cell_i =1:length(locsx)
        cell_frames = frames{cell_i};
        measure = [locsx{cell_i}, locsy{cell_i}];
        for t=1:length(cell_frames)
            inds_t = find(cell_frames<=(cell_frames(t)+DT+Y) & cell_frames>=(cell_frames(t)-DT-X));
            inds_t(inds_t==t)=[];
            if ~isempty(inds_t)
                %if length(inds_t)>1
                df = cell_frames(inds_t) - cell_frames(t);
                M = measure(inds_t,:) - measure(t,:);
                %M = M./[df,df]; %so, this makes step sizes per frame.
                M = M./[sign(df),sign(df)]; %this is step sizes
                %should we normalize each vector so that the "larger"
                %vectors do not carry more weight?
                M = M ./ sqrt(sum(M.^2,2)); %normalize each vector so they carry the same weight.
                %sometimes, there is no movement.
                inf_inds = isinf(M(:));
                M(inf_inds)=0.0; %the vector was zero initially.
                nan_inds = isnan(M(:));
                M(nan_inds)=0.0; %the vector was zero initially.
                
                %quick note: notice that some of the df is negative. This is
                %actually correct, since when we calculate the difference of
                %the measures, for time frames before t, this will end up with
                %x1-x4, when really, it should be x4-x1. Dividing by the
                %negative df makes this work. 

                %singular value decomposition. 
                %could just use [~,~,V]=svd(M), although the signs are
                %different.
                %That is, SVD points NP to the right, where the code below
                %points NP to the left of persistence.
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
                p_vecs{cell_i}(t,:) = reshape(v,1,[]);
                %else
                %    %df = cell_frames(inds_t) - cell_frames(t);
                %    %M = measure(inds_t,:) - measure(t,:);
                %    %M = M./[df,df];
                %    %p = M/norm(M);
                %    %%need nonpersistent direction perpedicular to this
                %    %%vector defining a plane %corresponding to a 90deg
                %    %%clockwise turn.
                %    %%np = [p(2), -p(1)];
                %    %%we want a counter clockwise turn. 
                %    %np = [-p(2), p(1)];
                %    %
                %    %p_vecs_before{cell_i}(t,:) = [p, np];
                %    p_vecs_before{cell_i}(t,:) = NaN(1,4);
                %end
            else
                p_vecs{cell_i}(t,:) = NaN(1,4);
                %p_vecs_before{cell_i}(t,:) = NaN(1,4);
            end
        end
    end

    %     %%%%%DOING AFTER PERSISTENCE VECTOR
    %     
    %     for cell_i =1:length(locsx)
    %         cell_frames = frames{cell_i};
    %         measure = [locsx{cell_i}, locsy{cell_i}];
    %         for t=1:length(cell_frames)
    %             inds_t = find(cell_frames<=(cell_frames(t)+DT+Y) & cell_frames>=cell_frames(t));
    %             inds_t(inds_t==t)=[];
    %             if ~isempty(inds_t)
    %                 %if length(inds_t)>1
    %                 df = cell_frames(inds_t) - cell_frames(t);
    %                 M = measure(inds_t,:) - measure(t,:);
    %                 M = M./[df,df];
    %                 %quick note: notice that some of the df is negative. This is
    %                 %actually correct, since when we calculate the difference of
    %                 %the measures, for time frames before t, this will end up with
    %                 %x1-x4, when really, it should be x4-x1. Dividing by the
    %                 %negative df makes this work. 
    %                 %singular value decomposition. 
    %                 %could just use [~,~,V]=svd(M), although the signs are
    %                 %different.
    %                 V = M' * M;
    %                 [v, d]= eig(V');
    %                 [~,idx]=sort(diag(d),'descend');
    %                 d=d(idx,idx);
    %                 v = v(:,idx);
    %                 %NOTE: THIS IS A COUNTER CLOCKWISE TURN - positive
    %                 %non-persistent direction is to the left of the
    %                 %persistent direction.
    %                 p_vecs_after{cell_i}(t,:) = reshape(v,1,[]);
    %                 %else
    %                 %    %df = cell_frames(inds_t) - cell_frames(t);
    %                 %    %M = measure(inds_t,:) - measure(t,:);
    %                 %    %M = M./[df,df];
    %                 %    %p = M/norm(M);
    %                 %    %%need nonpersistent direction perpedicular to this
    %                 %    %%vector defining a plane %corresponding to a 90deg
    %                 %    %%clockwise turn.
    %                 %    %%np = [p(2), -p(1)];
    %                 %    %%we want a counter clockwise turn. 
    %                 %    %np = [-p(2), p(1)];
    %                 %    %
    %                 %    %p_vecs_before{cell_i}(t,:) = [p, np];
    %                 %    p_vecs_after{cell_i}(t,:) = NaN(1,4);
    %                 %end
    %             else
    %                 p_vecs_after{cell_i}(t,:) = NaN(1,4);
    %             end
    %         end
    %     end
end