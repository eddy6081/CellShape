function [P_trans, P_class, D_time] = Compute_Prob_Trans(classes, frames, IM_number)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Calculate Probability of transition per unit time. Each row
    of transitions should be normalized to one. That is, the probability of
    P1 will transition to to P2:5

    Input Args:
        classes : cell array, shape 1 x NumCells. Classifications of
                  trajectory i in classes{i}
        
        frames:  cell array, 1 x NumCells. Contains the frame numbers of cell
                 i in frames{i}

        IM_number: specifies class for intermediate state. If you wish to
                   also see transition rates and dwell times for IM state,
                   make this an arbitrarily large number.

    Outputs:
        P_trans: double, shape phenotype x phenotype. 
                 Probability rates (probability per frame) of transition
                 matrices. Counts a transition such as AE-IM-IM-FP with a
                 transition rate from AE-FP in 1 / 3 frames.

        P_class: double, shape 1 x phenotype
                 Observed population fractions. Probability of observation.
        
        D_time: double, shape phenotype x 1
                Dwell time for phenotypes. Calculated as 1 /
                (1-P_trans(i,i)). Unit of frame time. 
    %}
    %% Compute transitions
    classes_flat = vertcat(classes{:}); %just determining biggest class number.
    [C,~] = max(unique(classes_flat));
    transitions=zeros(C,C);
    time_trans = cell(C,C);
    P_class = zeros(1,C);

    for i=1:length(frames)
        frs = frames{i}; %trajectory ith frames
        cs = classes{i}; %trajectory ith classes
        stop_inds=find(diff(frs)>1); %find continuous segments
        start_inds=[1;stop_inds+1];
        stop_inds(end+1)=length(frs);
        for si = 1:length(stop_inds)
            %we are now working with continuous segments
            this_cs = cs(start_inds(si):stop_inds(si));
            this_frs = frs(start_inds(si):stop_inds(si));
            this_frs(this_cs==IM_number)=[]; %delete IM for correct transition rate
            this_cs(this_cs==IM_number)=[];
            for tr = 1:length(this_cs)-1
                transitions(this_cs(tr),this_cs(tr+1)) = transitions(this_cs(tr),this_cs(tr+1)) + 1; %just record that there was a transition
                time_trans{this_cs(tr),this_cs(tr+1)}(end+1) = this_frs(tr+1) - this_frs(tr); %record how long the transition took. 
            end
        end
    end
    
    P_trans = zeros(C,C);
    for p1 = 1:C
        for p2 = 1:C
            P_trans(p1,p2) = sum(1./time_trans{p1,p2}) / sum(1./[time_trans{p1,:}]);
        end
    end
    %     %not correct, necessarily.
    %     %calculate 
    %     for p1=1:C
    %         for p2=1:C
    %             P_trans(p1,p2) = length(time_trans{p1,p2}) / mean(time_trans{p1,p2},'omitnan');
    %         end
    %     end
    %     %now, make each row sum to one.
    %     P_trans = P_trans ./ sum(P_trans,2,'omitnan');
    
    %the rate here is probability per frame (prob per 15 minutes).
    
    %class probabilities
    for c = 1:C
        P_class(c)=sum(classes_flat==c)/length(classes_flat);
    end
    
    %dwell times from rates.
    D_time = eye(C) ./ (1 - P_trans);
    D_time = diag(D_time);
    
    
    %What about the error?
end