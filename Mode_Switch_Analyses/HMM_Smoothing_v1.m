function [U, Z] = HMM_Smoothing_v1(P_trans, probs, CM)
    %{
    Implement hidden markov model = Upgrade to softmax.
    P_trans = transition probabilities
    probs = Observable probabilities - we descritize to a class as the
    observable.
    CM = emission probabilities. 
    find hidden states.
    
    Do not update the transition or dwell time matrix.
    Input Args:
            classes : cell array, shape 1 x NumCells. Classifications of
                      trajectory i in classes{i}
            
            CM : numeric array. Phenotype x Phenotype, no intermediate.
                 Confusion matrix from model training. See All_SVM.py.

    Output Args:
            U : max probability of observing sequence at each step

            Z : true class labels (hidden states)
    %}
    %from probs, get classes.
    classes = cell(1,length(probs));
    for traj=1:length(classes)
        classes{traj}=zeros(size(probs{traj},1),1);
        for frame=1:size(classes{traj},1)
            [~,classes{traj}(frame)] = max(probs{traj}(frame,:)); %softmax
        end 
    end 
    
    classes_flat = vertcat(classes{:}); %just determining biggest class number.
    %class probabilities
    for c = 1:length(unique(classes_flat))
        P_class(c)=sum(classes_flat==c)/length(classes_flat);
    end
    
    %given classes, use HMM
    %Y = observables = SVM classifications = classes 
    %Z = hidden states = True class labels
    Z = cell(1,length(classes));
    U = cell(1,length(classes));
    
    CM(CM<0.001)=0.001; %set minimum error rates
    
    for traj=1:length(classes)
        zs = zeros(length(classes{traj}),1);
        us = zeros(length(classes{traj}),1);
        for frame = 1:length(zs)
            if frame==1
                %compute first element
                zs(frame)=classes{traj}(frame); %set initial class as SVM's output.
                us(frame)=P_class(zs(frame))*CM(zs(frame),zs(frame));
            else
                %P_trans(zs(z_i-1),:) = P( z_k | z_(k-1) )
                %p_i = P(z_k^i | z_(k-1))
                %so now, given that the observed is classes{traj}(z_i)) = x_k
                %P( x_k | z_k) = emission prob. Confusion matrix. given
                %that given true is z_k, prob of observing x_k. so take
                %column of CM = C, with elements c_i = P(x_k | z_k^i).
                %so P(z_k^i | z_(k-1)) * P(x_k | z_k^i) = p_i * c_i
                [us(frame),zs(frame)] = max( P_trans(zs(frame-1),1:size(CM,2))' .* CM(:,classes{traj}(frame)) * us(frame-1)); %argmax 
            end
        end
        Z{traj}=zs;
        U{traj}=us;
    end   
end
    
    
    