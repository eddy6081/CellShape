function Make_Theta_Plots(theta_steps)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Show plots of theta angle, direction of steps taken from 
             persistent vector.
    
    Input Args: 
        theta_steps: cell array, phenotype x phenotype in size.
                     contains the theta angles, direction of steps taken
                     from persistent vectors.
    
    Outputs:
        Creates graph of theta angle for each phenotype. 
    %}
    %% Generate plot of theta.
    if size(theta_steps,1)==5
        ttls_transition={'AE-AE','AE-FP','AE-BB','AE-LA','AE-IM';...
                         'FP-AE','FP-FP','FP-BB','FP-LA','FP-IM';...
                         'BB-AE','BB-FP','BB-BB','BB-LA','BB-IM';...
                         'LA-AE','LA-FP','LA-BB','LA-LA','LA-IM';...
                         'IM-AE','IM-FP','IM-BB','IM-LA','IM-IM'};
    elseif size(theta_steps,1)==3
        ttls_transition={'AM-AM','AM-ME','AM-IM';...
                         'ME-AM','ME-ME','ME-IM';...
                         'IM-AM','IM-ME','IM-IM'};
    else
        disp("You must define title for transitions")
    end
    
    figure;
    count=0;
    for p1=1:size(theta_steps,1)
        for p2=1:size(theta_steps,2)
            count=count+1;
            subplot(size(theta_steps,1),size(theta_steps,2),count)
            h=histogram(theta_steps{p1,p2},'BinWidth',pi/10);
            xlim([-pi,pi])
            %title(ttls_transition{p1,p2})
            %xlabel(sprintf('$\theta$ between %s',ttls_transition{p1,p2}))
            %ylabel('counts')
            set(gca, 'XTick',-pi:pi:pi,'XTickLabel',{'-\pi','0','\pi'},'LineWidth',2)
            %h.LineWidth=2;
            
        end
    end
    
end