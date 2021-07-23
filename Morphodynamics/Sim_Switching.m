function [positions] = Sim_Switching(loc_data, loc_thresh, keep_persist, D_random, D_persist, mu_persist)
    %{
    Input Args:
        loc_data = double, N_cells x N_dim x N_steps
                   data from Constrained_Random_Walk
        
        loc_thresh = double, x value for where invasive state begins 
                     5 = halfway between
                     0 = all states persistent
                     >10 = all states not persistent

        keep_persist = boolean, True if persistent direction should be
                       constant for each cell, False if persistent
                       direction should be determined by the last velocity
                       vector.
                       Currently, kept constant over the entire trajectory.
                       See commented code lines near line 96.
        
        D = diffusion coefficient 
            0 if balistic
            sigma = standard deviation of step size = sqrt(D)/dt

        mu = average step size
    %}
    
    %using xlocation data.
    xlocs = loc_data(:,1,:);
    %xlocs = squeeze(xlocs);
    %if on left side of opening, make random. state = 0
    %if on right side of opening, make persistent. state = 1
    
    positions = zeros(size(loc_data,1),2,size(loc_data,3));
    
    xlocs = repmat(xlocs,1,size(positions,2),1); %reshape xlocs for ease.
    
    state = zeros(size(positions,1),size(positions,2),size(positions,3));
    state(xlocs>=loc_thresh) = 1; %Box width was set to 10. 
    %copy to same shape as positions.
    %state = reshape(repmat(state,1,size(positions,2)),size(positions,1),size(positions,2),size(positions,3));
    state = state>0; %make logical
    
    v = zeros(size(positions,1),size(positions,2),size(positions,3)); %velcity unit vectors
    %define initial velocity vectors anywhere
    %v(:,:,1)=rand(size(positions,1),size(positions,2));
    v_init = rand(size(positions,1),size(positions,2));
    %normalize?
    v_init = v_init ./ repmat(vecnorm(v_init,2,2),1,2);
    %v = v ./ repmat(vecnorm(v,2,2),1,size(v,2),1);
    v(:,:,1) = v_init;
    
    if keep_persist
        %if we wish for the persistent direction not to change
        v = repmat(v_init,1,1,size(positions,3)); 
    end
    
    %random, brownian motion:
    %X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)
    %D = 2; %diffusion coefficient.
    dt=1; %time difference between frames.
    sigma_random = sqrt(D_random*dt); 
    sigma_persist = sqrt(D_persist*dt);
    %simulate motion. let them all start at origin. 
    
    %create random steps for all points.
    steps = zeros(size(positions,1),size(positions,2),size(positions,3));%normrnd(0,sigma,size(loc_data,1),size(loc_data,2),size(loc_data,3)-1);
    
    %for all persistent points, we need to redo steps. It depends on
    %previous velocity unit vector 
    %steps(state(:)) %
    for traj=1:size(positions,1)
        for t=1:size(positions,3)-1
            %if state(traj,1,t)==1, then it was persistent, and we should
            %use the previous velocity unit vector v(traj,:,t) to tell us
            %direction of step.
            if state(traj,1,t)==1
                for dim=1:size(positions,2)
                    steps(traj,dim,t)=normrnd(v(traj,dim,t)*mu_persist,sigma_persist,1,1,1);
                end
                positions(traj,:,t+1) = positions(traj,:,t) + steps(traj,:,t);
            else
                steps(traj,:,t)=normrnd(0,sigma_random,1,size(positions,2),1);
                positions(traj,:,t+1) = positions(traj,:,t) + steps(traj,:,t);
            end
            if ~keep_persist
                %calculate the new velocity unit vector.
                v_temp = positions(traj,:,t+1) - positions(traj,:,t);
                %v_temp = v_temp ./ repmat(vecnorm(v_temp,2,2),1,size(v,2),1);
                v_temp = v_temp / vecnorm(v_temp,2,2);
                if any(isnan(v_temp))
                    %happens if we take a zero step.
                    v(traj,:,t+1) = v(traj,:,t);
                else
                    v(traj,:,t+1) = v_temp;
                end
            %{
            NOTE: One interesting way to calculate persistence is if the
            cell is persistent than it will hold over the persistent
            direction, but if not then it will recalculate the persistent
            direction. This is in contrast to what is done currently, where
            if keep_persist=true, then the same persistent direction is
            kept for the entire trajectory. Either way, this results in
            trajectories that do not closely resemble the undirected
            migration we see. In the method commented here, we get
            interesting levy flight type trajectories.
                
            else %just added CE
                if state(traj,1,t)==1
                    v(traj,:,t+1) = v(traj,:,t);
                else
                    %calculate the new velocity unit vector.
                    v_temp = positions(traj,:,t+1) - positions(traj,:,t);
                    %v_temp = v_temp ./ repmat(vecnorm(v_temp,2,2),1,size(v,2),1);
                    v_temp = v_temp / vecnorm(v_temp,2,2);
                    if any(isnan(v_temp))
                        %happens if we take a zero step.
                        v(traj,:,t+1) = v(traj,:,t);
                    else
                        v(traj,:,t+1) = v_temp;
                    end
                end
            %}
            end
        end
    end   
    
    %     %plot some of the trajectories!
    %     figure
    %     if size(positions,2)==3
    %         for traj=1:min([5,size(positions,1)])
    %             xs = squeeze(positions(traj,1,:));
    %             ys = squeeze(positions(traj,2,:));
    %             zs = squeeze(positions(traj,3,:));
    %             plot3(xs,ys,zs,'-x');
    %             hold on
    %         end
    %     elseif size(positions,2)==2
    %         for traj=1:min([5,size(positions,1)])
    %             xs = squeeze(positions(traj,1,:));
    %             ys = squeeze(positions(traj,2,:));
    %             plot(xs,ys,'-x');
    %             hold on
    %         end
    %     else
    %         disp("Cannot produce plot! Dimensionality too large")
    %     end
end