function [positions] = Sim_Switching_v2(loc_data, loc_thresh, D, Pt, Ps)
    %{
    PERSISTENT RANDOM WALK
    Input Args:
        loc_data = double, N_cells x N_dim x N_steps
                   data from Constrained_Random_Walk
        
        loc_thresh = double, x value for where invasive state begins 
                     5 = halfway between
                     0 = all states persistent
                     >10 = all states not persistent
        
        D = diffusion coefficient 
            0 if balistic
            sigma = standard deviation of step size = sqrt(D)/dt

        Pt = Persistance time (relative to time step size) >=1).
             the bigger this value, the more persistent.
             Pt = 1 is a random walk.

        Ps = Persistance speed %recommend less than Pt
             Scaling for random draw sqrt(Ps^3 * dt^2 / Pt)

        So doing, Pt=1, Ps = 1.13 should result in what we had in Sim_switching
        where the mu=1.2. No, because it is a random walk, not directional.

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
    
    %random, brownian motion:
    %X(t + dt) = X(t) + N(0, delta**2 * dt; t, t+dt)
    %D = 2; %diffusion coefficient.
    dt=1; %time difference between frames.
    sig = sqrt(D*dt); 
    %simulate motion. let them all start at origin. 
    
    %create random steps for all points.
    steps = zeros(size(positions,1),size(positions,2),size(positions,3));%normrnd(0,sigma,size(loc_data,1),size(loc_data,2),size(loc_data,3)-1);
    
    %fill the first time step. 
    if sig~=0
        steps(:,:,1) = normrnd(0,sig,size(steps,1),size(steps,2));
    else
        steps(:,:,1) = normrnd(1,sig,size(steps,1),size(steps,2));
    end
    
    %fill t=2 positions
    positions(:,:,2) = positions(:,:,1) + steps(:,:,1);
    %for all persistent points, we need to redo steps. It depends on
    %previous velocity unit vector 
    %steps(state(:)) %
    
    alpha_p = (1-dt/Pt); %unitless
    Fp = sqrt((Ps^2)* (dt^3) / Pt); %units of meters.
    
    if sig==0
        %set alpha_p to 1.
        %only for ballistic trajectory.
        alpha_p = 1;
    end
    
    for traj=1:size(positions,1)
        for t=2:size(positions,3)-1
            %if state(traj,1,t)==1, then it was persistent, and we should
            %use the previous velocity unit vector v(traj,:,t) to tell us
            %direction of step.
            if state(traj,1,t)==1
                for dim=1:size(positions,2)
                    steps(traj,dim,t)=alpha_p*steps(traj,dim,t-1) + Fp*normrnd(0,sig);
                end
                positions(traj,:,t+1) = positions(traj,:,t) + steps(traj,:,t);
            else
                steps(traj,:,t)=normrnd(0,sig,1,size(positions,2),1);
                positions(traj,:,t+1) = positions(traj,:,t) + steps(traj,:,t);
            end
        end
    end   
    
end