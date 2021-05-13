function [PosX_dyn, PosY_dyn, PosX_stat, PosY_stat]=Simulate_Steps_Persist_v2(P_trans, P_state, dR_params, theta_cells, dt, DT, IM_number, past)
    %{
    Purpose: This simulation works by defining an initial persistence
    direction that remains constant during the entire trajectory. 
    %}
    %% Simulate a persistent (random walk with two states)
    %simulate N cells, length L trajectories
    N = 1000;
    L = 100 / dt;

    %IMPORTANT NOTE: Note that the previous persistent step parameters were
    %obtained for time steps of dt. Therefore, this
    %simulation assumes dt time steps. 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %simulate classes
    
    %parameters
    sig = dR_params(:,:,2);
    mu = dR_params(:,:,1);
    
    %simulation
    %generate N numbers.
    rn = rand(N,1);
    num_end = min(IM_number-1,length(P_state)); %if we are excluding IM. 
    Psum = cumsum(P_state(1:num_end) / sum(P_state(1:num_end))); 
    C_sim = zeros(N,L);
    for traj=1:N
        C_sim(traj,1) = find(rn(traj)<=Psum, 1);
    end
    %previously, used randsample. But this is in the ML library now. 
    C_sim2 = repmat(C_sim(:,1),1,L);


    %%%% For each trajectory, generate a persistent direction (randomly at
    %%%% first)
    vec_simx = rand(N,1);
    vec_simy = rand(N,1);

    %also randomly assign first signs.
    rx_sgn = rand(N,1);
    ry_sgn = rand(N,1);
    vec_simx(rx_sgn<0.5)=-vec_simx(rx_sgn<0.5);
    vec_simy(ry_sgn<0.5)=-vec_simy(ry_sgn<0.5);

    % normalize the vectors.
    vec_sim = [vec_simx,vec_simy];
    
    vec_sim = vec_sim./sqrt(sum(vec_sim.^2,2)); %this makes it so the vector's MAGNITUDE is 1.
    
    %vec_sim = vec_sim./repmat(vecnorm(vec_sim,1,2),1,2); %this makes it so the vector ADDS to 1. 

    %so simulate the classes

    for t = 2:L
        %generate a random number between 0 and 1 for each cell
        r_num = rand(N,1);
        %find the index it corresponds to. 
        individual_Trans = cumsum(P_trans(C_sim(:,t-1),:),2);
        for traj = 1:N
            C_sim(traj,t) = find(r_num(traj)<=individual_Trans(traj,:), 1);%find(individual_Trans(traj,:)>=r_num(traj), 1);
        end
    end
    
    %determine the ecdf of angles using theta_p_AMME_cells %check that it exists.
    ecdf_theta = cell([size(P_trans),2]);
    for p1=1:size(P_trans,1)
        for p2=1:size(P_trans,2)
            %inds = abs(theta_cells{p1,p2})<pi/2;
            [ecdf_theta{p1,p2,1},ecdf_theta{p1,p2,2}]=ecdf(theta_cells{p1,p2});
        end
    end

    %from the classes, simulate the angles for each type of transition
    Angles_dyn = zeros(N,L-1); %really it should be L-1
    Angles_stat = zeros(N,L-1);

    for traj=1:N
        for t=2:L
            %generate a random number.
            rn = rand();
            %generate a second random number.
            rn2 = rand();
            
            %find where rn is in the ecdf of theta. This is angle AWAY from
            %persistent direction.
            %dynamic sim
            ind = find(rn<=ecdf_theta{C_sim(traj,t-1),C_sim(traj,t),1},1);
            angle_p = ecdf_theta{C_sim(traj,t-1),C_sim(traj,t),2}(ind);
            Angles_dyn(traj,t)= angle_p;
            %static sim
            ind = find(rn2<=ecdf_theta{C_sim2(traj,t-1),C_sim2(traj,t),1},1);
            angle_p = ecdf_theta{C_sim2(traj,t-1),C_sim2(traj,t),2}(ind);
            Angles_stat(traj,t)=angle_p;%angle_step + angle_p;
        end
    end
    
    %this is done in the simulation.
    %     Angles_dyn = cumsum(Angles_dyn, 2); %if the angle exceeds 2pi, subtract until not. 
    %     inds = find(Angles_dyn>2*pi);
    %     ang_ratio = floor(Angles_dyn(inds) / (2*pi));
    %     Angles_dyn(inds) = Anles_dyn(inds) - (ang_ratio * 2 * pi);
    
    %     Angles_stat = cumsum(Angles_stat, 2); %if the angle exceeds 2pi, subtract until not. 
    %     inds = find(Angles_stat>2*pi);
    %     ang_ratio = floor(Angles_stat(inds) / (2*pi));
    %     Angles_stat(inds) = Anles_stat(inds) - (ang_ratio * 2 * pi);


    %% SIMULATE STEPS FOR DYNAMIC SIMULATION
    %simulate steps depending on those classes.

    %use the classes to decide how big of a step to take. 
    temp=[reshape(C_sim(:,1:end-1)',1,[]);reshape(C_sim(:,2:end)',1,[])];
    lindex = sub2ind(size(mu),temp(1,:),temp(2,:));
    lindex = reshape(lindex,[],N);
    lindex = lindex';

    steps=lognrnd(mu(lindex), sig(lindex));

    % Angle is negative if to the left of the persistent direction, and
    % positive if it is to the right of the persistent direction.


    PosX_dyn=zeros(N,L);
    PosY_dyn=zeros(N,L);
    phi_dyn = zeros(N,L);
    for traj=1:N
        theta_p = acosfull(vec_sim(traj,1),vec_sim(traj,2));
        for t=2:L
            theta_step = Angles_dyn(traj,t-1); %this is the angle step AWAY from persistence angle.
            phi = theta_p - theta_step;
            phi_dyn(traj,t) = phi;
            PosX_dyn(traj,t) = steps(traj,t-1)*cos(phi);
            PosY_dyn(traj,t) = steps(traj,t-1)*sin(phi);
        end
    end
    PosX_dyn = cumsum(PosX_dyn,2);
    PosY_dyn = cumsum(PosY_dyn,2);


    %% SIMULATE STEPS FOR Static SIMULATION
    %simulate steps depending on those classes.

    %use the classes to decide how big of a step to take. 
    temp=[reshape(C_sim2(:,1:end-1)',1,[]);reshape(C_sim2(:,2:end)',1,[])];
    lindex = sub2ind(size(mu),temp(1,:),temp(2,:));
    lindex = reshape(lindex,[],N);
    lindex = lindex';

    steps=lognrnd(mu(lindex), sig(lindex));


    PosX_stat=zeros(N,L);
    PosY_stat=zeros(N,L);
    for traj=1:N
        theta_p = acosfull(vec_sim(traj,1),vec_sim(traj,2));
        for t=2:L
            theta_step = Angles_stat(traj,t-1);
            phi = theta_p - theta_step;
            %while phi>2*pi
            %    phi = phi - (2*pi);
            %end
            PosX_stat(traj,t) = steps(traj,t-1)*cos(phi);
            PosY_stat(traj,t) = steps(traj,t-1)*sin(phi);
        end
    end
    PosX_stat = cumsum(PosX_stat,2);
    PosY_stat = cumsum(PosY_stat,2);
end