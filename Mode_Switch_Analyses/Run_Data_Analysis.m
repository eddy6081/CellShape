function [] = Run_Data_Analysis(DT, use_softmax, delete_IM, past, future)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Put together all necessary functions to run analysis of step sizes and
    fit parameters.

    Results: So far, best results are using DT = 0, past = 6, future = 5.
             This illustrates that the cell integrates information about
             where it intends to go along with information about its
             momentum (where it has been).
             A physical model might include a force in the direction of
             persistence, and momentum calculated from its previous
             locations...
             %these two things together make the persistent vector. 
             Could also add a random force (which leads to weird
             displacements) caused by the gel structure. The magnitude of
             this dictates the displacements away from persistent vector.

    %use this code as a model for how to run the scripts.

    Input args:
        DT: double, This is a specified time step to use in which to calculate the
            persistent vector from SVD. Currently, time step t will use positions 
            at t - DT,... t - 2,  t - 1 and t + 1 , t + 2, ... t + DT as
            vectors from position at time t [past times are corrected to
            point towards time t position] to calculate persistent and
            non-persistent directions using SVD.

        past: double, This is a specified time step used to add additional PAST
              position difference vectors to calculate persistent and
              non-persistent directions using SVD. SVD will use positions
              at t - DT - past, ... t-2, t-1 and t + 1 , t + 2, ... t + DT
              as vectors (see DT above).

        future: double, This is a specified time step used to add additional PAST
                position difference vectors to calculate persistent and
                non-persistent directions using SVD. SVD will use positions
                at t - DT, ... t-2, t-1 and t + 1 , t + 2, ... t + DT + future
                as vectors (see DT above).

                - combine with "past" argument to use positions from 
                  t - DT - past, ... t-2, t-1 and t + 1 , t + 2, ... t + DT + future
                  as vectors in SVD. 

        use_softmax: boolean, Specifies whether classification should use softmax -
                     if true, highest classification probability is used to
                     determine class, regardless of probability value. Will
                     not have "intermediate" states. If false, a threshold
                     of 0.6 for probability classification - otherwise cell
                     is in "intermediate" state.

        delete_IM: boolean, Option during analysis that if "use_softmax" is
                   set to false, IM states can be "deleted" and not used in
                   analysis. This results in non-continuous trajectories
                   but uses data only where classification probability is
                   high.

    Outputs:
        Currently, I wasn't concerned about output, but feel free to put
        arguments in the brackets as necessary. 
    %}
    %%
    %load data
    [lx, ly, framescheck, numcellcheck, classes] = load_trajectories(use_softmax, delete_IM);
    
    ml=0;
    for i = 1:length(lx)
        if size(lx{i},1)>ml
            ml = size(lx{i},1);
        end
    end
    Xtemp = nan(length(lx),ml);
    Ytemp = nan(length(lx),ml);
    for i=1:length(lx)
        Xtemp(i,1:size(lx{i},1)) = lx{i}';
        Xtemp(i,:) = Xtemp(i,:) - Xtemp(i,1);
        Ytemp(i,1:size(lx{i},1)) = ly{i}';
        Ytemp(i,:) = Ytemp(i,:) - Ytemp(i,1);
    end
    Show_Trajectories(Xtemp, Ytemp, 5, 'Data')
    
    %calculate transition rates, observed class fractions, and dwell times
    %pass as argument to Simulation.
    [P_trans, P_class, D_time] = Compute_Prob_Trans(classes, framescheck, 5);
    
    %Calculate transition rates, but for AM, ME, IM states. 
    classes_temp = classes;
    for i=1:length(classes)
        classes_temp{i}(classes_temp{i}==4)=2; %change LA to FP, make ME state
        classes_temp{i}(classes_temp{i}==3)=1; %change BB to AE, make AM state
        classes_temp{i}(classes_temp{i}==5)=3; %change IM to class 3
    end
    [P_AMME_trans, P_AMME_class, D_AMME_time] = Compute_Prob_Trans(classes_temp, framescheck, 3);
    
    %Calculate MSD for data
    [MSD, MSD_error, MSD_AMME, MSD_AMME_error, MSD_All, MSD_All_error] = Calculate_MSD_Data(lx, ly, framescheck, classes);
    
    %Plot the AM-ME state MSD for dwelling types, all experiment data, and
    %a weighted average of the dwelling types. 
    Make_MSD_AMME_Plots(classes, MSD_AMME, MSD_AMME_error, MSD_All, MSD_All_error)

    %Estimate the persistent vectors using SVD
    [p_vecs]=Compute_Persist_Vectors(lx, ly, framescheck, DT, past, future);
    
    %Estimate persistent vectors using ALL vectors in each trajectory, so
    %each frame in a trajectory has the same persistent vector.
    %[p_vecs2]=Compute_Persist_Vectors2(lx, ly, framescheck);


    %Use persistent vectors to determine size of steps taken in the
    %persistent and non-persistent directions.
    [Persist_Steps, Deflect_Steps, Persist_AMME_Steps, Deflect_AMME_Steps]=Compute_Persist_Steps(1, lx, ly, framescheck, classes, p_vecs);

    %Feel free to uncomment these next two lines. I was looking at just the
    %AM ME data (see below) and therefore I didn't want to run these scripts. If you
    %want to know how large the steps taken are for each migration mode,
    %use these scripts. 
    
    %Determine magnitude of step sizes along with directionality away from
    %persistent direction for each migration mode.
    %[dR_steps, dTheta_steps] = Compute_dR_dTheta_steps(Persist_Steps, Deflect_Steps);

    %Use the step sizes to fit log-normal CDF to determine step size
    %parameters.
    %[dR_params_useful, dR_params, dR_params_error] = Fit_CDF_dR_steps(dR_steps);
    
    %{
    Three State Model -
        assign AE, BB to Amoeboid AM - assign LA, FP to Mesenchymal ME
        %switching probabilities are all contained within. 
    %}
    
    %Determine magnitude of step sizes along with directionality away from
    %persistent direction for each AM, ME, and IM states.
    [dR_AMME_steps, dTheta_AMME_steps] = Compute_dR_dTheta_steps(Persist_AMME_Steps, Deflect_AMME_Steps);
    %dTheta steps are the angle away from persistence direction
    
    %Show graphs of the angles of steps away from persistence direction.
    %{
        This gives you an idea of how persistent the migration is, which is
        determined by how you estimated the persistent vectors in line 76.
    %}
    Make_Theta_Plots(dTheta_AMME_steps)
    
    %Fit step sizes to Gaussian and log-normal distributions
    %with the goal of determining which is a better type of fit.
    [Gauss_AMME_params, Gauss_AMME_errors] = Fit_Gaussian_dR_steps(Persist_AMME_Steps);
    [Log_AMME_params, Log_AMME_errors] = Fit_LogNorm_dR_steps(dR_AMME_steps);
    
    %Note: Parameters above are subject to binning issues. The CDF is not.
    
    %Use the step sizes to fit log-normal CDF to determine step size
    %parameters.
    [dR_AMME_params_useful, dR_AMME_params, dR_AMME_params_error] = Fit_CDF_dR_steps(dR_AMME_steps);
    
    %Plot the parameters of the CDF of log-normal. 
    Make_CDF_Param_Plots(dR_AMME_params,dR_AMME_params_error,dR_AMME_params_useful,use_softmax)
    
    %Create simulations from step size parameters, theta angles, transition
    %probabilities. 
    
    %this function has a changing persistent direction that is determined
    %via SVD of position changes from t-DT-past, .. t-2, t-1 to t.
    [PosX_dyn, PosY_dyn, PosX_stat, PosY_stat]=Simulate_Steps_Persist(P_AMME_trans, P_AMME_class, dR_AMME_params, dTheta_AMME_steps, 1, DT, 3, past);
    Show_Trajectories(PosX_dyn, PosY_dyn, 5, 'moving p')
    
    %this function has a stationary persistent direction, randomly assigned
    %to each cell.
    [PosX_dyn2, PosY_dyn2, PosX_stat2, PosY_stat2]=Simulate_Steps_Persist_v2(P_AMME_trans, P_AMME_class, dR_AMME_params, dTheta_AMME_steps, 1, DT, 3, past);
    Show_Trajectories(PosX_dyn2, PosY_dyn2, 5, 'stationary p')
    
    %note: We ought to be able to run the same analysis on our simulated
    %data to get the same information back. Would need to convert to cell
    %arrays as needed, of course. 
    
    %Calculate MSD of these simulations. 
    [MSD_sim_dyn, MSD_sim_dyn_error, MSD_sim_stat, MSD_sim_stat_error]=Calculate_MSD_Sim(PosX_dyn2, PosY_dyn2, PosX_stat2, PosY_stat2);
    
    %Show plots of the simulations compared to actual data. 
    %{
        The experiment data shows a super-diffusive process and therefore,
        persistence plays a part in it.
        If theta plots showed peaks at 0 and +- pi, then you are likely
        looking at a 1-D random walk and therefore may observe a diffusive
        process (linear in MSD). 
        If theta angles showed no peaks but were uniform, then you are
        likely looking at a 2-D random walk, and again observe a diffusive
        process.
        You should be looking for data that dislays theta angles with
        likely a bell shaped distribution focused on 0.
    %}
    Make_MSD_Sim_Plots(MSD_sim_dyn, MSD_sim_dyn_error, MSD_sim_stat, MSD_sim_stat_error, MSD_All, MSD_All_error)
    
end
