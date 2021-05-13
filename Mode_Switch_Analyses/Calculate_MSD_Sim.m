function [MSD_sim_dyn, MSD_sim_dyn_error, MSD_sim_stat, MSD_sim_stat_error]=Calculate_MSD_Sim(PosX_dyn, PosY_dyn, PosX_stat, PosY_stat)
    %{
    Author: Chris Eddy
    Date: 03/18/21
    Purpose: Calculate the MSD of simulation data.
    
    Inputs: 
        PosX_dyn:  double, N x L matrix where N is the number of simulated cells
                   and L is the length of simulation. Contains the x-position of
                   trajectory i at frame f in X(i,f) of the dynamic
                   simulation where cells can transition between phenotype
                   states.

        PosY_dyn:  double, N x L matrix where N is the number of simulated cells
                   and L is the length of simulation. Contains the y-position of
                   trajectory i at frame f in X(i,f) of the dynamic
                   simulation where cells can transition between phenotype
                   states.

        PosX_stat: double, N x L matrix where N is the number of simulated cells
                   and L is the length of simulation. Contains the x-position of
                   trajectory i at frame f in X(i,f) of the dynamic
                   simulation where cell phenotype states are fixed.

        PosY_stat: double, N x L matrix where N is the number of simulated cells
                   and L is the length of simulation. Contains the y-position of
                   trajectory i at frame f in X(i,f) of the static
                   simulation where cell phenotype states are fixed.
    
    Outputs:
        MSD_sim_dyn: double, 1 x L. MSD of the dynamic simulation.

        MSD_sim_dyn_error: double, 1 x L. Standard deviation of the mean of
                           squared differences from the dynamic simulation.

        MSD_sim_stat: double, 1 x L. MSD of the static simulation.

        MSD_sim_stat_error: double, 1 x L. Standard deviation of the mean of
                            squared differences from the static simulation.
    %}
    %%MSD Code for simulations which are continuous.
    L = size(PosX_dyn,2);
    MSD_sim_dyn = zeros(1,L);
    MSD_sim_dyn_error = zeros(1,L);
    for tau=1:L-1
        xdiff = PosX_dyn(:,1+tau:L) - PosX_dyn(:,1:L-tau);
        ydiff = PosY_dyn(:,1+tau:L) - PosY_dyn(:,1:L-tau);
        sum_sq_diffs = (xdiff.^2) + (ydiff.^2);
        MSD_sim_dyn(tau+1) = mean(sum_sq_diffs(:));
        MSD_sim_dyn_error(tau+1) = std(sum_sq_diffs(:))/sqrt(length(sum_sq_diffs(:)));
    end
    
    L = size(PosX_stat,2);
    MSD_sim_stat = zeros(1,L);
    MSD_sim_stat_error = zeros(1,L);
    for tau=1:L-1
        xdiff = PosX_stat(:,1+tau:L) - PosX_stat(:,1:L-tau);
        ydiff = PosY_stat(:,1+tau:L) - PosY_stat(:,1:L-tau);
        sum_sq_diffs = (xdiff.^2) + (ydiff.^2);
        MSD_sim_stat(tau+1) = mean(sum_sq_diffs(:));
        MSD_sim_stat_error(tau+1) = std(sum_sq_diffs(:))/sqrt(length(sum_sq_diffs(:)));
    end

end