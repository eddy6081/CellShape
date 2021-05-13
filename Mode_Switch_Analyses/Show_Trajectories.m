function Show_Trajectories(X, Y, M, title_str)
    %{
    Author: Chris Eddy
    Date: 03/18/21
    Purpose: Display several trajectories of cell migration from simulation

    Inputs: 
        X:  double, N x L matrix where N is the number of simulated cells
            and L is the length of simulation. Contains the x-position of
            trajectory i at frame f in X(i,f).

        Y:  double, N x L matrix where N is the number of simulated cells
            and L is the length of simulation. Contains the y-position of
            trajectory i at frame f in X(i,f).

        M: double, specify the number of trajectories you wish to see.
        
        title_str: string, title of graph.

    Outputs: 
        creates plot of M random trajectories from simulation.
    %}
    N = size(X,1);
    idx = randperm(N);
    figure
    for i=1:M
        plot(X(idx(i),:), Y(idx(i),:), '-');
        hold on
    end
    title(title_str)
end