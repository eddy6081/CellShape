function Make_MSD_Sim_Plots(MSD_sim_dyn, MSD_sim_dyn_error, MSD_sim_stat, MSD_sim_stat_error, MSD_All, MSD_All_error)
    %{
    Author: Chris Eddy
    Date: 03/18/21
    Purpose: Show plots of simulation mean square displacement in
             comparison to MSD of data

    Input Args:

    Outputs:
        plots of MSD simulation in comparison to data.
    %}
    %% MSD Plots
    figure(80)%(1)
    leginfo={};
    msd_colors = {[0,1,0],[1,0,0],[0,0,1]};
    ttls_thisplot={'Stat','Dyn','RT'};
    count=0;
    handles = [];
    
    %Static simulation
    count=count+1;
    x=[0:length(MSD_sim_stat)-1];
    y=MSD_sim_stat;
    yuncert = MSD_sim_stat_error;

    CI_low = y - yuncert;
    CI_high = y + yuncert;
    CI_high = fliplr(CI_high);
    CI=[CI_low,CI_high];

    x_axis = [0:length(y)-1];
    x_plot = [x_axis, fliplr(x_axis)];

    h=fill(x_plot, CI, 1, 'facecolor', msd_colors{count}, 'edgecolor', 'none', 'facealpha', 0.4);
    handles = [handles,h(1)];
    hold on
    scatter(x,y,100,msd_colors{count},'o','filled')
    leginfo{end+1}=ttls_thisplot{count};
    hold on
    
    %Dynamic simulation
    count=count+1;
    x=[0:length(MSD_sim_dyn)-1];
    y=MSD_sim_dyn;
    yuncert = MSD_sim_dyn_error;

    CI_low = y - yuncert;
    CI_high = y + yuncert;
    CI_high = fliplr(CI_high);
    CI=[CI_low,CI_high];

    x_axis = [0:length(y)-1];
    x_plot = [x_axis, fliplr(x_axis)];

    h=fill(x_plot, CI, 1, 'facecolor', msd_colors{count}, 'edgecolor', 'none', 'facealpha', 0.4);
    handles = [handles,h(1)];
    hold on
    scatter(x,y,100,msd_colors{count},'o','filled')
    leginfo{end+1}=ttls_thisplot{count};
    hold on
    
    %Data MSD
    count=count+1;
    x=[0:length(MSD_All)];
    y=[0,MSD_All];
    yuncert = [0,MSD_All_error];

    CI_low = y - yuncert;
    CI_high = y + yuncert;
    CI_high = fliplr(CI_high);
    CI=[CI_low,CI_high];

    x_axis = [0:length(y)-1];
    x_plot = [x_axis, fliplr(x_axis)];

    h=fill(x_plot, CI, 1, 'facecolor', msd_colors{count}, 'edgecolor', 'none', 'facealpha', 0.4);
    handles = [handles,h(1)];
    hold on
    scatter(x,y,100,msd_colors{count},'o','filled')
    leginfo{end+1}=ttls_thisplot{count};
    
    [~,hobj,~,~]=legend(handles,leginfo,'FontSize',20, 'Location', 'northwest');

    xlabel('tau')
    ylabel('MSD (\mum^2)')
    xlim([0,25])

    set(gca, 'FontName','Helvetica','FontSize',28,...
        'YMinorTick', 'on', 'YTick',0:100:MSD_All(end), 'TickDir','out',...
        'YGrid','on','XMinorTick', 'off', 'XTick', 0:2:length(MSD_All), 'LineWidth',5)

end