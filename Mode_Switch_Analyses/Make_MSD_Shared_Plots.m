function Make_MSD_Shared_Plots(MSD_shape, MSD_shape_error, MSD_real, MSD_real_error)
    %{
    Author: Chris Eddy
    Date: 05/03/22
    Purpose: Display plots of shape-space and real-space MSD on one plot
    with shared x axes.

    Inputs Args:
        
        MSD_shape = double, shape 1 x tau. Contains MSD of all trajectories.
        
        MSD_shape_error = double, shape 1 x tau. Contains MSD errors.

        MSD_real = double, shape 1 x tau. Contains MSD of all trajectories.
        
        MSD_real_error = double, shape 1 x tau. Contains MSD errors.

    Outputs:
        Creates plots of MSD. 
    %}
    figure()
    leginfo={};
    handles = [];
    
        %%%%%%%%%%%%%%%%%%%%
    %PLOT REAL SPACE!
    %%%%%%%%%%%%%%%%%%%%
    
    x=[0:length(MSD_real)];
    y=[0,MSD_real];
    yuncert = [0,MSD_real_error];

    CI_low = y - yuncert;
    CI_high = y + yuncert;
    CI_high = fliplr(CI_high);
    CI=[CI_low,CI_high];

    x_axis = [0:length(y)-1];
    x_plot = [x_axis, fliplr(x_axis)];
    
    yyaxis left
    h=fill(x_plot, CI, 1, 'facecolor', [1,0,0], 'edgecolor', 'none', 'facealpha', 0.4);
    handles = [handles,h(1)];
    hold on
    scatter(x,y,100,[1,0,0],'o','filled')
    leginfo{end+1}='Real';
    ylim([0,600])
    ylabel('MSD (\mum^2)')
    
    %%%%%%%%%%%%%%%%%%%%
    %PLOT SHAPE SPACE!
    %%%%%%%%%%%%%%%%%%%%
    
    x=[0:length(MSD_shape)];
    y=[0,MSD_shape];
    yuncert = [0,MSD_shape_error];

    CI_low = y - yuncert;
    CI_high = y + yuncert;
    CI_high = fliplr(CI_high);
    CI=[CI_low,CI_high];
    
    yyaxis right
    h=fill(x_plot, CI, 1, 'facecolor', [0,0,1], 'edgecolor', 'none', 'facealpha', 0.4);
    handles = [handles,h(1)];
    hold on
    scatter(x,y,100,[0,0,1],'o','filled')
    leginfo{end+1}='Shape';
    ylim([0,15])
    ylabel('MSD')
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    [~,hobj,~,~]=legend(handles,leginfo,'FontSize',20, 'Location', 'northwest');

    xlabel('tau')
    
    xlim([0,20.2])

    %     set(gca, 'FontName','Helvetica','FontSize',28,...
    %         'YMinorTick', 'on', 'YTick',0:100:MSD_All(end), 'TickDir','out',...
    %         'YGrid','on','XMinorTick', 'off', 'XTick', 0:2:length(MSD_All), 'LineWidth',5)
    ax = gca;
    ax.YAxis(1).Color = [1,0,0];
    ax.YAxis(2).Color = [0,0,1];
    set(gca, 'FontName','Helvetica','FontSize',28,...
        'TickDir','out',...
        'XMinorTick', 'off', 'XTick', 0:2:length(MSD_real), 'LineWidth',5)
end