function Make_MSD_AMME_Plots(classes, MSD_AMME, MSD_AMME_error, MSD_All, MSD_All_error)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Display plots of AM, ME data 

    Inputs Args:
        %dt = time component from which steps were calculated. Use for
        %     simulation of MSD.
        
        classes = cell array, shape 1 x NumCells. Classifications of
                  trajectory i in classes{i}
        
        MSD_AMME = cell array, shape 1 x 3. MSD{1} is AM, MSD{2} is ME, and
                   MSD{3} is IM. 

        MSD_AMME_error = cell array, shape 1 x 3. Contains standard
                         deviation from the mean.

        MSD_All = double, shape 1 x tau. Contains MSD of all trajectories.
        
        MSD_All_error = double, shape 1 x tau. Contains MSD errors.

    Outputs:
        Creates plots of MSD. 
    %}
    figure(70)%(1)
    leginfo={};
    msd_colors = {[0,1,0],[1,0,0],[0,1,1],[0,0,1]};
    ttls_thisplot={'AM','ME','IM','RT'};
    count=0;
    handles = [];
    
    %Calculate fractions of cell phenotypes
    classes_flat = vertcat(classes{:});
    C = max(classes_flat);
    P_state = zeros(1,C); %calculate observed fractions of classes
    for c = 1:C
        P_state(c) = sum(classes_flat==c)/length(classes_flat);
    end
    %convert to AM, ME, IM states
    P_state(1) = P_state(1)+P_state(3);
    P_state(2) = P_state(2)+P_state(4);
    if length(P_state)==5
        P_state(3) = P_state(5);
    elseif length(P_state)==4
        P_state(3) = 0;
    else
        disp("Not sure on how many classes")
    end
    P_state = P_state(1:3);
    
    %Calculate Maximum MSD value. 
    
    
    for i=1:length(MSD_AMME)
        count=count+1;
        if i~=3
            x=[0:length(MSD_AMME{i})];
            y=[0,MSD_AMME{i}];
            yuncert = [0,MSD_AMME_error{i}];

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
        end
    end
    hold on
    
    %Add the data
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

    %Add a weighted average of Amoeboid and Mesenchymal.
    MSD_WA = zeros(1,max([length(MSD_AMME{1}),length(MSD_AMME{2}), length(MSD_AMME{3})]));
    for i=1:length(MSD_WA)
        if length(MSD_AMME{1})>=i
            MSD_WA(i) = MSD_WA(i) + P_state(1)*MSD_AMME{1}(i);
        end
        if length(MSD_AMME{2})>=i
            MSD_WA(i) = MSD_WA(i) + P_state(2)*MSD_AMME{2}(i);
        end
        if length(MSD_AMME{3})>=i
            MSD_WA(i) = MSD_WA(i) + P_state(3)*MSD_AMME{3}(i);
        end
    end

    x=[0:length(MSD_WA)];
    y=[0,MSD_WA];
    %find which is longer.
    E1 = MSD_AMME_error{1};
    E2 = MSD_AMME_error{2};
    E3 = MSD_AMME_error{3};
    ml=max([length(E1),length(E2), length(E3)]);
    if length(E1)<ml
        E1 = [E1,zeros(1,ml-length(E1))];
    end
    if length(E2)<ml
        E2 = [E2,zeros(1,ml-length(E2))];
    end
    if length(E3)<ml
        E3 = [E3,zeros(1,ml-length(E3))];
    end

    yuncert = [0,sqrt(P_state(1)*(E1.^2) + P_state(2)*(E2.^2) + P_state(3)*(E3.^2))];
    %yuncert = [0,sqrt((P_state(1)*(E1.^2) / (P_state(1)+P_state(2))) + (P_state(2)*(E2.^2) / (P_state(1)+P_state(2))))];


    CI_low = y - yuncert;
    CI_high = y + yuncert;
    CI_high = fliplr(CI_high);
    CI=[CI_low,CI_high];

    x_axis = [0:length(y)-1];
    x_plot = [x_axis, fliplr(x_axis)];

    %h=fill(x_plot, CI, 1, 'facecolor', [0,0,0], 'edgecolor', 'none', 'facealpha', 0.4);
    h = errorbar(x, y, yuncert,'vertical','ko');
    h.Marker = 'none';
    h.LineWidth = 2;
    handles = [handles,h(1)];
    hold on
    p = scatter(x,y,100,[0,0,0],'o','filled');
    leginfo{end+1}='WA';

    [~,hobj,~,~]=legend(handles,leginfo,'FontSize',20, 'Location', 'northwest');

    xlabel('tau')
    ylabel('MSD (\mum^2)')
    xlim([0,10.2])

    set(gca, 'FontName','Helvetica','FontSize',28,...
        'YMinorTick', 'on', 'YTick',0:100:MSD_All(end), 'TickDir','out',...
        'YGrid','on','XMinorTick', 'off', 'XTick', 0:2:length(MSD_All), 'LineWidth',5)
    ylim([0,300])
end