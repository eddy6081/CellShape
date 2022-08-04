function Make_CDF_Param_Plots(dR_params,dR_params_error,dR_params_useful,softmax)
    %% Generate bar-plot of the lognormal parameters 

    %Interpreting the parameters from dR_params is sort of tricky. See the
    %following plot:

    % % 
    % % figure
    % % subplot(2,2,1)
    % % h=histogram(lognrnd(0,.5,5000,1));
    % % title('\mu = 0, \sigma = 0.5')
    % % xlim([0,10])
    % % subplot(2,2,2)
    % % h2=histogram(lognrnd(1,0.5,5000,1));
    % % h2.BinWidth=h.BinWidth;
    % % title('\mu = 1, \sigma = 0.5')
    % % xlim([0,10])
    % % subplot(2,2,3)
    % % h3=histogram(lognrnd(0,0.2,5000,1));
    % % h3.BinWidth=h.BinWidth;
    % % title('\mu = 0, \sigma = 0.2')
    % % xlim([0,10])
    % % subplot(2,2,4)
    % % h4=histogram(lognrnd(0,0.8,5000,1));
    % % h4.BinWidth=h.BinWidth;
    % % title('\mu = 0, \sigma = 0.8')
    % % xlim([0,10])

    %a larger mean both shifts the peak and increases the spread of the
    %lognormal distribution. changing sigma widens or thins the curve, as
    %expected.

    %define which parameters to plot
    if softmax
        ph = 2;
    else
        ph = 3;
    end

    %see errorbar
    ttls_transition={'AM-AM','AM-ME','AM-IM';...
                     'ME-AM','ME-ME','ME-IM';...
                     'IM-AM','IM-ME','IM-IM'};
                 
    %{
    We would like to plot the data so that we have dwell first, then
    transitions. That means, we need to reshape the parameters
    %}
    if ~softmax
        inds = [1,4,7;9,3,6;5,2,8];
    else
        inds = [1,4,7;5,2,8;9,3,6];
    end
    ttls_transition = ttls_transition(inds);
    dR_params_useful=dR_params_useful(inds);
    for ch=1:size(dR_params,3)
        dR_params(:,:,ch)=dR_params(inds+((ch-1)*9));
    end
    for ch=1:size(dR_params_error,3)
        dR_params_error(:,:,ch)=dR_params_error(inds+((ch-1)*9));
    end
    
    %%%%%%%%%%SIGMA%%%%%%%%%%%%%
    figure
    %colors
    if softmax
        bar_colors = {[0,1,0],[0,.6,0.6],[1,0,0],[.93,.5,.93]}; %am-am, am-me, me-me, me-am
    else
        bar_colors = {[0,1,0], [0,.6,.6], [0, .6, .6],[0,1,1],[0,0.58,1],[0,0.58,1], [1,0,0],[.5,.5,.5],[.5,.5,.5]};
        %bar_colors = {[0,1,0], [0,.6,.6], [0, .6, .6], [.5,.5,.5],[1,0,0],[.5,.5,.5],[0,0.58,1],[0,0.58,1],[0,1,1]};
        %bar_colors = {[0,1,0], [0,.6,0.6], [0,1,0.5], [.93,.5,.93],  [1,0,0], [0.5,0.5,0.5], [0,0.58,1], [.415,1,0], [0,1,1]}; 
        %              %am-am,   am-me,     am-im,       me-am        me-me,      me-im,        im-am,      im-me,    im-im
    end

    bar_colors=cell2mat(bar_colors');
    %bar graphs of persistance diffusion.
    counter=0;
    yvals=reshape(dR_params(1:ph,1:ph,2)',1,[]);%[reshape(dR_params(1:2,1:2,2)',1,[]),Persist_params_cdf(end,2)];
    er_high = reshape(dR_params_error(1:ph,1:ph,4)',1,[]) - yvals;
    %[reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,4)',1,[]),Persist_params_cdf_error(end,4)]-yvals;
    er_low =  yvals - reshape(dR_params_error(1:ph,1:ph,3)',1,[]);
    %yvals - [reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,3)',1,[]),Persist_params_cdf_error(end,3)];
    b=bar([1:length(yvals)],yvals);
    b.CData = bar_colors;
    b.FaceColor='flat';

    hold on

    %for i=1:size(Persist_params_cdf,1)
    ylabel('\sigma')
    %ylim([0,7])

    set(gca, 'FontName','Helvetica','FontSize',28,...
        'YMinorTick', 'on', 'YTick',0:.2:2, 'TickDir','out',...
        'YGrid','on','XMinorTick', 'off', 'XTick',1:1:length(yvals),'XTickLabel',[reshape(ttls_transition(1:ph,1:ph)',1,[]),{'All'}],...
        'XTickLabelRotation',90,'LineWidth',5)

    er=errorbar([1:length(yvals)],yvals,er_low,er_high);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    er.LineWidth=5;

    %%%%%%%%%%MU%%%%%%%%%%%%%
    figure
    % %colors
    % if softmax
    %     bar_colors = {[0,1,0], [0,.6,0.6], [0,1,0.5], [.93,.5,.93],  [1,0,0], [0.5,0.5,0.5], [0,0.58,1], [.415,1,0], [0,1,1]}; 
    %     %              %am-am,   am-me,     am-im,       me-am        me-me,      me-im,        im-am,      im-me,    im-im
    % else
    %     bar_colors = {[0,1,0],[0,.6,0.6],[.93,.5,.93],[1,0,0]}; %am-am, am-me, me-me, me-am
    % end
    %     
    % bar_colors=cell2mat(bar_colors');
    %bar graphs of persistance diffusion.
    counter=0;
    yvals=reshape(dR_params(1:ph,1:ph,1)',1,[]);%[reshape(dR_params(1:2,1:2,2)',1,[]),Persist_params_cdf(end,2)];
    er_high = reshape(dR_params_error(1:ph,1:ph,2)',1,[]) - yvals;
    %[reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,4)',1,[]),Persist_params_cdf_error(end,4)]-yvals;
    er_low =  yvals - reshape(dR_params_error(1:ph,1:ph,1)',1,[]);
    %yvals - [reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,3)',1,[]),Persist_params_cdf_error(end,3)];
    b=bar([1:length(yvals)],yvals);
    b.CData = bar_colors;
    b.FaceColor='flat';

    hold on

    %for i=1:size(Persist_params_cdf,1)
    ylabel('\mu')
    %ylim([0,7])

    set(gca, 'FontName','Helvetica','FontSize',28,...
        'YMinorTick', 'on', 'YTick',0:.2:2, 'TickDir','out',...
        'YGrid','on','XMinorTick', 'off', 'XTick',1:1:length(yvals),'XTickLabel',[reshape(ttls_transition(1:ph,1:ph)',1,[]),{'All'}],...
        'XTickLabelRotation',90,'LineWidth',5)

    er=errorbar([1:length(yvals)],yvals,er_low,er_high);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    er.LineWidth=5;


    %% Generate bar-plot of the lognormal transformed parameters.

    %first calculate error on parameters using error propogation.
    %uncertainty in the mean back-calculated mean:
    %pull out uncertainty from necessary variables.
    mu_err = reshape(dR_params_error(1:ph,1:ph,1)',[],1);
    mus=reshape(dR_params(1:ph,1:ph,1)',[],1);
    mu_err = mus - mu_err; %the error in the 95 CI is symmetric. 

    sig_err = reshape(dR_params_error(1:ph,1:ph,3)',[],1);
    sigs=reshape(dR_params(1:ph,1:ph,2)',[],1);
    sig_err = sigs - sig_err; %the error in the 95 CI is symmetric. 

    %now, pull out the calculated mean and variance of the distribution.
    params=cell2mat(reshape(dR_params_useful(1:ph,1:ph)',[],1));
    m = params(:,1);
    v = params(:,2);

    %Now, calculate the erros. 
    %these were done by hand, could be good to verify. 

    m_err = sqrt( ((mu_err.^2).* (m.^2)) + ((sig_err.^2) .* (sigs.^2) .* (m.^2)) );

    v_err = sqrt( ((2*v).^2) .* ( (mu_err.^2) + ((sig_err.^2).*(sigs.^2).*((1 + (1./(1-exp(-(sigs.^2))))).^2)) ) );



    %%%%%%%%%%SIGMA%%%%%%%%%%%%%
    figure
    %colors
    %bar_colors = {[0,1,0],[0,.6,0.6],[.93,.5,.93],[1,0,0]};%,[0,0,1]}; %am-am, am-me, me-me, me-am, all
    %bar_colors=cell2mat(bar_colors');
    %bar graphs of persistance diffusion.
    counter=0;
    yvals=v;
    er_high = v_err;
    %[reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,4)',1,[]),Persist_params_cdf_error(end,4)]-yvals;
    er_low =  v_err;
    %yvals - [reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,3)',1,[]),Persist_params_cdf_error(end,3)];
    b=bar([1:length(yvals)],yvals);
    b.CData = bar_colors;
    b.FaceColor='flat';

    hold on

    %for i=1:size(Persist_params_cdf,1)
    ylabel('\sigma^2 (\mum^2)')
    %ylim([0,7])

    set(gca, 'FontName','Helvetica','FontSize',28,...
        'YMinorTick', 'on', 'YTick',0:3:15, 'TickDir','out',...
        'YGrid','on','XMinorTick', 'off', 'XTick',1:1:length(yvals),'XTickLabel',[reshape(ttls_transition(1:ph,1:ph)',1,[]),{'All'}],...
        'XTickLabelRotation',90,'LineWidth',5)

    er=errorbar([1:length(yvals)],yvals,er_low,er_high);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    er.LineWidth=5;

    %%%%%%%%%%MU%%%%%%%%%%%%%
    figure
    
    %bar graphs of persistance diffusion.
    counter=0;
    yvals=m;
    er_high = m_err;
    %[reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,4)',1,[]),Persist_params_cdf_error(end,4)]-yvals;
    er_low =  m_err;
    %yvals - [reshape(Persist_AMME_Transition_params_cdf_error(1:2,1:2,3)',1,[]),Persist_params_cdf_error(end,3)];
    b=bar([1:length(yvals)],yvals);
    b.CData = bar_colors;
    b.FaceColor='flat';

    hold on

    %for i=1:size(Persist_params_cdf,1)
    ylabel('m (\mum)')
    %ylim([0,7])

    set(gca, 'FontName','Helvetica','FontSize',28,...
        'YMinorTick', 'on', 'YTick',0:1:4, 'TickDir','out',...
        'YGrid','on','XMinorTick', 'off', 'XTick',1:1:length(yvals),'XTickLabel',[reshape(ttls_transition(1:ph,1:ph)',1,[]),{'All'}],...
        'XTickLabelRotation',90,'LineWidth',5)

    er=errorbar([1:length(yvals)],yvals,er_low,er_high);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none'; 
    er.LineWidth=5;
end
