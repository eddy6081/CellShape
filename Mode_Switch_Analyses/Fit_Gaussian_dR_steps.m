function [Gauss_params, Gauss_params_error] = Fit_Gaussian_dR_steps(Persist_steps)
    %{
    Author: Chris Eddy
    Date: 05/12/21
    Purpose: Compute parameters of step sizes for each type of phenotype
             transition. Show plots of fits of Gaussian distributions.
    
    Input Args:
        Persist_steps: cell array, shape phenotype x phenotype
                       Contains the magnitude of step sizes taken by cells in
                       each phenotype. 

    Outputs:
        Gauss_params: double, shape phenotype x phenotype x 3 
                      fitted parameters of step sizes in log-normal.

        Gauss_error: double, shape phenotype x phenotype x 6
                     95% CI of fitted parameters in Gauss_params.
    %}

    %% Fit the dR distributions with Gaussian
    if size(Persist_steps,1)==3
        ttls_transition={'AM-AM','AM-ME','AM-IM';...
                         'ME-AM','ME-ME','ME-IM';...
                         'IM-AM','IM-ME','IM-IM'};
    elseif size(Persist_steps,1)==5
        ttls_transition={'AE-AE','AE-FP','AE-BB','AE-LA','AE-IM';...
                         'FP-AE','FP-FP','FP-BB','FP-LA','FP-IM';...
                         'BB-AE','BB-FP','BB-BB','BB-LA','BB-IM';...
                         'LA-AE','LA-FP','LA-BB','LA-LA','LA-IM';...
                         'IM-AE','IM-FP','IM-BB','IM-LA','IM-IM'};
    else
        disp("you must define title types for transitions between phenotypes")
    end
    
    Gauss_params = nan([size(Persist_steps),3]);
    Gauss_params_error = nan([size(Persist_steps),6]);
    count=0;
    figure;
    
    for p1=1:size(Persist_steps,1)
        for p2=1:size(Persist_steps,2)
            count=count+1;
            subplot(size(Persist_steps,1),size(Persist_steps,2),count)
            %create 
            x = Persist_steps{p1,p2};
            if ~isempty(x)
                %histogram. Step sizes are in microns. Bin size 1
                %micrometer.
                %edges, -15,  to 15
                X = histogram(x, 'BinEdges', [-15.5:1:15.5]);
                hold on
                ft = fittype('A*exp(-((x-c)^2)/(2*(s^2)))','independent','x','dependent','y');
                %('1 - 0.5*erfc((x-a)/(b*sqrt(2)))','independent','x','dependent','y');
                xfit = [-15:1:15];
                yfit = X.Values;
                %fo = fitoptions('StartPoint',[max(yfit) 0 1]);
                fitresult = fit(xfit', yfit', ft,'StartPoint',[max(yfit) 0 1]);
                ci = confint(fitresult,0.95);
                
                Gauss_params(p1,p2,:) = [fitresult.A, fitresult.c, fitresult.s];
                Gauss_params_error(p1,p2,:)=reshape(ci,1,[]);
                xfitting = [-20:0.1:20];
                yfitting = fitresult.A * exp(-((xfitting - fitresult.c).^2) / (2*(fitresult.s^2)));
                plot(xfitting,yfitting,'r-','LineWidth',2)
                
            end
            %set axes 
            xlim([-15.5,15.5])
            xticks([-15:5:15])
            xlabel('Step Size (\mum)')
            ylabel('Counts')
            %set('gca','FontSize',12)
            set(gca, 'FontName','Times New Roman','FontSize',12, 'LineWidth',3)
        end
    end
end