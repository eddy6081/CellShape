function [dR_params, dR_params_error] = Fit_LogNorm_dR_steps(dR_steps)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Compute parameters of step sizes for each type of phenotype
             transition. Show plots of fits of log-normal distributions.
    
    Input Args:
        dR_steps: cell array, shape phenotype x phenotype
                  Contains the magnitude of step sizes taken by cells in
                  each phenotype. 

    Outputs:
        dR_params: double, shape phenotype x phenotype x 2 
                   fitted parameters of step sizes in log-normal.

        dR_params_error: double, shape phenotype x phenotype x 4
                         95% CI of fitted parameters in dR_params.
    %}
    
    %% Fit the dR distributions as log-normal cdf
    if size(dR_steps,1)==3
        ttls_transition={'AM-AM','AM-ME','AM-IM';...
                         'ME-AM','ME-ME','ME-IM';...
                         'IM-AM','IM-ME','IM-IM'};
    elseif size(dR_steps,1)==5
        ttls_transition={'AE-AE','AE-FP','AE-BB','AE-LA','AE-IM';...
                         'FP-AE','FP-FP','FP-BB','FP-LA','FP-IM';...
                         'BB-AE','BB-FP','BB-BB','BB-LA','BB-IM';...
                         'LA-AE','LA-FP','LA-BB','LA-LA','LA-IM';...
                         'IM-AE','IM-FP','IM-BB','IM-LA','IM-IM'};
    else
        disp("you must define title types for transitions between phenotypes")
    end
    
    dR_params = nan([size(dR_steps),2]);
    dR_params_error = nan([size(dR_steps),4]);
    count=0;
    figure;
    for p1=1:size(dR_steps,1)
        for p2=1:size(dR_steps,2)
            count=count+1;
            subplot(size(dR_steps,1),size(dR_steps,2),count)
            %create 
            x = dR_steps{p1,p2};
            if ~isempty(x)
                X = histogram(x,'BinEdges', [-0.5:1:20.5]);
                hold on
                %ft = fittype('A*exp(m + s*x)', 'independent','x','dependent','y');
                %ft = fittype('1 - 0.5*erfc((log(x)-a)/(sqrt(2)*b))','independent','x','dependent','y');
                %('1 - 0.5*erfc((x-a)/(b*sqrt(2)))','independent','x','dependent','y');
                xfit = [0:1:20];
                yfit = X.Values;
                [pHat, ci] = lognfit(x(x>0));

                dR_params(p1,p2,:) = [pHat(1), pHat(2)];
                dR_params_error(p1,p2,:)=reshape(ci,1,[]);
                xfitting = [0:0.1:21];
                yfitting = lognpdf(xfitting,pHat(1),pHat(2));%exp(pHat(1) + pHat(2)*xfitting);
                yfitting = yfitting * (max(yfit)/max(yfitting));
                plot(xfitting,yfitting,'r-','LineWidth',2)
                title(ttls_transition{p1,p2})
                ylabel('Counts')
                xlabel('Step Size (\mum)')

            else
                dR_params(p1,p2,:)=[NaN,NaN];
                dR_params_error(p1,p2,:)=[NaN,NaN,NaN,NaN];
                title(ttls_transition{p1,p2})
                ylabel('Cum. Prob')
                xlabel('Step Size (\mum)')
            end
            %set axes 
            xlim([0,20])
            xticks([0:5:20])
            set(gca, 'FontName','Times New Roman','FontSize',12, 'LineWidth',3)
        end
    end
end