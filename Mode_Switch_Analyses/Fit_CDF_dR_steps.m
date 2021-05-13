function [dR_params_useful, dR_params, dR_params_error] = Fit_CDF_dR_steps(dR_steps)
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
        dR_params_useful: double, shape phenotype x phenotype x 2 
                          conversion of parameters from log-normal to
                          pixel/micron scale.

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
    
    dR_params_useful = cell(size(dR_steps));
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
                %remember, CDF is [number of elements in sample <= t] / N
                %where N is the total number of samples.
                ordered = sort(x(~isnan(x)),'ascend');
                cdfy = zeros(1,length(ordered));
                for item = 1:length(ordered)
                    cdfy(item) = sum(ordered<=ordered(item)) / length(ordered);
                end
                %[0:length(x(~isnan(x)))-1];
                %cdfy = cumsum(cdfy)/sum(cdfy);
                %fit this. 
                ft = fittype('1 - 0.5*erfc((log(x)-a)/(sqrt(2)*b))','independent','x','dependent','y');
                %('1 - 0.5*erfc((x-a)/(b*sqrt(2)))','independent','x','dependent','y');

                ordered=unique(ordered);
                cdfy=unique(cdfy);

                fitresult = fit(ordered', cdfy', ft);
                ci = confint(fitresult,0.95);

                scatter(ordered,cdfy,60,'filled')
                hold on
                xlist = linspace(min(x),max(x));
                ylist = 1 - 0.5*erfc((log(xlist)-fitresult.a)/(sqrt(2)*fitresult.b));
                plot(xlist,ylist,'r-','LineWidth',4)
                title(ttls_transition{p1,p2})
                ylabel('Cum. Prob')
                xlabel('Step Size (\mum)')

                lnmu=fitresult.a;
                lnsigma=fitresult.b;

                finalmu = exp(lnmu)*exp(0.5*lnsigma^2);
                finalvar = exp((2*lnmu)+(lnsigma^2))*(exp(lnsigma^2)-1);

                dR_params_useful{p1,p2}=[finalmu,finalvar];
                dR_params(p1,p2,:)=[lnmu,lnsigma];
                dR_params_error(p1,p2,:)=reshape(ci,1,[]);

            else
                dR_params_useful{p1,p2}=[NaN,NaN];
                dR_params(p1,p2,:)=[NaN,NaN];
                dR_params_error(p1,p2,:)=[NaN,NaN,NaN,NaN];
                title(ttls_transition{p1,p2})
                ylabel('Cum. Prob')
                xlabel('Step Size (\mum)')
            end
        end
    end
    
    figure;

    ttls = cell2mat(reshape(ttls_transition',[],1));
    dR_flat = cell2mat(reshape(dR_params_useful,[],1));

end