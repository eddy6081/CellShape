function [dR_steps, dTheta_steps] = Compute_dR_dTheta_steps(Persist_Steps, Deflect_Steps)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Compute magnitude of step sizes and theta angle away from
    persistent direction.

    Input Args:
        Persist_Steps: cell array, phenotype x phenotype in shape.
                       Step sizes in the persistent direction taken from
                       initial phenotype to final phenotype as
                       Persist_Steps{initial, final}
        Deflect_Steps: cell array, phenotype x phenotype in shape.
                       Step sizes in the deflected direction taken from
                       initial phenotype to final phenotype as
                       Persist_Steps{initial, final}. Positive is to the
                       right of persistence direction.

    Outputs:
        dR_steps: Cell array, phenotype x phenotype in size.
                  Magnitude of step sizes taken from initial phenotype to
                  final phenotype as dR_steps{initial, final}
        
        dTheta_steps: Cell array, phenotype x phenotype in size.
                      Theta angle from persistence direction of step sizes 
                      taken from initial phenotype to final phenotype as 
                      dR_steps{initial, final}
    %}

    %%
    %calculate dr steps from Persist and deflection, also from these we could
    %calculate a theta (deflection from persistent direction).

    %given variables Persist_Steps, Deflect_Steps

    %calculate dR for each type

    dR_steps = cell(size(Persist_Steps));
    dTheta_steps = cell(size(Persist_Steps));
    %we can't really calculate dtheta. Really, it is dtheta away from persistent
    %direction, which is fine...
    %positive is to the LEFT of positive-persistence, negative is to the right
    %of positive persistence. 

    for p1=1:size(Persist_Steps,1)
        for p2=1:size(Persist_Steps,2)
            dR_steps{p1,p2}=sqrt(Persist_Steps{p1,p2}.^2 + Deflect_Steps{p1,p2}.^2);
            dTheta_steps{p1,p2} = acosfull(Persist_Steps{p1,p2}, Deflect_Steps{p1,p2});
            %atan(Deflect_Steps{p1,p2} ./ Persist_Steps{p1,p2});
            %dTheta_steps{p1,p2}(isnan(dTheta_steps{p1,p2})) = []; %simply delete. %this is the case 
            %where the persistence and deflection are both zero, so there is no
            %theta step. 
            %get rid of any that are nan 
            dR_steps{p1,p2}(isnan(dR_steps{p1,p2}))=[];
            dTheta_steps{p1,p2}(isnan(dTheta_steps{p1,p2}))=[];
        end
    end
end