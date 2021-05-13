function [MSD, MSD_error, MSD_AMME, MSD_AMME_error, MSD_All, MSD_All_error]=Calculate_MSD_Data(locsx, locsy, frames, classes)
    %{
    Author: Chris Eddy
    Date: 03/16/21
    Purpose: Calculate mean square displacements for data. Calculates MSD
    for each class in unique(classes), as well as all trajectories. 

    Input Arguments:

        locsx:   cell array, 1 x NumCells. Contains the x positions of cell
                 i in locsx{i}

        locsy:   cell array, 1 x NumCells. Contains the y positions of cell
                 i in locsy{i}

        frames:  cell array, 1 x NumCells. Contains the frame numbers of cell
                 i in frames{i}
        
        classes: cell array, 1 x NumCells. Contains the classifications of
                 cell i in classes{i}

    Outputs:
        
        MSD:        cell array, 1 x phenotype in size. MSD for each
                    phenotype.

        MSD_AMME:   cell array, 1 x 3 in size. MSD for classes restructured
                    as amoeboid and mesenchymal, and intermediate

        MSD_All:    double, 1 x t in size. MSD for all trajectories.

    %}
    
    %% Calculate MSD
    % we need to calculate the MSD for different sets:
    % (1) to find MSD when cells dwell (show if diffusivities are different
    %     for each phenotype)
    % (2) to find MSD for all trajectories of cells, regardless of
    %     phenotypes.
    % we also need (1) for AM, ME states.
    
    %% Determine longest time difference. 
    L = 0;
    for i=1:length(locsx)
        mdiff = frames{i}(end) - frames{i}(1);
        if mdiff > L
            L = mdiff;
        end
    end
    
    %% (1) MSD for phenotype dwelling 
    classes_flat = vertcat(classes{:}); %just determining biggest class number.
    [C,~] = max(unique(classes_flat));
    Squared_Displacements = cell(1,C);

    for i=1:length(locsx)
        frs = frames{i}; %trajectory ith frames
        cs = classes{i}; %trajectory ith classes
        x = locsx{i}; %trajectory ith x position
        y = locsy{i}; %trajectory ith y position
        for t = 1:length(frs)-1
            g = diff(frs(t:end)); %calculate consecutive frame differences
            h = (g==1); %find where frames are continuous
            ind = min([find(h==0,1),length(g)]);%find index of first instance where frames are not continuous, or set to length of frs.
            frs_cont = frs(t:t+ind); %these are the usable frames for this time point.
            cs_cont = cs(t:t+ind); % classes of continuous frames
            x_cont = x(t:t+ind); % x locations of continuous frames
            y_cont = y(t:t+ind); % y locations of continuous frames
            
            %now, for these continuous frames, classes, and locations, we
            %need to determine if the cell dwelled in phenotype. 
            cs_same = (cs_cont==cs_cont(1)); %may look like 1 1 1, or 1 0 1. 
            ind = min([find(cs_same==0,1),length(cs_same)]); %find first instance where not the same. 
            ind = ind-1; %end right before
            frs_cont = frs_cont(1:ind); %these are the usable frames and classes for this time point.
            cs_cont = cs_cont(1:ind); % classes of continuous frames and classes
            x_cont = x_cont(1:ind); % x locations of continuous frames and classes
            y_cont = y_cont(1:ind); % y locations of continuous frames and classes
            
            g = diff(frs_cont); %calculate differences in useable frames and classes
            taus = cumsum(g); %add differences of useable frames and classes to get tau points
            
            for tau_i = 1:length(taus)
                xdiff = x_cont(1+taus(tau_i)) - x_cont(1);
                ydiff = y_cont(1+taus(tau_i)) - y_cont(1);
                if length(Squared_Displacements{cs_cont(1)})<taus(tau_i)
                    Squared_Displacements{cs_cont(1)}{taus(tau_i)} = (xdiff^2) + (ydiff^2);
                else
                    Squared_Displacements{cs_cont(1)}{taus(tau_i)}(end+1)=(xdiff^2) + (ydiff^2);
                end
            end
        end
    end
    
    MSD = cell(1,C);
    MSD_error = cell(1,C);
    for i = 1:C
        MSD{i}=zeros(1,length(Squared_Displacements{i}));
        MSD_error{i}=zeros(1,length(Squared_Displacements{i}));
        for tau = 1:length(Squared_Displacements{i})
            MSD{i}(tau)=mean(Squared_Displacements{i}{tau});
            MSD_error{i}(tau) = std(Squared_Displacements{i}{tau})/sqrt(length(Squared_Displacements{i}{tau}));
        end
    end
    
    %% (1.5) MSD for AM ME states (same code as (1))
    Squared_AMME_Displacements = cell(1,3);

    for i=1:length(locsx)
        frs = frames{i}; %trajectory ith frames
        cs = classes{i}; %trajectory ith classes
        cs(cs==4)=2; %change lamellipodia to filopodia. This is mesenchymal.
        cs(cs==3)=1; %change blebbing to actin-edge. This is amoeboid.
        cs(cs==5)=3; %change intermediate to class 3.
        x = locsx{i}; %trajectory ith x position
        y = locsy{i}; %trajectory ith y position
        for t = 1:length(frs)-1
            g = diff(frs(t:end)); %calculate consecutive frame differences
            h = (g==1); %find where frames are continuous
            ind = min([find(h==0,1),length(g)]);%find index of first instance where frames are not continuous, or set to length of frs.
            frs_cont = frs(t:t+ind); %these are the usable frames for this time point.
            cs_cont = cs(t:t+ind); % classes of continuous frames
            x_cont = x(t:t+ind); % x locations of continuous frames
            y_cont = y(t:t+ind); % y locations of continuous frames
            
            %now, for these continuous frames, classes, and locations, we
            %need to determine if the cell dwelled in phenotype. 
            cs_same = (cs_cont==cs_cont(1)); %may look like 1 1 1, or 1 0 1. 
            ind = min([find(cs_same==0,1),length(cs_same)]); %find first instance where not the same. 
            ind = ind-1; %end right before
            frs_cont = frs_cont(1:ind); %these are the usable frames and classes for this time point.
            cs_cont = cs_cont(1:ind); % classes of continuous frames and classes
            x_cont = x_cont(1:ind); % x locations of continuous frames and classes
            y_cont = y_cont(1:ind); % y locations of continuous frames and classes
            
            g = diff(frs_cont); %calculate differences in useable frames and classes
            taus = cumsum(g); %add differences of useable frames and classes to get tau points
            
            for tau_i = 1:length(taus)
                xdiff = x_cont(1+taus(tau_i)) - x_cont(1);
                ydiff = y_cont(1+taus(tau_i)) - y_cont(1);
                if length(Squared_AMME_Displacements{cs_cont(1)})<taus(tau_i)
                    Squared_AMME_Displacements{cs_cont(1)}{taus(tau_i)} = (xdiff^2) + (ydiff^2);
                else
                    Squared_AMME_Displacements{cs_cont(1)}{taus(tau_i)}(end+1)=(xdiff^2) + (ydiff^2);
                end
            end
        end
    end
    
    MSD_AMME = cell(1,3);
    MSD_AMME_error = cell(1,3);
    for i = 1:3
        MSD_AMME{i}=zeros(1,length(Squared_AMME_Displacements{i}));
        MSD_AMME_error{i}=zeros(1,length(Squared_AMME_Displacements{i}));
        for tau = 1:length(Squared_AMME_Displacements{i})
            MSD_AMME{i}(tau)=mean(Squared_AMME_Displacements{i}{tau});
            MSD_AMME_error{i}(tau) = std(Squared_AMME_Displacements{i}{tau})/sqrt(length(Squared_AMME_Displacements{i}{tau}));
        end
    end
    
    %% (2) MSD for all trajectories
    Sq_diff_all = cell(1,L);
    for i=1:length(locsx)
        frs = frames{i};
        for t=1:length(frs)-1
            g = diff(frs(t:end));
            taus = cumsum(g);
            for tau_i = 1:length(taus)
                xdiff = locsx{i}(t+tau_i) - locsx{i}(t);
                ydiff = locsy{i}(t+tau_i) - locsy{i}(t);
                Sq_diff_all{taus(tau_i)}(end+1) = (xdiff.^2) + (ydiff.^2);
            end
        end
    end
    MSD_All = zeros(1, L);
    MSD_All_error = zeros(1, L);
    for tau=1:L
        MSD_All(tau) = mean(Sq_diff_all{tau});
        MSD_All_error(tau) = std(Sq_diff_all{tau})/sqrt(length(Sq_diff_all{tau}));
    end
end