function [MSD, MSD_error]=Compute_MSD(positions)
    %{
    Input Args:
        positions = double, N_cells x N_dim x N_steps
    %}
    MSD = zeros(size(positions,3),1);
    MSD_error = zeros(size(positions,3),1);
    for tau=1:size(positions,3)-1
        sq_diffs = (positions(:,:,1+tau:end) - positions(:,:,1:end-tau)).^2;
        sum_diffs = sum(sq_diffs,2);
        sum_diffs = squeeze(sum_diffs); %remove dimension
        MSD(tau+1) = mean(sum_diffs(:));
        MSD_error(tau+1) = std(sum_diffs(:))/sqrt(length(sum_diffs(:)));
    end
    
    
end