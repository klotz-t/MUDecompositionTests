function [sil,pnr] = compute_sil(ipts, mupulses, ignore_negative_ipts)
    % Calculate the Silhouette score (SIL)
    %
    % Parameters
    % ----------
    % ipts : numeric array
    %     The decomposed source (or pulse train) of the MU of interest.
    % mupulses : numeric array
    %     The time of firing of the MU of interest.
    % ignore_negative_ipts : logical
    %     If true, only use positive ipts during peak and noise clustering.
    %
    % Returns
    % -------
    % sil : float
    %     The SIL score.
    % pnr : float
    %     The PNR score.

    if nargin < 3
        ignore_negative_ipts = true;
    end

    % Handle case of no firings
    if isempty(mupulses)
        sil = NaN;
        return;
    end

    % Extract source
    source = ipts(:); % Ensure it's a column vector

    if ignore_negative_ipts
        % Ignore negative values
        source = source .* abs(source);
    end

    % Create clusters
    peak_cluster = source(mupulses);
    noise_cluster = source;
    noise_cluster(mupulses) = []; % Remove peaks from noise cluster

    % Create centroids for each cluster
    peak_centroid = mean(peak_cluster);
    noise_centroid = mean(noise_cluster);

    % Calculate within-cluster sums of point-to-centroid distances
    intra_sums = sum((peak_cluster - peak_centroid) .^ 2);

    % Calculate between-cluster sums of point-to-centroid distances
    inter_sums = sum((peak_cluster - noise_centroid) .^ 2);

    % Calculate sil and pnr
    sil = (inter_sums - intra_sums) / max(intra_sums, inter_sums);
end