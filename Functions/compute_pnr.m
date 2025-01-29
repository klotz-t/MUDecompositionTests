%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to compute pulse-to-noise ratio (PNR) for a MU source
% Translated from Python code obtained from openhdemg
% See: https://github.com/GiacomoValliPhD/openhdemg
%
% Input:    ipts = MU source
%           mupulses = time instants of identified MU source peaks
%           constrain_pulses = 1 means the times of firing are considered
%               those in mupulses +/- the number of samples specified in
%               constrain_pulses (e.g., [1 3]). Otherwise, the times of
%               firing are estimated via a heuristic penalty function.
%           ignore_negative_ipts = 1 means normalising by ignoring negative
%               values
%           separate_paired_firings = Whether to treat differently paired
%               and non-paired firings during the estimation of the
%               signal/noise threshold
%
% Output:   pnr = Pulse-to-noise ratio (dB)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function pnr = compute_pnr(ipts, mupulses, fsamp, constrain_pulses, ignore_negative_ipts, separate_paired_firings)

    % Default values
    if nargin < 4
        constrain_pulses = [true, 3];
    end

    if nargin < 6
        separate_paired_firings = false;
    end

    if ignore_negative_ipts
        % Ignore negative values
        ipts = ipts .* abs(ipts);
    end

    % Handle case of no firings
    if isempty(mupulses)
        pnr = NaN;
        return;
    end

    % Extract and normalize the source
    source = ipts(:); % Ensure it's a column vector
    source = source / mean(source(mupulses));

    % Check how to estimate PNR
    if constrain_pulses(1) == true
        % Estimate by mupulses
        start = -constrain_pulses(2);
        stop = constrain_pulses(2);
        extended_mupulses = [];
        for t = start:stop
            extended_mupulses = [extended_mupulses; mupulses + t];
        end
        extended_mupulses = unique(extended_mupulses);

        % Identify noise times
        noise_times = setdiff(mupulses(1):mupulses(end), extended_mupulses);

        % Create clusters
        peak_cluster = source(mupulses);
        noise_cluster = source(noise_times);
        noise_cluster = noise_cluster(~isnan(noise_cluster) & noise_cluster >= 0);

    elseif constrain_pulses(1) == false
        % Estimate by Pi
        idi = diff(mupulses);

        % Remove outlier values (> 500 ms)
        idi = idi(idi <= (fsamp * 0.5));

        % Calculate Pi
        if ~separate_paired_firings
            % Calculate Pi on all IDI
            CoV_all_IDI = std(idi) / mean(idi);
            if isnan(CoV_all_IDI)
                CoV_all_IDI = 0;
            end
            Pi = CoV_all_IDI;
        else
            % Separate paired and non-paired firings
            idinonp = idi(idi >= (fsamp * 0.05));
            idip = idi(idi < (fsamp * 0.05));

            if numel(idinonp) > 1
                CoVIDI = std(idinonp) / mean(idinonp);
            else
                CoVIDI = 0;
            end

            if numel(idip) > 1
                CoVpIDI = std(idip) / mean(idip);
            else
                CoVpIDI = 0;
            end

            % Calculate Pi
            Pi = CoVIDI + CoVpIDI;
        end

        % Create clusters
        peak_cluster = source(source >= Pi);
        noise_cluster = source(source < Pi);
    else
        error('constrain_pulses(1) can only be true or false.');
    end

    % Calculate PNR
    peak_cluster = peak_cluster .^ 2;
    noise_cluster = noise_cluster .^ 2;
    pnr = 10 * log10(mean(peak_cluster) / mean(noise_cluster));
end