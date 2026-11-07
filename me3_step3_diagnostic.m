%% Diagnostic only: identify candidate tracks from the high-frequency envelope
baseDir = fileparts(mfilename('fullpath'));
load(fullfile(baseDir,'me3_step2_results','step2_envelope.mat'),'t','filteredRms');
dt = median(diff(t));
smoothed = movmean(filteredRms,3);
thresholds = [0.00025 0.00035 0.00045 0.00060];
fprintf('Envelope baseline median %.6g; 99th percentile %.6g\n',median(smoothed),prctile(smoothed,99));
for cutoff=thresholds
    active = smoothed>cutoff;
    starts = find(diff([false;active])==1);
    ends = find(diff([active;false])==-1);
    keep = (ends-starts+1)*dt>=0.008;
    starts = starts(keep);ends = ends(keep);
    % Tracks within a set are close in time; an inter-set gap is much longer.
    groupStart = [1;find(t(starts(2:end))-t(ends(1:end-1))>1)+1];
    counts = diff([groupStart;numel(starts)+1]);
    fprintf('Threshold %.6g: %d track candidates, %d sets; group counts: ',cutoff,numel(starts),numel(groupStart));
    fprintf('%d ',counts);fprintf('\n');
end
