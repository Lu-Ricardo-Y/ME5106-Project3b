%% Export candidate boundaries from the validated group-count setting
baseDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(baseDir,'me3_step3_results');
if ~isfolder(resultDir), mkdir(resultDir); end
load(fullfile(baseDir,'me3_step2_results','step2_envelope.mat'),'t','filteredRms');
dt = median(diff(t));
smoothed = movmean(filteredRms,3);
threshold = 0.00035;
active = smoothed>threshold;
starts = find(diff([false;active])==1);
ends = find(diff([active;false])==-1);
keep = (ends-starts+1)*dt>=0.008;
starts = starts(keep);ends = ends(keep);
groupStart = [1;find(t(starts(2:end))-t(ends(1:end-1))>1)+1];
counts = diff([groupStart;numel(starts)+1]);
assert(numel(starts)==245 && numel(groupStart)==49 && all(counts==5), ...
    'Segmentation count mismatch. Review thresholds before exporting.');
sampleID = repelem((1:49)',5);
trackID = repmat((1:5)',49,1);
startSeconds = t(starts)-dt/2;
endSeconds = t(ends)+dt/2;
durationSeconds = endSeconds-startSeconds;
peakRms = zeros(245,1);
for k=1:245
    peakRms(k)=max(filteredRms(starts(k):ends(k)));
end
tracks = table(sampleID,trackID,startSeconds,endSeconds,durationSeconds,peakRms);
writetable(tracks,fullfile(resultDir,'candidate_245_tracks.csv'));
save(fullfile(resultDir,'candidate_245_tracks.mat'),'tracks','threshold','dt');
fprintf('Saved 245 candidate tracks. First and last rows:\n');
disp(tracks([1:5,241:245],:));
fprintf('Duration range: %.3f to %.3f s; median %.3f s\n', ...
    min(durationSeconds),max(durationSeconds),median(durationSeconds));

f=figure('Color','w','Name','49 detected groups');
plot(t,filteredRms,'Color',[0.45 0.65 0.85]); hold on;
for j=1:49
    idx=groupStart(j);
    text(startSeconds(idx),max(filteredRms)*0.94,sprintf('%d',j), ...
        'FontSize',6,'Rotation',90,'HorizontalAlignment','center');
end
xlim([0 t(end)]);grid on;
xlabel('Time (s)');ylabel('RMS >8 kHz');
title('Candidate locations of 49 five-track groups');
exportgraphics(f,fullfile(resultDir,'01_all_49_groups.png'),'Resolution',200);
savefig(f,fullfile(resultDir,'01_all_49_groups.fig'));
