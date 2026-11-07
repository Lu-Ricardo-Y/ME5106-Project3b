%% Verify candidate track boundaries and attach labels from the project PDF
baseDir = fileparts(mfilename('fullpath'));
resultDir = fullfile(baseDir,'me3_step4_results');
if ~isfolder(resultDir), mkdir(resultDir); end
load(fullfile(baseDir,'me3_step2_results','step2_envelope.mat'),'t','filteredRms');
ref = readtable(fullfile(baseDir,'me3_step3_results','candidate_245_tracks.csv'));
dt = median(diff(t));
smoothRms = movmean(filteredRms,3);
cutoffs = [0.00025 0.00035 0.00045 0.00060];
startGrid = zeros(245,numel(cutoffs));
endGrid = zeros(245,numel(cutoffs));
for c=1:numel(cutoffs)
    mask = smoothRms>cutoffs(c);
    starts = find(diff([false;mask])==1);
    ends = find(diff([mask;false])==-1);
    keep=(ends-starts+1)*dt>=0.008;
    starts=starts(keep);ends=ends(keep);
    assert(numel(starts)==245,'Threshold %.6g produced %d tracks',cutoffs(c),numel(starts));
    startGrid(:,c)=t(starts)-dt/2;
    endGrid(:,c)=t(ends)+dt/2;
end
startSpread=max(startGrid,[],2)-min(startGrid,[],2);
endSpread=max(endGrid,[],2)-min(endGrid,[],2);
fprintf('All thresholds: 245 tracks. Median/max start spread: %.3f/%.3f ms.\n',median(startSpread)*1000,max(startSpread)*1000);
fprintf('Median/max end spread: %.3f/%.3f ms.\n',median(endSpread)*1000,max(endSpread)*1000);

%% Explicit PDF page 3 transcription, with unknowns kept as NaN
power = nan(49,1); speed = nan(49,1);
power([1:18,19:21,25:27,31:33])=225;
power([22:24,28:30,34:36,37,45])=175;
power(38)=275; power(39)=325; power(40)=275; power(42)=370;
power(44)=425; power(46)=275; power(47)=325; power(48)=370;
speed([1:3,7:9,13:15,45:49])=850;
speed([4:6,10:12,16:18,37:39])=550;
speed([19:36,40,42,44])=700;
assert(sum(isnan(power))==3 && sum(isnan(speed))==2,'Unexpected missing labels');
assert(isequal(find(isnan(power))',[41 43 49]));
assert(isequal(find(isnan(speed))',[41 43]));
labels=table((1:49)',power,speed,'VariableNames',{'sampleID','power_W','speed_mm_s'});
writetable(labels,fullfile(resultDir,'sample_labels_with_unknowns.csv'));
tracks=join(ref,labels,'Keys','sampleID');
tracks.startThresholdSpread_ms=startSpread*1000;
tracks.endThresholdSpread_ms=endSpread*1000;
writetable(tracks,fullfile(resultDir,'tracks_with_labels_and_boundary_qa.csv'));

%% Six panels covering the start, middle, uncommon and unknown samples
checkIDs=[1 19 37 41 43 49];
fig=figure('Color','w','Name','Selected segmentation quality checks');
for panel=1:numel(checkIDs)
    id=checkIDs(panel);
    rows=tracks(tracks.sampleID==id,:);
    left=max(0,rows.startSeconds(1)-0.05);
    right=rows.endSeconds(5)+0.05;
    subplot(3,2,panel);
    plot(t,filteredRms,'Color',[0.3 0.5 0.8]);hold on;
    xlim([left right]);grid on;
    for k=1:5
        xline(rows.startSeconds(k),'Color',[0.1 0.6 0.2],'LineWidth',0.8);
        xline(rows.endSeconds(k),'Color',[0.8 0.3 0.1],'LineWidth',0.8);
    end
    title(sprintf('Sample %d: P=%g W, V=%g mm/s',id,rows.power_W(1),rows.speed_mm_s(1)));
    xlabel('Time (s)');ylabel('RMS >8 kHz');
end
exportgraphics(fig,fullfile(resultDir,'01_six_sample_boundary_checks.png'),'Resolution',180);
savefig(fig,fullfile(resultDir,'01_six_sample_boundary_checks.fig'));

%% Check whether measured track duration supports speed inference
groupDuration=splitapply(@median,tracks.durationSeconds,tracks.sampleID);
known=~isnan(speed);
fprintf('Known-speed median durations by speed:\n');
for v=[550 700 850]
    d=groupDuration(speed==v);
    fprintf('%d mm/s: n=%d, median=%.4f s, range=%.4f-%.4f s\n',v,numel(d),median(d),min(d),max(d));
end
durationTable=table((1:49)',power,speed,groupDuration, ...
    'VariableNames',{'sampleID','power_W','speed_mm_s','medianTrackDuration_s'});
writetable(durationTable,fullfile(resultDir,'sample_duration_summary.csv'));
fig2=figure('Color','w','Name','Duration versus scan speed');
scatter(speed(known),groupDuration(known)*1000,55,power(known),'filled');grid on;
xlabel('Known scan speed (mm/s)');ylabel('Median detected track duration (ms)');
title('Does acoustic track duration reflect scan speed?');
cb=colorbar;ylabel(cb,'Laser power (W)');
exportgraphics(fig2,fullfile(resultDir,'02_duration_vs_speed.png'),'Resolution',180);
savefig(fig2,fullfile(resultDir,'02_duration_vs_speed.fig'));
fprintf('Step 4 complete: %s\n',resultDir);
