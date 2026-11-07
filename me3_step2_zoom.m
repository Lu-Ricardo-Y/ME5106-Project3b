resultDir = fullfile(fileparts(mfilename('fullpath')),'me3_step2_results');
load(fullfile(resultDir,'step2_envelope.mat'),'t','filteredRms');
centers = [6.4 12.1 18.0];
f = figure('Color','w','Name','Step 2 - Event zoom');
for j=1:3
    subplot(3,1,j);
    plot(t,filteredRms); grid on;
    xlim(centers(j)+[-0.3 0.3]);
    xlabel('Time (s)'); ylabel('RMS >8 kHz');
    title(sprintf('Candidate event near %.1f s',centers(j)));
end
exportgraphics(f,fullfile(resultDir,'04_event_zoom.png'),'Resolution',180);
savefig(f,fullfile(resultDir,'04_event_zoom.fig'));
