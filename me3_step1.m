%% ME5106 project 3: read-only inspection of the original recording
% Run this entire script. No additional MATLAB toolboxes are required.
sourceFile = fullfile(fileparts(mfilename('fullpath')), 'data', 'Mini Project 3 Data.mat');
outputDir = fullfile(fileparts(mfilename('fullpath')), 'me3_step1_results');
if ~isfolder(outputDir), mkdir(outputDir); end
assert(isfile(sourceFile), 'Cannot find input file: %s', sourceFile);
info = whos('-file', sourceFile);
disp(info);
loaded = load(sourceFile, 'IEPE0_Raw');
assert(isfield(loaded, 'IEPE0_Raw'), 'IEPE0_Raw is missing.');
x = double(loaded.IEPE0_Raw(:));
clear loaded;
fs = 160000; % Stated in the project description
nSamples = numel(x);
durationSeconds = nSamples/fs;
nInvalid = sum(~isfinite(x));
fprintf('\nSamples: %d\nSampling rate: %d Hz\nDuration: %.3f s\nNonfinite samples: %d\n', ...
    nSamples, fs, durationSeconds, nInvalid);
assert(nInvalid == 0, 'Nonfinite samples found. Inspect before processing.');
dcOffset = mean(x);
x = x-dcOffset;
fprintf('Removed DC offset: %.8g\n', dcOffset);

%% One-millisecond RMS and min/max envelopes, using every sample
blockSize = round(0.001*fs);
nBlocks = floor(nSamples/blockSize);
rmsEnvelope = zeros(nBlocks,1);
minEnvelope = zeros(nBlocks,1);
maxEnvelope = zeros(nBlocks,1);
batchBlocks = 10000;
for firstBlock = 1:batchBlocks:nBlocks
    lastBlock = min(firstBlock+batchBlocks-1,nBlocks);
    firstSample = (firstBlock-1)*blockSize+1;
    lastSample = lastBlock*blockSize;
    chunk = reshape(x(firstSample:lastSample),blockSize,[]);
    rmsEnvelope(firstBlock:lastBlock) = sqrt(mean(chunk.^2,1)).';
    minEnvelope(firstBlock:lastBlock) = min(chunk,[],1).';
    maxEnvelope(firstBlock:lastBlock) = max(chunk,[],1).';
end
clear chunk;
tEnvelope = ((0:nBlocks-1)'+0.5)*blockSize/fs;

%% Global map: min/max retains spikes that display decimation could miss
overviewFigure = figure('Name','ME5106 Step 1 - Full recording','Color','w');
subplot(2,1,1);
plot(tEnvelope,maxEnvelope,'Color',[0.2 0.45 0.8]); hold on;
plot(tEnvelope,minEnvelope,'Color',[0.2 0.45 0.8]);
grid on; xlim([0 durationSeconds]);
xlabel('Time (s)'); ylabel('Amplitude (recorded units)');
title('Full recording: minimum and maximum in each 1 ms block');
subplot(2,1,2);
plot(tEnvelope,rmsEnvelope); grid on; xlim([0 durationSeconds]);
xlabel('Time (s)'); ylabel('RMS (recorded units)');
title('Full recording: 1 ms RMS, before filtering');
drawnow;
savefig(overviewFigure,fullfile(outputDir,'01_full_recording.fig'));
exportgraphics(overviewFigure,fullfile(outputDir,'01_full_recording.png'),'Resolution',160);

%% Inspect the beginning without assuming where the first valid track starts
localRange = [0 min(15,durationSeconds)];
localMask = tEnvelope>=localRange(1) & tEnvelope<=localRange(2);
localFigure = figure('Name','ME5106 Step 1 - First 15 seconds','Color','w');
subplot(2,1,1);
plot(tEnvelope(localMask),maxEnvelope(localMask)); hold on;
plot(tEnvelope(localMask),minEnvelope(localMask)); grid on;
xlim(localRange); xlabel('Time (s)'); ylabel('Amplitude (recorded units)');
title('First 15 seconds: 1 ms min/max envelope');
subplot(2,1,2);
plot(tEnvelope(localMask),rmsEnvelope(localMask)); grid on;
xlim(localRange); xlabel('Time (s)'); ylabel('RMS (recorded units)');
title('First 15 seconds: 1 ms RMS, before filtering');
drawnow;
savefig(localFigure,fullfile(outputDir,'02_first_15_seconds.fig'));
exportgraphics(localFigure,fullfile(outputDir,'02_first_15_seconds.png'),'Resolution',160);
save(fullfile(outputDir,'step1_envelopes.mat'), 'sourceFile','fs','nSamples', ...
    'durationSeconds','dcOffset','blockSize','tEnvelope','rmsEnvelope', ...
    'minEnvelope','maxEnvelope');
fprintf('\nSTEP 1 COMPLETE. Figures and envelopes saved to:\n%s\n',outputDir);
fprintf('These are unfiltered inspection plots, not final track boundaries.\n');
