%% ME5106 project 3, step 2: exploratory spectrum and high-frequency envelope
sourceFile = fullfile(fileparts(mfilename('fullpath')), 'data', 'Mini Project 3 Data.mat');
outputDir = fullfile(fileparts(mfilename('fullpath')), 'me3_step2_results');
if ~isfolder(outputDir), mkdir(outputDir); end
if ~exist('x','var') || numel(x) ~= 41120000
    loaded = load(sourceFile,'IEPE0_Raw');
    x = double(loaded.IEPE0_Raw(:));
    clear loaded;
end
x = x(:);
fs = 160000;
N = numel(x);
durationSeconds = N/fs;
fprintf('Step 2: %d samples, %.3f seconds.\n',N,durationSeconds);

%% Compare two explicitly chosen regions; their physical interpretation is provisional
quietRange = [0.5 3.0];
activeRange = [6.0 9.0];
quiet = x(floor(quietRange(1)*fs)+1:floor(quietRange(2)*fs));
active = x(floor(activeRange(1)*fs)+1:floor(activeRange(2)*fs));
[f,pQuiet] = localSpectrum(quiet,fs);
[~,pActive] = localSpectrum(active,fs);
spectrumFigure = figure('Name','Step 2 - Spectrum comparison','Color','w');
plot(f/1000,10*log10(pQuiet+eps),'LineWidth',1); hold on;
plot(f/1000,10*log10(pActive+eps),'LineWidth',1);
xlabel('Frequency (kHz)'); ylabel('Power per frequency bin (dB, relative)');
title('Candidate quiet (0.5-3 s) vs active (6-9 s) spectrum');
legend('Candidate quiet','Candidate active','Location','best');
xlim([0 80]); grid on;
exportgraphics(spectrumFigure,fullfile(outputDir,'01_spectrum_comparison.png'),'Resolution',160);
savefig(spectrumFigure,fullfile(outputDir,'01_spectrum_comparison.fig'));

%% Exploratory 8 kHz high-pass in one-second FFT blocks
% 8 kHz follows the handout's lowest-10%-of-0-80-kHz example.
% This sharp FFT cutoff is for locating activity, not final feature extraction.
cutoffHz = 8000;
blockSamples = fs;
rmsBlockSamples = 160; % 1 ms
nSeconds = ceil(N/blockSamples);
nRmsBlocks = floor(N/rmsBlockSamples);
filteredRms = zeros(nRmsBlocks,1);
for sec = 1:nSeconds
    first = (sec-1)*blockSamples+1;
    last = min(sec*blockSamples,N);
    chunk = x(first:last);
    chunk = chunk-mean(chunk);
    n = numel(chunk);
    bins = (0:n-1)'*fs/n;
    keep = bins>=cutoffHz & bins<=fs-cutoffHz;
    spectrum = fft(chunk);
    spectrum(~keep) = 0;
    filtered = real(ifft(spectrum));
    blockFirst = floor((first-1)/rmsBlockSamples)+1;
    blockCount = floor(n/rmsBlockSamples);
    if blockCount>0
        y = reshape(filtered(1:blockCount*rmsBlockSamples),rmsBlockSamples,[]);
        filteredRms(blockFirst:blockFirst+blockCount-1) = sqrt(mean(y.^2,1)).';
    end
end
t = ((0:nRmsBlocks-1)'+0.5)*rmsBlockSamples/fs;
overviewFigure = figure('Name','Step 2 - High-frequency activity','Color','w');
plot(t,filteredRms); grid on;
xlabel('Time (s)'); ylabel('RMS above 8 kHz (recorded units)');
title('Full recording: candidate high-frequency activity');
xlim([0 durationSeconds]);
exportgraphics(overviewFigure,fullfile(outputDir,'02_filtered_full.png'),'Resolution',160);
savefig(overviewFigure,fullfile(outputDir,'02_filtered_full.fig'));
localFigure = figure('Name','Step 2 - First 15 seconds filtered','Color','w');
plot(t,filteredRms); grid on;
xlabel('Time (s)'); ylabel('RMS above 8 kHz (recorded units)');
title('First 15 seconds: candidate high-frequency activity');
xlim([0 15]);
exportgraphics(localFigure,fullfile(outputDir,'03_filtered_first_15s.png'),'Resolution',160);
savefig(localFigure,fullfile(outputDir,'03_filtered_first_15s.fig'));
save(fullfile(outputDir,'step2_envelope.mat'),'fs','cutoffHz','t','filteredRms', ...
    'quietRange','activeRange','f','pQuiet','pActive');
fprintf('Step 2 complete. Results: %s\n',outputDir);

function [frequency,power] = localSpectrum(signal,fs)
    signal = signal(:);
    nfft = 4096;
    hop = nfft/2;
    window = 0.5-0.5*cos(2*pi*(0:nfft-1)'/(nfft-1));
    power = zeros(nfft/2+1,1);
    count = 0;
    for start = 1:hop:(numel(signal)-nfft+1)
        frame = signal(start:start+nfft-1);
        frame = (frame-mean(frame)).*window;
        spectrum = fft(frame);
        power = power+abs(spectrum(1:nfft/2+1)).^2;
        count = count+1;
    end
    power = power/max(count,1)/sum(window.^2);
    frequency = (0:nfft/2)'*fs/nfft;
end
