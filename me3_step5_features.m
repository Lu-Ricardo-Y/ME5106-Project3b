%% Physical, compact acoustic features; one row per printed track and sample
baseDir=fileparts(mfilename('fullpath'));
resultDir=fullfile(baseDir,'me3_step5_results');
if ~isfolder(resultDir), mkdir(resultDir); end
sourceFile=fullfile(fileparts(mfilename('fullpath')), 'data', 'Mini Project 3 Data.mat');
tracks=readtable(fullfile(baseDir,'me3_step4_results','tracks_with_labels_and_boundary_qa.csv'));
fs=160000;
data=load(sourceFile,'IEPE0_Raw');
x=double(data.IEPE0_Raw(:));clear data;
n=height(tracks);
rawRms=zeros(n,1); highRms=zeros(n,1);
band8_20=zeros(n,1);band20_40=zeros(n,1);band40_70=zeros(n,1);
centroidHz=zeros(n,1);
for k=1:n
    a=max(1,floor(tracks.startSeconds(k)*fs)+1);
    b=min(numel(x),ceil(tracks.endSeconds(k)*fs));
    y=x(a:b);y=y-mean(y);
    rawRms(k)=sqrt(mean(y.^2));
    N=numel(y);
    w=0.5-0.5*cos(2*pi*(0:N-1)'/(N-1));
    F=fft(y.*w);
    power=abs(F(1:floor(N/2)+1)).^2;
    freq=(0:numel(power)-1)'*fs/N;
    energyHigh=sum(power(freq>=8000));
    % Windowed high-band RMS is a comparable relative energy measure.
    highRms(k)=sqrt(2*energyHigh/max(N*sum(w.^2),eps));
    if energyHigh>0
        band8_20(k)=sum(power(freq>=8000 & freq<20000))/energyHigh;
        band20_40(k)=sum(power(freq>=20000 & freq<40000))/energyHigh;
        band40_70(k)=sum(power(freq>=40000 & freq<70000))/energyHigh;
        centroidHz(k)=sum(freq(freq>=8000).*power(freq>=8000))/energyHigh;
    end
end
features=tracks(:,{'sampleID','trackID','power_W','speed_mm_s','durationSeconds'});
features.rawRms=rawRms;features.highRms=highRms;
features.band8_20=band8_20;features.band20_40=band20_40;
features.band40_70=band40_70;features.centroidHz=centroidHz;
sampleID=(1:49)';
group=@(v) splitapply(@median,v,features.sampleID);
groupFeatures=table(sampleID,group(features.power_W),group(features.speed_mm_s), ...
    group(features.durationSeconds),group(rawRms),group(highRms), ...
    group(band8_20),group(band20_40),group(band40_70),group(centroidHz), ...
    'VariableNames',{'sampleID','power_W','speed_mm_s','duration_s', ...
    'rawRms','highRms','band8_20','band20_40','band40_70','centroidHz'});
% median() of all-NaN groups remains NaN; unknown labels stay unknown.
save(fullfile(resultDir,'acoustic_features.mat'),'features','groupFeatures');

%% Controlled comparisons: same speed, changing power; same power, changing speed
fig=figure('Color','w','Name','Physical feature checks');
subplot(2,2,1);
sel=groupFeatures.speed_mm_s==700 & ~isnan(groupFeatures.power_W);
scatter(groupFeatures.power_W(sel),groupFeatures.highRms(sel),40,'filled');grid on;
xlabel('Power (W), speed = 700 mm/s');ylabel('High-band energy measure');
subplot(2,2,2);
sel=groupFeatures.speed_mm_s==850 & ~isnan(groupFeatures.power_W);
scatter(groupFeatures.power_W(sel),groupFeatures.highRms(sel),40,'filled');grid on;
xlabel('Power (W), speed = 850 mm/s');ylabel('High-band energy measure');
subplot(2,2,3);
sel=groupFeatures.power_W==225 & ~isnan(groupFeatures.speed_mm_s);
scatter(groupFeatures.speed_mm_s(sel),groupFeatures.duration_s(sel)*1000,40,'filled');grid on;
xlabel('Speed (mm/s), power = 225 W');ylabel('Detected duration (ms)');
subplot(2,2,4);
scatter(groupFeatures.speed_mm_s(sel),groupFeatures.centroidHz(sel)/1000,40,'filled');grid on;
xlabel('Speed (mm/s), power = 225 W');ylabel('Spectral centroid above 8 kHz (kHz)');
exportgraphics(fig,fullfile(resultDir,'01_controlled_feature_checks.png'),'Resolution',170);
savefig(fig,fullfile(resultDir,'01_controlled_feature_checks.fig'));
fprintf('Feature extraction complete: %d tracks and %d groups.\n',height(features),height(groupFeatures));
fprintf('Known labels: %d complete; unknown power groups: ',sum(~isnan(groupFeatures.power_W)&~isnan(groupFeatures.speed_mm_s)));
fprintf('%d ',groupFeatures.sampleID(isnan(groupFeatures.power_W)));fprintf('\n');
