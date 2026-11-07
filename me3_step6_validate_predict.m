%% Group-level validation and unknown parameter prediction
base = fileparts(mfilename('fullpath'));
data = load(fullfile(base,'me3_step5_results','acoustic_features.mat'),'groupFeatures');
g = data.groupFeatures;
known = isfinite(g.power_W) & isfinite(g.speed_mm_s);
assert(nnz(known)==46 && isequal(g.sampleID(~known)',[41 43 49]));
% Use one median feature vector per group. Band fractions capture spectral
% shape, while the log RMS captures relative acoustic amplitude.
X = [g.duration_s, log(g.highRms), g.band8_20, g.band20_40, g.centroidHz];
assert(all(isfinite(X),'all'));
Y = [g.power_W(known),g.speed_mm_s(known)];
trainX = X(known,:);
[~,~,combo] = unique(Y,'rows');
count = accumarray(combo,1);
n = size(Y,1);
lambda = 5;
% Fixed ridge penalty. Centering and scaling are estimated within each fold.
oof = nan(n,2,3,2);
for scheme = 1:2
    for k = 1:n
        fit = (1:n)'~=k;
        if scheme==2
            fit = combo~=combo(k);
        end
        oof(k,:,1,scheme) = median(Y(fit,:),1);
        oof(k,:,2,scheme) = ridgePredict(trainX(fit,1),Y(fit,:),trainX(k,1),lambda);
        oof(k,:,3,scheme) = ridgePredict(trainX(fit,:),Y(fit,:),trainX(k,:),lambda);
    end
end
unknownID = g.sampleID(~known);
prediction = ridgePredict(trainX,Y,X(~known,:),lambda);
prediction(isfinite(g.speed_mm_s(~known)),2) = g.speed_mm_s(g.sampleID==49);
outOfRange = X(~known,:) < min(trainX,[],1) | X(~known,:) > max(trainX,[],1);

dest = fullfile(base,'me3_step6_results');
if ~isfolder(dest), mkdir(dest); end
save(fullfile(dest,'validation_and_predictions.mat'),'oof','Y','combo','count', ...
    'unknownID','prediction','outOfRange','lambda','g','X','known');
fid = fopen(fullfile(dest,'validation_and_predictions.txt'),'w');
assert(fid>=0);
cleanup = onCleanup(@() fclose(fid));
schemes = {'Leave one sample group out','Leave one power-speed combination out'};
models = {'Median baseline','Duration-only ridge','Acoustic-feature ridge'};
sets = {true(n,1),count(combo)>1,count(combo)==1};
names = {'All','Repeated combinations','Singleton combinations'};
fprintf(fid,'Validation: 46 fully labelled groups, %d distinct settings\n',numel(count));
for scheme = 1:2
    fprintf(fid,'\n%s\n',schemes{scheme});
    fprintf(fid,'Subset | Model | N | Power MAE W | Power RMSE W | Speed MAE mm/s | Speed RMSE mm/s\n');
    for s = 1:3
        idx = sets{s};
        for m = 1:3
            e = oof(idx,:,m,scheme)-Y(idx,:);
            mae = mean(abs(e),1);
            rmse = sqrt(mean(e.^2,1));
            fprintf(fid,'%s | %s | %d | %.2f | %.2f | %.2f | %.2f\n', ...
                names{s},models{m},nnz(idx),mae(1),rmse(1),mae(2),rmse(2));
        end
    end
end
fprintf(fid,'\nPredictions from the acoustic-feature ridge (lambda=%.1f)\n',lambda);
fprintf(fid,'Group | Power W | Speed mm/s | Speed source | Features outside training range\n');
for k=1:numel(unknownID)
    source = 'predicted';
    if unknownID(k)==49, source = 'given'; end
    fprintf(fid,'%d | %.2f | %.2f | %s | %d / %d\n',unknownID(k), ...
        prediction(k,1),prediction(k,2),source,nnz(outOfRange(k,:)),size(X,2));
end
clear cleanup;
fprintf('%s',fileread(fullfile(dest,'validation_and_predictions.txt')));

function prediction = ridgePredict(x,y,query,lambda)
    center = mean(x,1);
    scale = std(x,0,1);
    scale(scale==0) = 1;
    z = (x-center)./scale;
    meanY = mean(y,1);
    coefficients = (z'*z + lambda*eye(size(z,2))) \ (z'*(y-meanY));
    prediction = (query-center)./scale*coefficients + meanY;
end
