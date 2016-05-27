function LSCAN_concatenateResults_CAAX( dirname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% then get all the subfolders in this folder (that have LACT, in this case)
listing = dir(strcat(dirname,'/*CAAX*'));
subfoldersi = [listing.isdir];
subfolderNames = {listing(subfoldersi).name};
% get all the .mat files in each of the subfolders
for f = 1:length(subfolderNames)
    mat_files{f} = dir(strcat(dirname,'/',subfolderNames{f},'/*.mat'));
    % concatenate data from all the mat files in a subfolder
    if ~isempty(mat_files{f})
        shiftsAndPeriod{f} = nan(length(mat_files{f})*2,3);
        maxRates{f} = nan(length(mat_files{f})*2,3);
        maxChanges{f} = nan(length(mat_files{f})*2,3);
        for n=1:2:length(mat_files{f})*2
            mat_files{f}((n+1)/2).name
            data = load(strcat(dirname,'/',subfolderNames{f},'/',mat_files{f}((n+1)/2).name));
            % load linescan data:
            % pole 1                
            shiftsAndPeriod{f}(n,:)=[data.Parameters.Linescan.ImaxShiftRelativeToFullContraction(1,1) data.Parameters.Linescan.WidthShiftRelativeToFullContraction(1,1) data.Parameters.Linescan.period(1,1)/10];
            maxRates{f}(n,:)=[data.Parameters.Linescan.maximumRufflingRate(1,1) data.Parameters.Linescan.maximumFoldingRate(1,1) -data.Parameters.Linescan.maximumContractionRate(1,1)];
            maxChanges{f}(n,:) = [data.Parameters.Linescan.maximumRuffling(1,1) data.Parameters.Linescan.maximumWidthChange(1,1) data.Parameters.Linescan.maximumSurfaceChange(1,1)];
            % pole 2
            shiftsAndPeriod{f}(n+1,:)=[data.Parameters.Linescan.ImaxShiftRelativeToFullContraction(2,1) data.Parameters.Linescan.WidthShiftRelativeToFullContraction(2,1) data.Parameters.Linescan.period(2,1)/10];
            maxRates{f}(n+1,:)=[data.Parameters.Linescan.maximumRufflingRate(2,1) data.Parameters.Linescan.maximumFoldingRate(2,1) -data.Parameters.Linescan.maximumContractionRate(2,1)];
            maxChanges{f}(n+1,:) = [data.Parameters.Linescan.maximumRuffling(2,1) data.Parameters.Linescan.maximumWidthChange(2,1) data.Parameters.Linescan.maximumSurfaceChange(2,1)];
            clear data
        end
    else
        shiftsAndPeriod{f}(1,:)=[NaN NaN NaN];
        maxRates{f}(1,:)=[NaN NaN NaN];
        maxChanges{f}(1,:) = [NaN NaN NaN];
    end
end

% make all the same length
for f=1:length(subfolderNames)
    shiftsAndPeriod{f}(end+1:20,:)=NaN;
    maxRates{f}(end+1:20,:)=NaN;
    maxChanges{f}(end+1:20,:)=NaN;
end
shiftsAndPeriodS = nan(20,length(subfolderNames)*3);
maxRatesS = nan(20,length(subfolderNames)*3);
maxChangesS = nan(20,length(subfolderNames)*3);
for f=1:3:length(subfolderNames)*3
    % different genotypes in different columns
    shiftsAndPeriodS(:,f:f+2)=shiftsAndPeriod{(f+2)/3};
    maxRatesS(:,f:f+2)=maxRates{(f+2)/3};
    maxChangesS(:,f:f+2) = maxChanges{(f+2)/3};
end

h=figure('Position',[20 20 1000 700]);
file_title='results_lscan';
style={'b','r','g','k','m','y','y','c','b','r','g','k','m'};
for n=1:size(shiftsAndPeriodS,2)
    labels{n} = subfolderNames{floor((n-1)/3)+1};
end


for n=1:3:size(shiftsAndPeriodS,2)
    labels1{n} = strcat(labels{n},'_INTENSITY SHIFT');
    labels1{n+1} = strcat(labels{n+1},'_WIDTH SHIFT');
    labels1{n+2} = strcat(labels{n+2},'_PERIOD/10');
end
subplot(1,3,1)
boxplot(shiftsAndPeriodS,'labels',labels1,'labelorientation','inline');
for n=1:3:size(shiftsAndPeriodS,2)
    hold on;
    plot((n-1)+0.9+0.2*rand(1,length(shiftsAndPeriodS(:,n))),shiftsAndPeriodS(:,n),strcat(style{(n+2)/3},'.'),...
        (n-1)+1+0.9+0.2*rand(1,length(shiftsAndPeriodS(:,n+1))),shiftsAndPeriodS(:,n+1),strcat(style{(n+2)/3},'.'),...
        (n-1)+2+0.9+0.2*rand(1,length(shiftsAndPeriodS(:,n+2))),shiftsAndPeriodS(:,n+2),strcat(style{(n+2)/3},'.'))
    axis([0 length(subfolderNames)*3+1 -60 60])
    title('shifts','FontSize',18);
end

for n=1:3:size(maxRatesS,2)
    labels2{n} = strcat(labels{n},'_RUFFLING');
    labels2{n+1} = strcat(labels{n+1},'_FOLDING');
    labels2{n+2} = strcat(labels{n+2},'_CONTRACTION');
end
subplot(1,3,2)
boxplot(maxRatesS,'labels',labels2,'labelorientation','inline');
for n=1:3:size(maxRatesS,2)
    hold on;
    plot((n-1)+0.9+0.2*rand(1,length(maxRatesS(:,n))),maxRatesS(:,n),strcat(style{(n+2)/3},'.'),...
        (n-1)+1+0.9+0.2*rand(1,length(maxRatesS(:,n+1))),maxRatesS(:,n+1),strcat(style{(n+2)/3},'.'),...
        (n-1)+2+0.9+0.2*rand(1,length(maxRatesS(:,n+2))),maxRatesS(:,n+2),strcat(style{(n+2)/3},'.'))
    axis([0 length(subfolderNames)*3+1 0 0.015])
    title('maximum rates','FontSize',18);
end


for n=1:3:size(maxChangesS,2)
    labels3{n} = strcat(labels{n},'_INTENSITY');
    labels3{n+1} = strcat(labels{n+1},'_WIDTH');
    labels3{n+2} = strcat(labels{n+2},'_SURFACE AREA');
end
subplot(1,3,3)
boxplot(maxChangesS,'labels',labels3,'labelorientation','inline');
for n=1:3:size(maxChangesS,2)
    hold on;
    plot((n-1)+0.9+0.2*rand(1,length(maxChangesS(:,n))),maxChangesS(:,n),strcat(style{(n+2)/3},'.'),...
        (n-1)+1+0.9+0.2*rand(1,length(maxChangesS(:,n))),maxChangesS(:,n+1),strcat(style{(n+2)/3},'.'),...
        (n-1)+2+0.9+0.2*rand(1,length(maxChangesS(:,n))),maxChangesS(:,n+2),strcat(style{(n+2)/3},'.'))
    axis([0 length(subfolderNames)*3+1 0 10])
    title('maxium fold changes','FontSize',18);
end

saveas(h,[dirname,'_',file_title,'.fig']);
saveas(h,[dirname,'_',file_title,'.eps'],'psc2');


end

