function LSCAN_concatenateResults_LACT( dirname )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% then get all the subfolders in this folder (that have LACT, in this case)
listing = dir(strcat(dirname,'/*LACT*'));
subfoldersi = [listing.isdir];
subfolderNames = {listing(subfoldersi).name};
% get all the .mat files in each of the subfolders
for f = 1:length(subfolderNames)
    mat_files{f} = dir(strcat(dirname,'/',subfolderNames{f},'/*.mat'));
    % concatenate data from all the mat files in a subfolder
    if ~isempty(mat_files{f})
        TRK{f} = nan(length(mat_files{f}),2);
        alpha2K{f} = nan(length(mat_files{f}),2);
        tauActin{f} = nan(length(mat_files{f}),2);
        for n=1:length(mat_files{f})
            mat_files{f}(n).name
            data = load(strcat(dirname,'/',subfolderNames{f},'/',mat_files{f}(n).name));
            % load TTE data:
            % parameters(1:6) = [T,stdT,alpha,stdalpha,K3,stdK3];
            % parameters(7:12) = [kd=1/tauActin,stdkd,c0,stdc0,Cv,stdCv];
%            if (data.Parameters.TTE(1)<1.7) && (data.Parameters.TTE(4)<5)
listdata(n,:) = [data.Parameters.TTE]; 
                TRK{f}(n,:)=data.Parameters.TTE(1:2);
                alpha2K{f}(n,:)=data.Parameters.TTE(3:4);
                tauActin{f}(n,:) = [1./data.Parameters.TTE(7) sqrt((1./data.Parameters.TTE(7)^2)^2*data.Parameters.TTE(8)^2)];
                clear data
%             else
%                 TRK{f}(n,:)=NaN;
%                 alpha2K{f}(n,:)=NaN;
%                 tauActin{f}(n,:) = NaN;
%             end
        end
    else
        listdata(1,:) = NaN; 
        TRK{f}(1,:)=NaN;
        alpha2K{f}(1,:)=NaN;
        tauActin{f}(1,:) = NaN;
    end
    save(strcat(dirname,'/',subfolderNames{f},'/results.txt'),'listdata','-ASCII');
    clear listdata
end

% make all the same length
for f=1:length(subfolderNames)
    TRK{f}(end+1:10,:)=NaN;
    alpha2K{f}(end+1:10,:)=NaN;
    tauActin{f}(end+1:10,:)=NaN;
end
for f=1:length(subfolderNames)
    % different genotypes in different columns
    TRKs(:,f)=TRK{f}(:,1);
    alpha2Ks(:,f)=alpha2K{f}(:,1);
    tauActins(:,f) = tauActin{f}(:,1);
end

h=figure('Position',[20 20 1000 700]);
file_title='results_tte';
style={'b','m','g','k','r','y','c','b','r','g','k','m'};
for n=1:size(TRKs,2)
    labels{n} = subfolderNames{n};
end

subplot(1,4,1)
boxplot(TRKs,'labels',labels,'labelorientation','inline');
for n=1:size(TRKs,2)
    hold on;plot((n-1)+0.9+0.2*rand(1,length(TRKs(:,n))),TRKs(:,n),strcat(style{n},'.'))
    axis([0 length(subfolderNames)+1 0 3])
    title('T/(RK)','FontSize',18);
end

subplot(1,4,2)
boxplot(alpha2Ks,'labels',labels,'labelorientation','inline');
for n=1:size(alpha2Ks,2)
    hold on;plot((n-1)+0.9+0.2*rand(1,length(alpha2Ks(:,n))),alpha2Ks(:,n),strcat(style{n},'.'))
    axis([0 length(subfolderNames)+1 0 60])
    title('\alpha /(2K)','FontSize',18);
end

subplot(1,4,3)
boxplot(tauActins,'labels',labels,'labelorientation','inline');
for n=1:size(tauActins,2)
    hold on;plot((n-1)+0.9+0.2*rand(1,length(tauActins(:,n))),tauActins(:,n),strcat(style{n},'.'))
    axis([0 length(subfolderNames)+1 0 100])
    title('\tau_{cortex}','FontSize',18);
end

% make phase diagram with means/stds of all points
subplot(1,4,4)
LSCAN_mech_phasediagram([NaN NaN NaN NaN],'blue');
for n=1:size(TRKs,2)
    hold on;
    % plot means and stdeviations
    plot(nanmean(alpha2Ks(:,n)./tauActins(:,n)),nanmean(TRKs(:,n)),strcat(style{n},'o'),'LineWidth',3)
    plot([nanmean(alpha2Ks(:,n)./tauActins(:,n))-nanstd(alpha2Ks(:,n)./tauActins(:,n)) nanmean(alpha2Ks(:,n)./tauActins(:,n))+nanstd(alpha2Ks(:,n)./tauActins(:,n))],[nanmean(TRKs(:,n)) nanmean(TRKs(:,n))],strcat(style{n},'-'),...
        [nanmean(alpha2Ks(:,n)./tauActins(:,n)) nanmean(alpha2Ks(:,n)./tauActins(:,n))],[nanmean(TRKs(:,n))-nanstd(TRKs(:,n)) nanmean(TRKs(:,n))+nanstd(TRKs(:,n))],strcat(style{n},'-'),'LineWidth',3)
    % plot individuals
%     plot(alpha2Ks(:,n)./tauActins(:,n),TRKs(:,n),strcat(style{n},'.'))
%     % plot ellipse
%     eccentricity = axes2ecc(nanstd(alpha2Ks(:,n)./tauActins(:,n)),nanstd(TRKs(:,n)));
%     [lat,lon] = ellipse1(nanmean(alpha2Ks(:,n)./tauActins(:,n)),nanmean(TRKs(:,n)),[nanstd(alpha2Ks(:,n)./tauActins(:,n)),eccentricity]);
%     plot(lat,lon,strcat(style{n},'-'));
end
title('Phasediagram','FontSize',18);

saveas(h,[dirname,'_',file_title,'.fig']);
saveas(h,[dirname,'_',file_title,'.eps'],'psc2');


end

