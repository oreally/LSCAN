function [Parameters] = LSCAN_getMembraneParametersBlebs(Variables, filenametext)
smoothy = 3;
SurfaceArea = Variables.SurfaceArea;
Linescan = Variables.Linescan;
dt = Variables.dt;

% estimate total CAAX amount in the bleb:
S1 = smooth(SurfaceArea(:,1),smoothy)./nanmean(smooth(SurfaceArea(:,1),smoothy));
Mtot1 = smooth(smooth(Linescan.IntegratedIntensity(:,1),smoothy).*smooth(SurfaceArea(:,1),smoothy),smoothy)./nanmean(smooth(smooth(Linescan.IntegratedIntensity(:,1),smoothy).*smooth(SurfaceArea(:,1),smoothy),smoothy));

% period of bleb retraction:
[maxS, from] = max(S1);
[c, to] = min(abs(S1(from:end)-0.1*maxS));
to = from+to(1)-1;

folding = Mtot1(to)/Mtot1(from);
flow = 1- folding;

h1=figure;
file_title='_membrane';

subplot(1,2,1)
plot([0:length(SurfaceArea)-1]*dt,S1,'k',...
    [0:length(SurfaceArea)-1]*dt,Mtot1,'b')
legend('Bleb: Apparent Surface','Bleb: Total PM-CAAX');
xlabel('time (s)');
hold on
plot([from from]*dt-dt,[0 3],'k--',...
    [to to]*dt-dt,[0 3],'k--')

subplot(1,2,2)
bar([folding,flow,0;0 0 0],'stacked');
axis([0 2 0 1.2])

saveas(h1,[filenametext,file_title,'.fig'])
saveas(h1,[filenametext,file_title,'.eps'],'psc2')

% 
% %divide by dt to get rates in s-1
% IdotI1=1/dt*smooth(diff(smooth(Linescan.MaximumIntensity(:,1,1),smoothy))./smooth(Linescan.MaximumIntensity(1:end-1,1,1),smoothy),smoothy);
% wdotw1=1/dt*smooth(diff(smooth(Linescan.Width(:,1),smoothy))./smooth(Linescan.Width(1:end-1,1),smoothy),smoothy);
% SdotS1=1/dt*smooth(diff(smooth(SurfaceArea(:,1),smoothy))./smooth(SurfaceArea(1:end-1,1),smoothy),smoothy);
% IdotI2=1/dt*smooth(diff(smooth(Linescan.MaximumIntensity(:,2,1),smoothy))./smooth(Linescan.MaximumIntensity(1:end-1,2,1),smoothy),smoothy);
% wdotw2=1/dt*smooth(diff(smooth(Linescan.Width(:,2),smoothy))./smooth(Linescan.Width(1:end-1,2),smoothy),smoothy);
% SdotS2=1/dt*smooth(diff(smooth(SurfaceArea(:,2),smoothy))./smooth(SurfaceArea(1:end-1,2),smoothy),smoothy);
% 
% % width Ch1
% w1 = smooth(Linescan.Width(:,1),smoothy)./nanmean(smooth(Linescan.Width(:,1),smoothy));
% w2 = smooth(Linescan.Width(:,2),smoothy)./nanmean(smooth(Linescan.Width(:,2),smoothy));
% % Imax Ch1
% Imax1 = smooth(Linescan.MaximumIntensity(:,1,1),smoothy)./nanmean(smooth(Linescan.MaximumIntensity(:,1,1),smoothy));
% Imax2 = smooth(Linescan.MaximumIntensity(:,2,1),smoothy)./nanmean(smooth(Linescan.MaximumIntensity(:,2,1),smoothy));
% % surface area
% S1 = smooth(SurfaceArea(:,1),smoothy)./nanmean(smooth(SurfaceArea(:,1),smoothy));
% S2 = smooth(SurfaceArea(:,2),smoothy)./nanmean(smooth( SurfaceArea(:,2),smoothy));
% % total membrane
% Mtot1 = smooth(smooth(Linescan.IntegratedIntensity(:,1),smoothy).*smooth(SurfaceArea(:,1),smoothy),smoothy)./nanmean(smooth(smooth(Linescan.IntegratedIntensity(:,1),smoothy).*smooth(SurfaceArea(:,1),smoothy),smoothy));
% Mtot2 = smooth(smooth(Linescan.IntegratedIntensity(:,2),smoothy).*smooth(SurfaceArea(:,2),smoothy),smoothy)./nanmean(smooth(smooth(Linescan.IntegratedIntensity(:,2),smoothy).*smooth(SurfaceArea(:,2),smoothy),smoothy));
% 
% 
% h1=figure;
% % x in frame units, not absolute time:
% x = 1:length(w1);
% plot(x,w1,'g',x,Imax1,'r',x,S1,'b')
% % get maxima and minima positions LHS
% [LHS_Imax_maxima,LHS_I_max] = ginput; % IMax maxima
% [LHS_Imax_minima,LHS_I_min] = ginput; % Imax minima
% [LHS_w_maxima,LHS_w_max] = ginput; % w maxima
% [LHS_w_minima,LHS_w_min] = ginput; % w minima
% [LHS_S_maxima,LHS_S_max] = ginput; % S maxima
% [LHS_S_minima,LHS_S_min] = ginput; % S minima
% LHS_Mtot_cV = Mtot1(round(LHS_S_minima));% Mtot in contracted pole
% LHS_Mtot_eV = Mtot1(round(LHS_S_maxima)); % Mtot in expanded pole
% 
% plot(x,w2,'g',x,Imax2,'r',x,S2,'b')
% % get maxima and minima positions RHS
% [RHS_Imax_maxima,RHS_I_max] = ginput; % IMax maxima
% [RHS_Imax_minima,RHS_I_min] = ginput; % Imax minima
% [RHS_w_maxima,RHS_w_max] = ginput; % w maxima
% [RHS_w_minima,RHS_w_min] = ginput; % w minima
% [RHS_S_maxima,RHS_S_max] = ginput; % S maxima
% [RHS_S_minima,RHS_S_min] = ginput; % S minima
% RHS_Mtot_cV = Mtot2(round(RHS_S_minima));% Mtot in contracted pole
% RHS_Mtot_eV = Mtot2(round(RHS_S_maxima)); % Mtot in expanded pole
% 
% xp = 1:length(wdotw1);
% plot(xp,wdotw1,'g',xp,IdotI1,'r',xp,SdotS1,'b')
% axis([0 length(xp) -0.04 0.04])
% % get maxima and minima dot positions LHS
% [LHS_IdotI_maxima,LHS_IdotI_max] = ginput; % IMax maxima
% [LHS_IdotI_minima,LHS_IdotI_min] = ginput; % Imax minima
% [LHS_wdotw_maxima,LHS_wdotw_max] = ginput; % w maxima
% [LHS_wdotw_minima,LHS_wdotw_min] = ginput; % w minima
% [LHS_SdotS_maxima,LHS_SdotS_max] = ginput; % S maxima
% [LHS_SdotS_minima,LHS_SdotS_min] = ginput; % S minima
% 
% plot(xp,wdotw2,'g',xp,IdotI2,'r',xp,SdotS2,'b')
% axis([0 length(xp) -0.04 0.04])
% % get maxima and minima dot positions RHS
% [RHS_IdotI_maxima,RHS_IdotI_max] = ginput; % IMax maxima
% [RHS_IdotI_minima,RHS_IdotI_min] = ginput; % Imax minima
% [RHS_wdotw_maxima,RHS_wdotw_max] = ginput; % w maxima
% [RHS_wdotw_minima,RHS_wdotw_min] = ginput; % w minima
% [RHS_SdotS_maxima,RHS_SdotS_max] = ginput; % S maxima
% [RHS_SdotS_minima,RHS_SdotS_min] = ginput; % S minima
% 
% close(h1);
% 
% % parameters:
% % peakshift of I and w around full contraction:
% Parameters.Linescan.ImaxShiftRelativeToFullContraction(1,:) = dt*[mean(LHS_Imax_maxima - LHS_S_minima), std(LHS_Imax_maxima - LHS_S_minima)];
% Parameters.Linescan.WidthShiftRelativeToFullContraction(1,:) = dt*[mean(LHS_w_maxima - LHS_S_minima), std(LHS_w_maxima - LHS_S_minima)];
% Parameters.Linescan.ImaxShiftRelativeToFullContraction(2,:) = dt*[mean(RHS_Imax_maxima - RHS_S_minima), std(RHS_Imax_maxima - RHS_S_minima)];
% Parameters.Linescan.WidthShiftRelativeToFullContraction(2,:) = dt*[mean(RHS_w_maxima - RHS_S_minima), std(RHS_w_maxima - RHS_S_minima)];
% % peakshift of I and w around full expansion:
% Parameters.Linescan.ImaxShiftRelativeToFullExpansion(1,:) = dt*[mean(LHS_Imax_minima - LHS_S_maxima), std(LHS_Imax_minima - LHS_S_maxima)];
% Parameters.Linescan.WidthShiftRelativeToFullExpansion(1,:) = dt*[mean(LHS_w_minima - LHS_S_maxima), std(LHS_w_minima - LHS_S_maxima)];
% Parameters.Linescan.ImaxShiftRelativeToFullExpansion(2,:) = dt*[mean(RHS_Imax_minima - RHS_S_maxima), std(RHS_Imax_minima - RHS_S_maxima)];
% Parameters.Linescan.WidthShiftRelativeToFullExpansion(2,:) = dt*[mean(RHS_w_minima - RHS_S_maxima), std(RHS_w_minima - RHS_S_maxima)];
% % percentage flow into contracted / expanded pole
% Parameters.Linescan.percTotalMembraneWhenContracted(1,:) = [mean(LHS_Mtot_cV) std(LHS_Mtot_cV)];
% Parameters.Linescan.percTotalMembraneWhenContracted(2,:) = [mean(RHS_Mtot_cV) std(RHS_Mtot_cV)];
% Parameters.Linescan.percTotalMembraneWhenExpanded(1,:) = [mean(LHS_Mtot_eV) std(LHS_Mtot_eV)];
% Parameters.Linescan.percTotalMembraneWhenExpanded(2,:) = [mean(RHS_Mtot_eV) std(RHS_Mtot_eV)];
% 
% % maximum and minimum rates of folding (i.e. folding = maximum and unfolding = minimum)
% Parameters.Linescan.maximumFoldingRate(1,:) = [mean(LHS_wdotw_max), std(LHS_wdotw_max)];
% Parameters.Linescan.maximumUnFoldingRate(1,:) = [mean(LHS_wdotw_min), std(LHS_wdotw_min)];
% Parameters.Linescan.maximumFoldingRate(2,:) = [mean(RHS_wdotw_max), std(RHS_wdotw_max)];
% Parameters.Linescan.maximumUnFoldingRate(2,:) = [mean(RHS_wdotw_min), std(RHS_wdotw_min)];
% % maximum ruffling and unruffling rates
% Parameters.Linescan.maximumRufflingRate(1,:) = [mean(LHS_IdotI_max), std(LHS_IdotI_max)];
% Parameters.Linescan.maximumUnRufflingRate(1,:) = [mean(LHS_IdotI_min), std(LHS_IdotI_min)];
% Parameters.Linescan.maximumRufflingRate(2,:) = [mean(RHS_IdotI_max), std(RHS_IdotI_max)];
% Parameters.Linescan.maximumUnRufflingRate(2,:) = [mean(RHS_IdotI_min), std(RHS_IdotI_min)];
% % maximum contraction and expansion rates
% Parameters.Linescan.maximumContractionRate(1,:) = [mean(LHS_SdotS_min), std(LHS_SdotS_min)];
% Parameters.Linescan.maximumExpansionRate(1,:) = [mean(LHS_SdotS_max), std(LHS_SdotS_max)];
% Parameters.Linescan.maximumContractionRate(2,:) = [mean(RHS_SdotS_min), std(RHS_SdotS_min)];
% Parameters.Linescan.maximumExpansionRate(2,:) = [mean(RHS_SdotS_max), std(RHS_SdotS_max)];
% 
% % maximum ruffling / folding / surface change
% Parameters.Linescan.maximumRuffling(1,:) = [mean(LHS_I_max./LHS_I_min) std(LHS_I_max./LHS_I_min)];
% Parameters.Linescan.maximumRuffling(2,:) = [mean(RHS_I_max./RHS_I_min) std(RHS_I_max./RHS_I_min)];
% Parameters.Linescan.maximumWidthChange(1,:) = [mean(LHS_w_max./LHS_w_min) std(LHS_w_max./LHS_w_min)];
% Parameters.Linescan.maximumWidthChange(2,:) = [mean(RHS_w_max./RHS_w_min) std(RHS_w_max./RHS_w_min)];
% Parameters.Linescan.maximumSurfaceChange(1,:) = [mean(LHS_S_max./LHS_S_min) std(LHS_S_max./LHS_S_min)];
% Parameters.Linescan.maximumSurfaceChange(2,:) = [mean(RHS_S_max./RHS_S_min) std(RHS_S_max./RHS_S_min)];
% 
% % period
% if length(LHS_S_maxima)>1
%     Parameters.Linescan.period(1,:) = dt*[mean([diff(LHS_S_maxima); diff(LHS_S_minima)]), std([diff(LHS_S_maxima); diff(LHS_S_minima)])];
%     Parameters.Linescan.period(2,:) = dt*[mean([diff(RHS_S_maxima); diff(RHS_S_minima)]), std([diff(RHS_S_maxima); diff(RHS_S_minima)])];
% else
%     Parameters.Linescan.period(1,:)=dt*abs(LHS_S_maxima-LHS_S_minima)*2;
%     Parameters.Linescan.period(2,:)=dt*abs(RHS_S_maxima-RHS_S_minima)*2;
% end
% 
% file_title = '_linescan_plots';
% h = figure('name',file_title,'numbertitle','off');
% 
% subplot(2,2,1)
% x = [0:length(SurfaceArea)-1]*dt;
% plot(x,smooth(SurfaceArea(:,1),smoothy)./mean(smooth(SurfaceArea(:,1),smoothy)),'b',...
%     x,smooth(Linescan.MaximumIntensity(:,1,1),smoothy)./mean(smooth(Linescan.MaximumIntensity(:,1,1),smoothy)),'r',...
%     x,smooth(Linescan.Width(:,1),smoothy)./mean(smooth(Linescan.Width(:,1),smoothy)),'g',...
%     x,smooth(Linescan.IntegratedIntensity(:,1),smoothy)./mean(smooth(Linescan.IntegratedIntensity(:,1),smoothy)),'k',...
%     x,smooth(Linescan.IntegratedIntensity(:,1).*SurfaceArea(:,1),smoothy)./mean(smooth(Linescan.IntegratedIntensity(:,1).*SurfaceArea(:,1),smoothy)),'k--')
% axis([0 length(SurfaceArea(:,1))*dt 0 2])
% xlabel('time (s)');
% legend('SurfaceArea 1','maximumIntensity','width','integratedIntensity','total membrane');
% %plot selected peaks on top
% hold on
% plot(LHS_Imax_maxima*dt,LHS_I_max,'ro',LHS_Imax_minima*dt,LHS_I_min,'ro',...
%     LHS_w_maxima*dt,LHS_w_max,'go',LHS_w_minima*dt,LHS_w_min,'go',...
%     LHS_S_maxima*dt,LHS_S_max,'bo',LHS_S_minima*dt,LHS_S_min,'bo')
% 
% 
% subplot(2,2,2) % pole 2
% x = [0:length(SurfaceArea(:,2))-1]*dt;
% plot(x,smooth(SurfaceArea(:,2),20)./mean(smooth(SurfaceArea(:,2),20)),'b',...
%     x,smooth(Linescan.MaximumIntensity(:,2,1),20)./mean(smooth(Linescan.MaximumIntensity(:,2,1),20)),'r',...
%     x,smooth(Linescan.Width(:,2),20)./mean(smooth(Linescan.Width(:,2),20)),'g',...
%     x,smooth(Linescan.IntegratedIntensity(:,2),20)./mean(smooth(Linescan.IntegratedIntensity(:,2),20)),'k',...
%     x,smooth(Linescan.IntegratedIntensity(:,2).*SurfaceArea(:,2),20)./mean(smooth(Linescan.IntegratedIntensity(:,2).*SurfaceArea(:,2),20)),'k--')
% axis([0 length(SurfaceArea(:,2))*dt 0 2])
% xlabel('time (s)');
% legend('SurfaceArea 2','maximumIntensity','width','integratedIntensity','total membrane');
% %plot selected peaks on top
% hold on
% plot(RHS_Imax_maxima*dt,RHS_I_max,'ro',RHS_Imax_minima*dt,RHS_I_min,'ro',...
%     RHS_w_maxima*dt,RHS_w_max,'go',RHS_w_minima*dt,RHS_w_min,'go',...
%     RHS_S_maxima*dt,RHS_S_max,'bo',RHS_S_minima*dt,RHS_S_min,'bo')
% 
% 
% subplot(2,2,3)
% x = [0:length(IdotI1)-1]*dt;
% plot(x,IdotI1,'r',...
%     x,wdotw1,'g',...
%     x,-SdotS1,'b',...
%     x,IdotI1+wdotw1,'k')
% axis([0 length(IdotI1)*dt -0.02 0.02])
% xlabel('time (s)');
% legend('Imax./Imax','width./width','-S./S','Imax./Imax + width./width');
% %plot selected peaks on top
% hold on
% plot(LHS_IdotI_maxima*dt,LHS_IdotI_max,'ro',LHS_IdotI_minima*dt,LHS_IdotI_min,'ro',...
%     LHS_wdotw_maxima*dt,LHS_wdotw_max,'go',LHS_wdotw_minima*dt,LHS_wdotw_min,'go',...
%     LHS_SdotS_maxima*dt,-LHS_SdotS_max,'bo',LHS_SdotS_minima*dt,-LHS_SdotS_min,'bo')
% 
% subplot(2,2,4)
% x = [0:length(IdotI2)-1]*dt;
% plot(x,IdotI2,'r',...
%     x,wdotw2,'g',...
%     x,-SdotS2,'b',...
%     x,IdotI2+wdotw2,'k')
% axis([0 length(IdotI2)*dt -0.02 0.02])
% xlabel('time (s)');
% legend('Imax./Imax','width./width','-S./S','Imax./Imax + width./width');
% %plot selected peaks on top
% hold on
% plot(RHS_IdotI_maxima*dt,RHS_IdotI_max,'ro',RHS_IdotI_minima*dt,RHS_IdotI_min,'ro',...
%     RHS_wdotw_maxima*dt,RHS_wdotw_max,'go',RHS_wdotw_minima*dt,RHS_wdotw_min,'go',...
%     RHS_SdotS_maxima*dt,-RHS_SdotS_max,'bo',RHS_SdotS_minima*dt,-RHS_SdotS_min,'bo')
% 
% saveas(h,[filenametext,file_title,'.fig'])
% saveas(h,[filenametext,file_title,'.eps'],'psc2')
% 
% Parameters.Linescan
% end