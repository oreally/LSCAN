function parameters = LSCAN_mech_analyseOsciCells2(File,LHSactinCortexMeanNorm,RHSactinCortexMeanNorm,geometry,fromFrame,toFrame,pix,dt)
%
% version date: Fenuary 19, 2013 / June 28, 2013
% created by: Jochen Schneider 
% (modified slightly for presentation purposes by Ortrud Wartlick)
%
% This function performs the actual analysis of cells. 
% There needs to be a folder named 'images'. In this folder it will store figures for
% turnover, tension, friction and the phase diagram. The figures are stored
% in both types 'eps' and 'fig'. 


% saturation value in myosin function
cstar = 0.8;

% cross section areas; spatial units = um
LHSRadius = geometry.Radius(:,1);
RHSRadius = geometry.Radius(:,2);
LHSVolume = geometry.Volume(:,1);
RHSVolume = geometry.Volume(:,2);
LHSSurfaceArea = geometry.SurfaceArea(:,1);
RHSSurfaceArea = geometry.SurfaceArea(:,2);
v=(LHSVolume-RHSVolume)./(LHSVolume+RHSVolume);
smoothy = 20;



% initialize DATA matrix
parameters = zeros(1,16);

% calculating actin mean values
actinCortexMean = mean(LHSactinCortexMeanNorm+RHSactinCortexMeanNorm)/2;


%%%%%%%%%%%%%%%
%%% FITTING %%%
%%%%%%%%%%%%%%%

% calculate derivative of volume
dv=1/dt*[diff(smooth(v,smoothy));0];%divide by dt to get velocity in s-1
v = v(fromFrame:toFrame); v = v(:);
dv = dv(fromFrame:toFrame);

% already in um
R0=(3/4/pi*(LHSVolume(fromFrame:toFrame)./2+RHSVolume(fromFrame:toFrame)./2)).^(1/3);
Rc = geometry.FurrowWidth(fromFrame:toFrame,1)./(2*R0);
nanmean(Rc)

% function giving myosin activity as a function of the intensity of actin
f= @(x) x./(x+cstar)*(1+cstar);



%%% Fitting volume equation %%%

% pressure P = T / R
ActinActivePressure=R0.*(f(LHSactinCortexMeanNorm(fromFrame:toFrame)./actinCortexMean))./(LHSRadius(fromFrame:toFrame))...
    -R0.*(f(RHSactinCortexMeanNorm(fromFrame:toFrame)./actinCortexMean))./(RHSRadius(fromFrame:toFrame));
ActinActivePressure = ActinActivePressure(:);

X=[ActinActivePressure v v.^3];

% perform regression
stats=regstats(dv,X);
P0m=stats.beta;
varP0m=stats.covb;
Rsquare=stats.rsquare;

% tension divided by R0 and K
T=P0m(2)/P0m(3);stdT=sqrt(varP0m(2,2)/P0m(3)^2+varP0m(3,3)*P0m(2)^2/P0m(3)^4);

% friction divided by 2K
alpha=-1/P0m(3);stdalpha=sqrt(varP0m(3,3)/P0m(3)^4);

% K3 divided by K
K3=P0m(4)/P0m(3);stdK3=sqrt(varP0m(4,4)/P0m(3)^2+varP0m(3,3)*P0m(4)^2/P0m(3)^4);

% print results to Command Window
disp(['T/R0K=' num2str(T) ' pm ' num2str(stdT) ', alpha/2K=' num2str(alpha) ' pm ' ...
    num2str(stdalpha) ' K3/K= ' num2str(K3) ' pm ' num2str(stdK3)]);

% save result in variable DATA
parameters(1:6) = [T,stdT,alpha,stdalpha,K3,stdK3];

% returns a straight line with the slope -T
T_data_x = ActinActivePressure;
T_data_y = alpha*dv + v + K3*v.^3;

% returns a straight line with the slope alpha
alpha_data_x = dv(:);
alpha_data_y = -T*ActinActivePressure-v-K3*v.^3;

% elastic curve
elastic_data_x = v(:);
elastic_data_y = alpha*dv+T*ActinActivePressure;

maxv=max(T_data_x);minv=min(T_data_x);
pdistrib=[minv:abs(maxv-minv)/30:maxv];
T_fit_x = pdistrib;
T_fit_y = -T*pdistrib;

maxv=max(alpha_data_x);minv=min(alpha_data_x);
dvdistrib=[minv:abs(maxv-minv)/30:maxv];
alpha_fit_x = dvdistrib;
alpha_fit_y = alpha*dvdistrib + nanmean(alpha_data_y);

maxv=max(elastic_data_x);minv=min(elastic_data_x);
vdistrib=[minv:abs(maxv-minv)/30:maxv];
elastic_fit_x = vdistrib;
elastic_fit_y = -(vdistrib+K3*vdistrib.^3);


%%% Fitting actin equation %%%

% LHS of actin equation for both actin and myosin
LHSactinCortexTotal=LHSactinCortexMeanNorm.*LHSSurfaceArea;
RHSactinCortexTotal=RHSactinCortexMeanNorm.*RHSSurfaceArea;
dtLHSactinCortexMeanNorm=(1/dt)*[diff(smooth(LHSactinCortexTotal,smoothy));0]./LHSSurfaceArea/actinCortexMean;%divide by dt to get velocity in 1/s
dtRHSactinCortexMeanNorm=(1/dt)*[diff(smooth(RHSactinCortexTotal,smoothy));0]./RHSSurfaceArea/actinCortexMean;

dT=[dtLHSactinCortexMeanNorm(fromFrame:toFrame);dtRHSactinCortexMeanNorm(fromFrame:toFrame)];
T=[LHSactinCortexMeanNorm(fromFrame:toFrame)/actinCortexMean,RHSactinCortexMeanNorm(fromFrame:toFrame)/actinCortexMean];
T = T(:);

% perform regression
LHRHv=[v(:); -v(:)];
X=[T LHRHv LHRHv.^2];
stats=regstats(dT,X);
P0m=stats.beta;
varP0m=stats.covb;
Rsquare=stats.rsquare;

% c0
c0=-P0m(1)/P0m(2);stdc0=sqrt(varP0m(1,1)/P0m(2)^2+varP0m(2,2)*P0m(1)^2/P0m(2)^4);

% inverse turnover
kd=-P0m(2);stdkd=sqrt(varP0m(2,2));%a positive value indicates negative relaxation

Cv=P0m(3);stdCv=sqrt(varP0m(3,3));
Cv2=P0m(4);stdCv2=sqrt(varP0m(4,4));

% print results to Command Window
disp(['1/tau=' num2str(kd) ' pm ' num2str(stdkd) ', c0=' num2str(c0) ' pm ' ...
    num2str(stdc0) ' Cv=' num2str(Cv) ' pm ' num2str(stdCv)]);
% save results to DATA variable
parameters(7:12) = [kd,stdkd,c0,stdc0,Cv,stdCv];

turn_data_x = T;
turn_data_y = dT-Cv*LHRHv-Cv2*LHRHv.^2;
maxT=max(turn_data_x);minT=min(turn_data_x);
Tdistrib=[minT:abs(maxv-minv)/30:maxT];
turn_fit_x = Tdistrib;
turn_fit_y = kd*(c0-Tdistrib);

clear X stats P0m varP0m Rsquare vdistrib dvdistribut pdistribut dT T LHRHv Tdistrib

% phase diagram parameter
parameters(13) = parameters(3).*parameters(7);
parameters(14) = abs(parameters(3)).*parameters(8) + abs(parameters(7)).*parameters(4);
parameters(15) = parameters(1);
parameters(16) = parameters(2);


%%%%%%%%%%%%%
%%% PLOTS %%%
%%%%%%%%%%%%%

file_title = 'TTEfits';

h = figure('name',file_title,'numbertitle','off');

subplot(2,2,1)
endfr = length(LHSactinCortexMeanNorm);
c = linspace(1+fromFrame/endfr*10,toFrame/endfr*10,length(T_data_x));
c(1)=1; c(end)=10;
scatter(T_data_x,T_data_y,16,c,'o','filled');hold on;
xmax = max(abs(T_data_x));
ymax = max(abs(T_data_y));
xlim([-1.2*xmax 1.2*xmax]);
ylim([-1.2*ymax,1.2*ymax]);
plot(T_fit_x,T_fit_y,'red','LineWidth',2);hold on;
title('Tension', 'FontWeight', 'bold','FontSize',16);
set(gca,'FontSize',14);
xlabel('pressure','FontWeight','bold','FontSize',14);
%ylabel('f(v)','FontWeight','bold','FontSize',14);

subplot(2,2,2)
endfr = length(LHSactinCortexMeanNorm);
c = linspace(1+fromFrame/endfr*10,toFrame/endfr*10,length(alpha_data_x));
c(1)=1; c(end)=10;
scatter(alpha_data_x,alpha_data_y,16,c,'o','filled');hold on;
xmax = max(abs(alpha_data_x));
ymax = max(abs(alpha_data_y));
xlim([-1.2*xmax 1.2*xmax]);
ylim([-1.2*ymax,1.2*ymax]);
plot(alpha_fit_x,alpha_fit_y,'red','LineWidth',2);hold on;
title('Friction', 'FontWeight', 'bold','FontSize',16);
set(gca,'FontSize',14);
xlabel('dv/dt','FontWeight','bold','FontSize',14);
%ylabel('f(v)','FontWeight','bold','FontSize',14);

subplot(2,2,3)
endfr = length(LHSactinCortexMeanNorm);
c = linspace(1+fromFrame/endfr*10,toFrame/endfr*10,length(elastic_data_x));
c(1)=1; c(end)=10;
scatter(elastic_data_x,elastic_data_y,16,c,'o','filled');hold on;
maxabsv = max(abs(elastic_data_x));
maxdv = max(abs(elastic_data_y));
xlim([-1.2*maxabsv 1.2*maxabsv]);
ylim([-1.2*maxdv,1.2*maxdv]);
plot(elastic_fit_x,elastic_fit_y,'red','LineWidth',2);hold on;
title('Cell Elastic Response', 'FontWeight', 'bold','FontSize',16);
set(gca,'FontSize',14);
xlabel('v','FontWeight','bold','FontSize',14);
ylabel('f(v)','FontWeight','bold','FontSize',14);

subplot(2,2,4)
endfr = length(LHSactinCortexMeanNorm);
c = linspace(1+fromFrame/endfr*10,toFrame/endfr*10,length(turn_data_x));
c(1)=1; c(end)=10;
scatter(turn_data_x,turn_data_y,16,c,'o','filled');hold on;
maxT=max(turn_data_x);
minT=min(turn_data_x);
xlim([0.7*minT 1.2*maxT]);
maxY = max(abs(turn_data_y));
ylim([-1.2*maxY 1.2*maxY]);
plot(turn_fit_x,turn_fit_y,'red','LineWidth',2);hold on;
title('Evolution of actin density', 'FontWeight', 'bold','FontSize',16);
set(gca,'FontSize',14);
xlabel('ca','FontWeight','bold','FontSize',14);
ylabel('1/S d(S ca)/dt','FontWeight','bold','FontSize',14);

saveas(h,[File.pathName,'figures/',File.fileName(1:end-4),'_',file_title,'.fig'])
saveas(h,[File.pathName,'figures/',File.fileName(1:end-4),'_',file_title,'.eps'],'psc2')
%close(h)

% plot phase diagrams
file_title = 'phasediagram';

h1 = figure('name',file_title,'numbertitle','off');
LSCAN_mech_phasediagram([parameters(13) parameters(14) parameters(15) parameters(16)],'blue');
saveas(h1,[File.pathName,'figures/',File.fileName(1:end-4),'_',file_title,'.fig'])
saveas(h1,[File.pathName,'figures/',File.fileName(1:end-4),'_',file_title,'.eps'],'psc2')
%close(h1)
