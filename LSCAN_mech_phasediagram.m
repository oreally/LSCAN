function h = LSCAN_mech_phasediagram(varargin)
% This function plots an arbitrary large number of points to the phase
% diagram obtained by the cell oscillation model introduced in
%
% doi:10.1038/nature10286
%
% The function requires an input in the following way
%
% phasediagram([x1 x1_e y1 y1_e],<'color1'>,[x2 x2_e y2 y2_e],<'color2'>,...)
%
% where [x1 x1_e y1 y1_e] is the vector of the specific point which should
% be plotted to the phase diagram (here the components marked with 'e' are 
% the errors) and <'color1'> are optional parameters which set the color in
% which this point should be plotted. 'color' can be any MATLAB built-in
% color.
%
% author: Jochen Schneider
% last update: Feb. 24, 2012


count3 = 1;
for count1 =1:length(varargin)
if ischar(varargin{count1}) == 1
  point_color{count3-1} = varargin{count1};
else
  point_color{count3} = 'black';
  point(count3,1:4) = varargin{count1};
  count3 = count3 + 1;
end
end

toscillant1 = @(x)((4.99283e-63-2.57867e-47*1i)*((3.12951e30+3.23262e46*1i)+...
    (5.53697e30+2.23154e45*1i)*x-(2.97526e23-2.32301e23*1i)*sqrt((-3.80309e45...
    -1.52124e46*1i)*x+(1.77028e45+7.08111e45*1i)*x.^2)));
toscillant2 = @(x)((4.99283e-63-2.57867e-47*1i)*((3.12951e30+3.23262e46*1i)+...
    (5.53697e30+2.23154e45*1i)*x + (2.97526e23-2.32301e23*1i)*...
    sqrt((-3.80309e45-1.52124e46*1i)*x+(1.77028e45+7.08111e45*1i)*x.^2)));
tinstable = @(x)((0.833587+8.06997e-17*1i)*(1+x));
tau1 = [0:0.01:1.15];
tau2 = [0:0.01:1.15];
tau3 = [0:0.01:1.15];
tau4 = [1.15:0.01:2.1];

h=area([0 2.1],[3 3],'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
hold on
area(tau1,real(toscillant1(tau1)),'EdgeColor','none','FaceColor','white')
hold on
area(tau3,real(tinstable(tau3)),'EdgeColor','none','FaceColor',[1 1 0.7])
hold on
area(tau4,1.79*ones(1,length(tau4)),'EdgeColor','none','FaceColor',[1 1 0.7])
hold on
plot(tau1,real(toscillant1(tau1)),'--',tau3,real(tinstable(tau3)),tau4,1.79*ones(1,length(tau4)),'Color','black');
hold on
if length(point_color) > 1
 for count2 = 2:length(point_color)
  LSCAN_mech_errorxy([point(count2,1),point(count2,2),point(count2,3),point(count2,4)],'ColX',1,'ColY',3,'ColXe',2,'ColYe',4,'EdgeColor',point_color{count2});
  hold on
 end
end
LSCAN_mech_errorxy([point(1,1),point(1,2),point(1,3),point(1,4)],'ColX',1,'ColY',3,'ColXe',2,'ColYe',4,'EdgeColor',point_color{1});
set(gca,'Layer','top')
hold off
xlim([0 2.1])
ylim([0 3])
axis square
title('Phase Diagram', 'FontWeight', 'bold','FontSize',16);
xlabel('\tau_{cell}/\tau_{cortex}','FontWeight','bold','FontSize',14);
ylabel('T/(R_0K)','FontWeight','bold','FontSize',14);
end
