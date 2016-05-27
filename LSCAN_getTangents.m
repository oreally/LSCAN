function [tangentpoints] = LSCAN_getTangents(fit1,theta)

tangentpoints = cell([length(theta),1]);

for n = 1:length(theta)
    tangentpoints{n}(:,1) = [theta(n)-0.1/(2*pi) theta(n)+0.1/(2*pi)];
    tangentpoints{n}(:,2) = [fit1(theta(n)-0.1/(2*pi)) fit1(theta(n)+0.1/(2*pi))];
end
