function [xlinescan, ylinescan] = LSCAN_getLinescanPoints(intersectionX,intersectionY,innerFit, from,to,distance)
% gets points on a polar coordinate curve at fixed distances from each other
% and returns the points in cartesian coordinates
% Ortrud Wartlick 2013

xlinescan(1) = intersectionX(from);
ylinescan(1) = intersectionY(from);
[tlinescan, rlinescan] = cart2pol(xlinescan(1),ylinescan(1));
[tlinescanEnd, rlinescanEnd] = cart2pol(intersectionX(to),intersectionY(to));

% account for matlab's weird way of doing angles
if tlinescan<0,
    tlinescan = pi+(pi+tlinescan);
end
if tlinescanEnd<0
    tlinescanEnd = pi+(pi+tlinescanEnd);
end

if tlinescan >= tlinescanEnd 
    tlinescanEnd = tlinescanEnd+2*pi;
end


% make sure the start and end are right
if abs(tlinescanEnd-tlinescan)<pi % then it's the small pie
    temp = tlinescan;
    tlinescan = tlinescanEnd;
    tlinescanEnd = temp;
end

if abs(tlinescanEnd-tlinescan)<pi % then it's the small pie
    tlinescanEnd = tlinescanEnd+2*pi;
end

xlinescan = nan(10000,1);
ylinescan = nan(10000,1);
for n=2:10000
    % calculate d for different possible thetas
    t = tlinescan(n-1)+[0:0.00001:pi/12]'; % march a little distance from tlinescan(n-1) around the fit
    r = innerFit(t);
    [xtest, ytest]=pol2cart(t,r);
    
    % find position on curve where d = distance
    xx = diff(xtest); yy = diff(ytest); % in preparation for d calculation
    d = cumsum(sqrt((xx).^2+(yy).^2))-sqrt((xx(1)).^2+(yy(1)).^2)-distance;
    [m, pos] = min(abs(d)); % d = distance
    
    % store new point
    tlinescan(n) = t(pos);
    xlinescan(n) = xtest(pos);
    ylinescan(n) = ytest(pos);
    if tlinescan(n)>=tlinescanEnd
        break
    end
    clear t r d xtest ytest pos
end
xlinescan(isnan(xlinescan))=[];
ylinescan(isnan(ylinescan))=[];
