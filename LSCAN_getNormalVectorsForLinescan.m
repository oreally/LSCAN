function [ vectorX, vectorY ] = LSCAN_getNormalVectorsForLinescan(fitresult,xPoints,yPoints,L)
%calculates normal vectors of length L for points (x,y) on a polar coordinate fitted curve and returns x y
%coordinates for those vectors in vectorX (2 x numberOfPoints) and VectorY (2 x numberOfPoints)
%   Ortrud Wartlick, 2013

theta = [0:0.0001:2*pi+0.0001];
[THETAL, RHOL] = cart2pol(xPoints,yPoints);    
[xOutline,yOutline] = pol2cart(theta(:),fitresult(theta(:))); 

% get tangent vector points
[tangentpoints] = LSCAN_getTangents(fitresult,THETAL);

% get normal vector points
vectorX = nan(2,length(THETAL));
vectorY = nan(2,length(THETAL));
for n = 1:length(THETAL)
    % get tangent vector of length L
    [x1, y1] = pol2cart(tangentpoints{n}(:,1),tangentpoints{n}(:,2),'k');
    tangentvectorY = L*diff(y1)/sqrt(diff(x1)^2+diff(y1)^2);
    tangentvectorX = L*diff(x1)/sqrt(diff(x1)^2+diff(y1)^2);

    % get normal vector of length L
    normalvectorY=tangentvectorX;
    normalvectorX=-tangentvectorY;
    
    % check if endpoint-to-be is inside the cell (1) or outside (0)
    if inpolygon(xPoints(n)+normalvectorX*0.1,yPoints(n)+normalvectorY*0.1,xOutline,yOutline)==1
        normalvectorY=-tangentvectorX;
        normalvectorX=tangentvectorY;
    end
    
    % store vector points
    vectorX(:,n) = [xPoints(n)-normalvectorX*0.1 xPoints(n)+normalvectorX];
    vectorY(:,n) = [yPoints(n)-normalvectorY*0.1 yPoints(n)+normalvectorY];
end


