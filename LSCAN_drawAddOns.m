function LSCAN_drawAddOns(show,addons, current)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
currentZoom = current.currentZoom;
currentC = current.currentC;
currentZ = current.currentZ;
currentT = current.currentT;
if show
    % if segmented, draw segmentation
    CellMask = addons.CellMask(:,:,1,currentZ,currentT);
    if sum(CellMask(:))>0
        style = {'r','m'};
        for k = 1:max(unique(CellMask))
            mask = zeros(size(CellMask));
            mask(CellMask==k)=1;
            hold on
            contour(mask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3)),[0 0],style{k});
            hold off
        end
    end
    
    % if outside (cortex height) is activated, draw outside
    CellOutside = addons.CellOutside(:,:,1,currentZ,currentT);
    if sum(CellOutside(:))>0
        style = {'r','m'};
        for k = 1:max(unique(CellOutside))
            mask = zeros(size(CellOutside));
            mask(CellOutside==k)=1;
            hold on
            contour(mask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3)),[0 0],style{k});
            hold off
        end
    end
    
    % if the background ROI has been selected
    backgPos = addons.backgroundROI;
    if ~isempty(backgPos)
        hold on
        plot([backgPos(1)-currentZoom(1),backgPos(1)-currentZoom(1)+backgPos(3),backgPos(1)-currentZoom(1)+backgPos(3),backgPos(1)-currentZoom(1),backgPos(1)-currentZoom(1)],...
            [backgPos(2)-currentZoom(2)+backgPos(4),backgPos(2)-currentZoom(2)+backgPos(4),backgPos(2)-currentZoom(2),backgPos(2)-currentZoom(2),backgPos(2)-currentZoom(2)+backgPos(4)],...
            'g');
        hold off
    end
    
    % if furrow has been selected
    centerpoints = addons.FurrowPlane{currentT};
    if ~isempty(centerpoints)
        try
            meanlength = centerpoints(currentZ,4)/2;
            x = [centerpoints(currentZ,1)-cos(atan(centerpoints(currentZ,3)))*meanlength centerpoints(currentZ,1)+cos(atan(centerpoints(currentZ,3)))*meanlength];
            y = [centerpoints(currentZ,2)-sin(atan(centerpoints(currentZ,3)))*meanlength centerpoints(currentZ,2)+sin(atan(centerpoints(currentZ,3)))*meanlength];
            hold on
            plot(x,y,'b','LineWidth',3)
            hold off
        catch
            'furrow drawing problem'
        end
    end
    
    % if inner circumference has been fit:
    if ~isempty(addons.innerCircFit{currentT,1})
        t = currentT;
        theta = [0:0.01*pi:2*pi,0];
        % Pole 1
        [X1 Y1] = pol2cart(theta(:),addons.innerCircFit{t,1}(theta(:)));
        BW1 = poly2mask(X1+addons.innerCircFit{t,3}(1),Y1+addons.innerCircFit{t,3}(2),current.imageHeightWidth(1),current.imageHeightWidth(2));
        % Pole 1
        [X2 Y2] = pol2cart(theta(:),addons.innerCircFit{t,2}(theta(:)));
        BW2 = poly2mask(X2+addons.innerCircFit{t,4}(1),Y2+addons.innerCircFit{t,4}(2),current.imageHeightWidth(1),current.imageHeightWidth(2));
        BW = double(BW1)+double(BW2);
        BW(BW>1)=1;
        
        hold on
        contour(BW(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3)),[0 0],'g');
        hold off
    end

end


