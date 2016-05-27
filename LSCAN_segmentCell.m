function BW = LSCAN_segmentCell(image,segmentationThreshold)


image = medfilt2(image,[3 3]);
BW = im2bw(image,segmentationThreshold); % or one segmentation parameter
BW=bwmorph(BW,'dilate',3);
BW=~BW;
BW = imclearborder(BW);
BW=bwmorph(BW,'dilate',4);
BW=bwmorph(BW,'thicken',4); %5
BW=bwmorph(BW,'bridge',1);
BW = imfill(BW, 'holes');
BW=bwmorph(BW,'erode',4); %4
BW=bwareaopen(BW,300);

% add cell pole bridge
% BW1 = bwmorph(BW,'dilate',15);
% BW1 = bwmorph(BW1,'erode',17);
% BW(BW1==1)=1;

% when signal is not so good, try this:
% image = medfilt2(image,[3 3]);
% BW = im2bw(image,segmentationThreshold); % or one segmentation parameter
% BW=bwmorph(BW,'dilate',5);
% BW=~BW;
% BW=bwmorph(BW,'dilate',2);
% BW = imclearborder(BW);
% BW = imfill(BW, 'holes');
% BW=bwmorph(BW,'dilate',1);
% BW = imfill(BW, 'holes');


%active contour segmentation: (not as fast)
% hold on
% contour(BW,[0 0],'r');
% pause(0.1);
% mask = BW;
% BW = activecontour(image,mask,300,'edge');
% BW=bwmorph(BW,'erode',2); %4


