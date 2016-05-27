function varargout = LSCAN_analysis(varargin)
% LSCAN_ANALYSIS MATLAB code for LSCAN_analysis.fig
%      LSCAN_ANALYSIS, by itself, creates a new LSCAN_ANALYSIS or raises the existing
%      singleton*.
%
%      H = LSCAN_ANALYSIS returns the handle to a new LSCAN_ANALYSIS or the handle to
%      the existing singleton*.
%
%      LSCAN_ANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LSCAN_ANALYSIS.M with the given input arguments.
%
%      LSCAN_ANALYSIS('Property','Value',...) creates a new LSCAN_ANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LSCAN_analysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LSCAN_analysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LSCAN_analysis

% Last Modified by GUIDE v2.5 13-May-2014 13:22:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LSCAN_analysis_OpeningFcn, ...
    'gui_OutputFcn',  @LSCAN_analysis_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before LSCAN_analysis is made visible.
function LSCAN_analysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LSCAN_analysis (see VARARGIN)

% Choose default command line output for LSCAN_analysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LSCAN_analysis wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LSCAN_analysis_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



% --------------------------------------------------------------------
function file_Callback(hObject, eventdata, handles)
% hObject    handle to file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function loadStack_Callback(hObject, eventdata, handles)
% hObject    handle to loadStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile({'.mat';'.tif'},'Select a multitif or a .mat file');

% if it is a tif file, load the images and initialise variables
if strcmp(FileName(length(FileName)-2:length(FileName)),'tif')
    % get image info in cell array:
    %     Channels: 2
    %     Slices: 30
    %     Frames: 41
    %     Width
    %     Height
    %     PixelX: 0.2412076
    %     PixelY: 0.2412076
    %     PixelZ: 1
    %     FrameInterval: 20
    %     Channel1: name
    %     Channel2: name
    fileID = fopen(strcat(PathName,FileName(1:end-4),'.txt'));
    ImageInfo = textscan(fileID, '%s %s');
    fclose(fileID);
    
    % load xyczt stack:
    numChannels = str2num(ImageInfo{2}{1});
    numZSlices = str2num(ImageInfo{2}{2});
    numTFrames = str2num(ImageInfo{2}{3});
    imageWidth = str2num(ImageInfo{2}{4});
    imageHeight = str2num(ImageInfo{2}{5});
    ImageStack = zeros([imageHeight,imageWidth,numChannels,numZSlices,numTFrames],'uint8');
    
    % load the images: xyczt
    for t = 1:numTFrames
        for z = 1:numZSlices
            for c = 1:numChannels
                ImageStack(:,:,c,z,t) = imread(strcat(PathName,FileName),(t-1)*numChannels*numZSlices+(c*z+(numChannels-c)*(z-1)));
            end
        end
    end
    
    handles.File.fileName = FileName;
    handles.File.pathName = PathName;
    handles.File.ImageStack = ImageStack;
    handles.File.ImageInfo = ImageInfo;
    
    handles.AddOns.excludedFrames = zeros(1,numTFrames);
    handles.AddOns.CellMask = zeros([size(ImageStack,1),size(ImageStack,2),1,size(ImageStack,4),size(ImageStack,5)]);
    handles.AddOns.CellOutside = zeros([size(ImageStack,1),size(ImageStack,2),1,size(ImageStack,4),size(ImageStack,5)]);
    handles.AddOns.ExternalMask = zeros([size(ImageStack,1),size(ImageStack,2),1,size(ImageStack,4),size(ImageStack,5)]);
    handles.AddOns.FurrowPlane = cell([numZSlices,numTFrames]);
    handles.AddOns.backgroundROI = [];
    handles.AddOns.innerCircFit = cell([numTFrames,4]);
    
    handles.Variables.dt = str2num(ImageInfo{2}{9});
    handles.Variables.px = [str2num(ImageInfo{2}{6}) str2num(ImageInfo{2}{7}) str2num(ImageInfo{2}{8})] ;
    handles.Variables.BGFI = nan([numZSlices,numChannels,numTFrames]);
    handles.Variables.MeanFI = nan([numTFrames,numChannels,2]);
    handles.Variables.Bleaching= nan([numTFrames,numChannels]);
    % all spatial units = um!
    handles.Variables.Radius = nan([numTFrames,2]);
    handles.Variables.Volume = nan([numTFrames,2]);
    handles.Variables.SurfaceArea = nan([numTFrames,2]);
    handles.Variables.FurrowWidth = nan([numTFrames,1]);
    
    handles.Variables.Linescan.InnerCircLength = nan([numTFrames,2]);
    handles.Variables.Linescan.InnerCircPoints = cell([numTFrames,2]);
    handles.Variables.Linescan.LinescanIntensityMatrix = cell([numTFrames,numChannels,2]);
    
    handles.Parameters = struct([]);
    
    handles.current.currentImageStack = handles.File.ImageStack;
    handles.current.imageHeightWidth = [imageHeight,imageWidth];
    handles.current.numTFrames = numTFrames;
    handles.current.numChannels = numChannels;
    handles.current.numZSlices = numZSlices;
    handles.current.currentT = 1;
    for c=1:numChannels
        handles.current.maximumIntensity(c) = 1.2*max(max(ImageStack(:,:,c,round(numZSlices/2),1)));
    end
    handles.current.currentZoom = [0 0 imageWidth imageHeight];
    
elseif strcmp(FileName(length(FileName)-2:length(FileName)),'mat')
    % load images and variables from mat file
    vars = load(strcat(PathName,FileName));
    handles.File = vars.File;
    handles.File.pathName = PathName;
    handles.File.fileName = FileName;
    handles.Variables = vars.Variables;
    handles.Parameters = vars.Parameters;
    handles.AddOns = vars.AddOns;
   
    % and set some additional params
    handles.current.currentImageStack = handles.File.ImageStack;
    handles.current.imageHeightWidth = [str2num(handles.File.ImageInfo{2}{5}) str2num(handles.File.ImageInfo{2}{4})];
    handles.current.numChannels = str2num(handles.File.ImageInfo{2}{1});
    handles.current.numZSlices = str2num(handles.File.ImageInfo{2}{2});
    handles.current.numTFrames = str2num(handles.File.ImageInfo{2}{3});
    for c=1:handles.current.numChannels
        handles.current.maximumIntensity(c) = 1.2*max(max(handles.File.ImageStack(:,:,c,round(str2num(handles.File.ImageInfo{2}{2})/2),1)));
    end
    handles.current.currentZoom = [0 0 handles.current.imageHeightWidth(2) handles.current.imageHeightWidth(1)];
end

handles.current.currentC = 1;
handles.current.currentZ = 1;
handles.current.currentT = 1;
handles.current.showAddOns = 1;
set(handles.sliderT,'Min',1,'Max',handles.current.numTFrames,'SliderStep',[1/(handles.current.numTFrames-1) 5/(handles.current.numTFrames-1)],'Value',1);
set(handles.sliderZ,'Min',1,'Max',handles.current.numZSlices+0.01,'SliderStep',[0.01 0.1],'Value',1);

imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,handles.current.currentT),[0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(1),' / ',num2str(handles.current.numTFrames)));
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);




% --- Executes on slider movement.
function zSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.current.currentZ = round(get(hObject,'Value'));
set(hObject,'Value',handles.current.currentZ);
currentZoom = handles.current.currentZoom;
imshow(handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT),[0 handles.current.maximumIntensity(handles.current.currentC)]);
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function zSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.sliderZ = hObject;
guidata(hObject, handles);


% --- Executes on slider movement.
function tSlider_Callback(hObject, eventdata, handles)
% hObject    handle to tSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.current.currentT = round(get(hObject,'Value'));
set(hObject,'Value',handles.current.currentT);
currentZoom = handles.current.currentZoom;
imshow(handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT),[0 handles.current.maximumIntensity(handles.current.currentC)]);
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
handles.sliderT = hObject;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function showChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to showChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
guidata(hObject, handles);

% --- Executes when selected object is changed in showChannels.
function showChannels_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in showChannels
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'showCH1'
        handles.current.currentC=1;
        set(eventdata.NewValue,'String',strcat('CH1: ',handles.File.ImageInfo{2}{9+handles.current.currentC}));
    case 'showCH2'
        handles.current.currentC=2;
        if handles.current.currentC>handles.current.numChannels
            handles.current.currentC=handles.current.numChannels;
            set(eventdata.NewValue,'String',strcat('CH3: ',handles.File.ImageInfo{2}{9+handles.current.currentC}));
        else
            set(eventdata.NewValue,'String',strcat('CH2: ',handles.File.ImageInfo{2}{9+handles.current.currentC}));
        end
    case 'showCH3'
        handles.current.currentC=3;
        if handles.current.currentC>handles.current.numChannels
            handles.current.currentC=handles.current.numChannels;
            set(eventdata.NewValue,'String',strcat('CH3: ',handles.File.ImageInfo{2}{9+handles.current.currentC}));
        else
            set(eventdata.NewValue,'String',strcat('CH3: ',handles.File.ImageInfo{2}{9+handles.current.currentC}));
        end
end
currentZoom = handles.current.currentZoom;
imshow(handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT),[0 handles.current.maximumIntensity(handles.current.currentC)]);
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);









% --------------------------------------------------------------------
function view_Callback(hObject, eventdata, handles)
% hObject    handle to view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function zoomIn_Callback(hObject, eventdata, handles)
% hObject    handle to zoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = imrect;
position = round(getPosition(h)); % [xmin ymin width height]
delete(h)
% for multiple zooms, correct xmin and ymin
currentZoom = handles.current.currentZoom;
currentZoom = [position(1)+currentZoom(1), position(2)+currentZoom(2), position(3),position(4)];
% show zoomed image
imshow(handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT),...
    [0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
handles.current.currentZoom = position;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);

% --------------------------------------------------------------------
function zoomOut_Callback(hObject, eventdata, handles)
% hObject    handle to zoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current.currentZoom = [0 0 handles.current.imageHeightWidth(2) handles.current.imageHeightWidth(1)];
imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,handles.current.currentT),...
    [0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);

% --------------------------------------------------------------------
function show_Callback(hObject, eventdata, handles)
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function showAddOnsON_Callback(hObject, eventdata, handles)
% hObject    handle to showAddOnsON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.current.showAddOns
    handles.current.showAddOns = 0;
else
    handles.current.showAddOns = 1;
end
imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,handles.current.currentT),...
    [0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);







% --------------------------------------------------------------------
function processing_Callback(hObject, eventdata, handles)
% hObject    handle to processing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on slider movement.
function segmentationThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to segmentationThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.current.showAddOns = 1;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
% get threshold
handles.current.segmentationThreshold = get(hObject,'Value');
set(hObject,'Value',handles.current.segmentationThreshold);
% and current image
currentZoom = handles.current.currentZoom;
currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT);
% segment cell
BW = LSCAN_segmentCell(currentImage, handles.current.segmentationThreshold);
imshow(currentImage,[0 handles.current.maximumIntensity(handles.current.currentC)])
% store segmentation
handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT) = BW;
% plot segmentation
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function segmentationThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to segmentationThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', 0,'Max',1,'Value',0,'SliderStep',[0.01 0.1]);
guidata(hObject, handles);

% --------------------------------------------------------------------
function segmentCurrentCZT_Callback(hObject, eventdata, handles)
% hObject    handle to segmentCurrentCZT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current.showAddOns = 1;
currentZoom = handles.current.currentZoom;
currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT);
BW = LSCAN_segmentCell(currentImage, handles.current.segmentationThreshold);
% store segmentation
handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT) = BW;
% plot segmentation
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);

% --------------------------------------------------------------------
function segmentAll_Callback(hObject, eventdata, handles)
% hObject    handle to segmentAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current.showAddOns = 1;
currentZoom = handles.current.currentZoom;
for t = 1: handles.current.numTFrames
    handles.current.currentT = t;
    for z = handles.current.currentZ %1: handles.current.numZSlices
        handles.current.currentZ = z;
        currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,z,t);
        imshow(currentImage,[0 handles.current.maximumIntensity(handles.current.currentC)]);
        title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
        BW = LSCAN_segmentCell(currentImage, handles.current.segmentationThreshold);
        % store segmentation
        handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,z,t) = BW;
        % plot segmentation
        LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
        set(handles.sliderZ, 'Value', z);
        pause(0.1);
    end
    set(handles.sliderT, 'Value', t);
    guidata(hObject, handles);
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function correctSegmentationClear_Callback(hObject, eventdata, handles)
% hObject    handle to correctSegmentationClear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% remove segmentation for this frame:
handles.current.showAddOns = 1;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
currentZoom = handles.current.currentZoom;
handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT) = 0;
currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT);
imshow(currentImage, [0 handles.current.maximumIntensity(handles.current.currentC)]);
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);

% --------------------------------------------------------------------
function correctSegmentationEraser_Callback(hObject, eventdata, handles)
% hObject    handle to correctSegmentationEraser (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% draw rois which are then removed
handles.current.showAddOns = 1;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
currentZoom = handles.current.currentZoom;
currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT);
h_im = imshow(currentImage, [0 handles.current.maximumIntensity(handles.current.currentC)]);
% draw roi
h = imfreehand;
BW = createMask(h,h_im);
delete(h)
currentMask = handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT);
% delete region from mask
currentMask(BW==1)=0;
% store new mask
handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT) = currentMask;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);

% --------------------------------------------------------------------
function correctSegmentationSelector_Callback(hObject, eventdata, handles)
% hObject    handle to correctSegmentationSelector (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% select those objects which are real
handles.current.showAddOns = 1;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
currentZoom = handles.current.currentZoom;
currentMask = handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT);
currentMask = bwselect(currentMask,4);
% store new mask
handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT) = currentMask;
currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT);
imshow(currentImage, [0 handles.current.maximumIntensity(handles.current.currentC)]);
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);

% --------------------------------------------------------------------
function correctAddSubregion_Callback(hObject, eventdata, handles)
% hObject    handle to correctAddSubregion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current.showAddOns = 1;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
currentZoom = handles.current.currentZoom;
currentImage = handles.current.currentImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),handles.current.currentC,handles.current.currentZ,handles.current.currentT);
h_im = imshow(currentImage, [0 handles.current.maximumIntensity(handles.current.currentC)]);
% draw roi
h = imfreehand;
BW = createMask(h,h_im);
currentMask = handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT);
delete(h)
% delete region from mask
currentMask(BW==1)=1;
currentMask = imfill(currentMask,'holes');
% store new mask
handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT) = currentMask;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));
guidata(hObject, handles);










% --------------------------------------------------------------------
function correctIntensities_Callback(hObject, eventdata, handles)
% hObject    handle to correctIntensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function selectExternalROI_Callback(hObject, eventdata, handles)
% hObject    handle to selectExternalROI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentZoom = handles.current.currentZoom;
h = imrect;
position = getPosition(h); % N-by-2 array [X1 Y1;...;XN YN]
% correct for current zoom
position(1) = position(1)+currentZoom(1);
position(2) = position(2)+currentZoom(2);
backgPos = position;
handles.AddOns.backgroundROI = backgPos;
delete(h)

backgroundmask = poly2mask([backgPos(1) backgPos(1)+backgPos(3) backgPos(1)+backgPos(3) backgPos(1)],...
    [backgPos(2) backgPos(2) backgPos(2)+backgPos(4) backgPos(2)+backgPos(4)],...
    handles.current.imageHeightWidth(1),handles.current.imageHeightWidth(2));

% get mean intensity values in this ROI on every c,z,t
for t = 1:handles.current.numTFrames
    handles.current.currentT = t;
    set(handles.sliderT, 'Value', t);
    for c=1:handles.current.numChannels
        handles.current.currentC = c;
        for z=handles.current.currentZ %1:handles.current.numZSlices
            handles.current.currentZ = z;
            set(handles.sliderZ, 'Value', z);
            currentImage = handles.File.ImageStack(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),c,z,t);
            handles.Variables.BGFI(z,c,t) = mean(double(currentImage(backgroundmask==1)));
        end
    end
    guidata(hObject, handles);
end

guidata(hObject, handles);
















% --------------------------------------------------------------------
function getFurrowPlanes_Callback(hObject, eventdata, handles)
% hObject    handle to getFurrowPlanes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function furrowCurrentT_Callback(hObject, eventdata, handles)
% hObject    handle to furrowCurrentT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current.currentZoom = [0 0 handles.current.imageHeightWidth(2) handles.current.imageHeightWidth(1)];
z = handles.current.currentZ;
currentZoom = handles.current.currentZoom;
imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,handles.current.currentT),...
    [0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(handles.current.currentT),' / ',num2str(handles.current.numTFrames)));

% find center of cell object
cellmask = handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,handles.current.currentT);
cellmask(cellmask>0)=1;
% find center and major axis
stats = regionprops(cellmask,'Centroid','Orientation'); % Centroid = x,y
Centroid = [stats.Centroid];
Orientation = -[stats.Orientation]*2*pi/360;
% plot major axis
hold on
slope = tan(Orientation);
intercept = Centroid(2)-slope*Centroid(1);
plot([Centroid(1)-200,Centroid(1)+200],slope*[Centroid(1)-200,Centroid(1)+200]+intercept,'b')

% get distance between points on boundary of object and major axis line
% slope of the line: -a/b; intercept: - c/b ; points: x0, y0
% --> distance = abs(a/bx0+y0+c/b)/sqrt((a/b)^2+1)
B = bwboundaries(cellmask); % row,col = y,x
boundary = B{1};
boundary(:,1) = smooth(boundary(:,1),10);
boundary(:,2) = smooth(boundary(:,2),10);
distance = smooth(abs(-slope*boundary(:,2)+boundary(:,1)-intercept)./(sqrt(slope^2+1)),20);
[pks,locs] = findpeaks(-distance,'NPEAKS',4,'MINPEAKHEIGHT',-50);
peaks=[locs,-pks];

% remove smallest 2 (cell)
peaks = sortrows(peaks,2); peaks(1:2,:)=[];
plot(boundary(peaks(:,1),2),boundary(peaks(:,1),1),'b')
% store furrow center (x,y) and slope and length
handles.AddOns.FurrowPlane{handles.current.currentT}(z,:) = [mean(boundary(peaks(:,1),2)), mean(boundary(peaks(:,1),1)), diff(boundary(peaks(:,1),1))./diff(boundary(peaks(:,1),2)), sqrt((diff(boundary(peaks(:,1),1)))^2+(diff(boundary(peaks(:,1),2)))^2)];
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);

% --------------------------------------------------------------------
function furrowAllT_Callback(hObject, eventdata, handles)
% hObject    handle to furrowAllT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.current.currentZoom = [0 0 handles.current.imageHeightWidth(2) handles.current.imageHeightWidth(1)];
z = handles.current.currentZ;
currentZoom = handles.current.currentZoom;

for t = handles.current.currentT:handles.current.numTFrames
    set(handles.sliderT,'Value',t)
    handles.current.currentT=t;
    imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,t),...
        [0 handles.current.maximumIntensity(handles.current.currentC)]);
    title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
    %LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
    
    % find center of cell object
    cellmask = handles.AddOns.CellMask(currentZoom(2)+1:currentZoom(2)+currentZoom(4),currentZoom(1)+1:currentZoom(1)+currentZoom(3),1,handles.current.currentZ,t);
    cellmask(cellmask>0)=1;
    % find center and major axis
    stats = regionprops(cellmask,'Centroid','Orientation'); % Centroid = x,y
    Centroid = [stats.Centroid];
    Orientation = -[stats.Orientation]*2*pi/360;
    % plot major axis
    slope = tan(Orientation);
    intercept = Centroid(2)-slope*Centroid(1);
    
    % test
    % hold on
    % plot([Centroid(1)-100,Centroid(1)+100],slope*[Centroid(1)-100,Centroid(1)+100]+intercept,'b');
    
    % get distance between points on boundary of object and major axis line
    % slope of the line: -a/b; intercept: - c/b ; points: x0, y0
    % --> distance = abs(a/bx0+y0+c/b)/sqrt((a/b)^2+1)
    B = bwboundaries(cellmask); % row,col = y,x
    boundary = B{1};
    boundary(:,1) = smooth(boundary(:,1),10);
    boundary(:,2) = smooth(boundary(:,2),10);
    distance = smooth(abs(-slope*boundary(:,2)+boundary(:,1)-intercept)./(sqrt(slope^2+1)),20);
    [pks,locs] = findpeaks(-distance,'NPEAKS',4,'MINPEAKHEIGHT',-50);
    peaks=[locs,-pks];
    
    % test
    %     figure
    %     plot(-distance);
    %     hold on
    %     plot(locs,pks,'o')
    
    try
        % remove smallest 2 (cell)
        peaks = sortrows(peaks,2); peaks(1:2,:)=[];
        % plot furrow
        hold on
        plot(boundary(peaks(:,1),2),boundary(peaks(:,1),1),'b','LineWidth',3)
        % store furrow center (x,y) and slope
        handles.AddOns.FurrowPlane{t}(z,:) = [mean(boundary(peaks(:,1),2)), mean(boundary(peaks(:,1),1)), diff(boundary(peaks(:,1),1))./diff(boundary(peaks(:,1),2)), sqrt((diff(boundary(peaks(:,1),1)))^2+(diff(boundary(peaks(:,1),2)))^2)];
        LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
        guidata(hObject, handles);
    catch
    end
    pause(0.1);
    
end
guidata(hObject, handles);

% --------------------------------------------------------------------
function correctFurrowPosition_Callback(hObject, eventdata, handles)
% hObject    handle to correctFurrowPosition (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t = handles.current.currentT;
z = handles.current.currentZ;
imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,t),...
    [0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);

% draw line
h = imline;
position = wait(h); % [X1 Y1; X2 Y2].
delete(h);

% for t =1:handles.current.numTFrames

% furrow coordinates:
cellmask = handles.AddOns.CellMask(:,:,1,z,t);
B = bwboundaries(cellmask);
boundary = B{1};
[X,Y]=LSCAN_curveintersect(position(:,1),position(:,2),boundary(:,2),boundary(:,1));

handles.AddOns.FurrowPlane{t}(z,1:4) = [mean(X), mean(Y), diff(Y)./diff(X), sqrt((diff(Y))^2+(diff(X))^2)];

% end

% draw it
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
clear cellmask B boundary

guidata(hObject, handles);







% --- Executes on slider movement.
function cortexHeight_Callback(hObject, eventdata, handles)
% hObject    handle to cortexHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.AddOns.CortexHeight = get(hObject,'Value');
z = handles.current.currentZ;
t = handles.current.currentT;
imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,t),...
    [0 handles.current.maximumIntensity(handles.current.currentC)]);
title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
cellmask = handles.AddOns.CellMask(:,:,1,z,t);
cellmask(cellmask>0)=1;
outside = bwmorph(cellmask,'dilate',handles.AddOns.CortexHeight);
handles.AddOns.CellMask(:,:,1,z,t) = cellmask;
handles.AddOns.CellOutside(:,:,1,z,t) = outside;
LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function cortexHeight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cortexHeight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
set(hObject, 'Min', 1,'Max',100,'Value',1,'SliderStep',[0.01 1]);
guidata(hObject, handles);

% --------------------------------------------------------------------
function separatePoles_Callback(hObject, eventdata, handles)
% hObject    handle to separatePoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = waitbar(0,'Separating Poles');

% reset in case the function has been run before
CellMask = handles.AddOns.CellMask;
CellMask(CellMask>0)=1;
handles.AddOns.CellMask=CellMask;

% separate poles
for t = 1:handles.current.numTFrames
    set(handles.sliderT,'Value',t)
    handles.current.currentT=t;
    % use separation lines calculated from furrow plane
    furrowplane = handles.AddOns.FurrowPlane{t};
    meanlength = max(handles.current.imageHeightWidth);
    for z = handles.current.currentZ %1:handles.current.numZSlices
        set(handles.sliderZ,'Value',z)
        handles.current.currentZ=z;
        % apply cortex height addition in case it hasn't been done yet
        CellMask = handles.AddOns.CellMask(:,:,1,z,t);
        outside = bwmorph(CellMask,'dilate',handles.AddOns.CortexHeight);
        handles.AddOns.CellOutside(:,:,1,z,t) = outside;
        % on each plane, set segmented cell to the left ==1 and to the right ==2
        % first, get line coordinates at the edge of the image
        slope = furrowplane(z,3); intercept = furrowplane(z,2)-slope*furrowplane(z,1);
        centerx = furrowplane(z,1); centery = furrowplane(z,2);
        % if intercept is negative,
        if intercept<=0 || intercept>=handles.current.imageHeightWidth(1)
            x = [centerx-cos(atan(slope))*meanlength centerx+cos(atan(slope))*meanlength];
            y = [centery-sin(atan(slope))*meanlength centery+sin(atan(slope))*meanlength];
            yl = 1:1:handles.current.imageHeightWidth(1);
            xl = round(interp1(y,x,yl,'linear','extrap'));
            % make a mask with all pixels in one half == 1
            yl(length(yl)+1:length(yl)+2) = [handles.current.imageHeightWidth(1) 1];
            xl(length(xl)+1:length(xl)+2) = [handles.current.imageHeightWidth(2) handles.current.imageHeightWidth(2)];
            mask = poly2mask(xl, yl, handles.current.imageHeightWidth(1), handles.current.imageHeightWidth(2));
            % now set the pixel values in this region in the CellMask == 2
            currentMask = handles.AddOns.CellMask(:,:,1,z,t);
            currentMask(currentMask==1 & mask==1)=2;
            currentMaskO = handles.AddOns.CellOutside(:,:,1,z,t);
            currentMaskO(currentMaskO==1 & mask==1)=2;
            handles.AddOns.CellMask(:,:,1,z,t) = currentMask;
            handles.AddOns.CellOutside(:,:,1,z,t) = currentMaskO;
            handles.AddOns.ExternalMask(:,:,1,z,t) = mask;
        else
            % on each plane, set segmented cell to the left ==1 and to the right ==2
            % first, get line coordinates at the edge of the image
            x = [centerx-cos(atan(slope))*meanlength centerx+cos(atan(slope))*meanlength];
            y = [centery-sin(atan(slope))*meanlength centery+sin(atan(slope))*meanlength];
            xl = 1:1:handles.current.imageHeightWidth(2);
            yl = round(interp1(x,y,xl,'linear','extrap'));
            % make a mask with all pixels in one half == 1
            yl(length(yl)+1:length(yl)+2) = [handles.current.imageHeightWidth(1) handles.current.imageHeightWidth(1)];
            xl(length(xl)+1:length(xl)+2) = [handles.current.imageHeightWidth(2) 1];
            mask = poly2mask(xl, yl, handles.current.imageHeightWidth(1), handles.current.imageHeightWidth(2));
            if slope>0
                % now set the pixel values outside this region in the CellMask == 2
                currentMask = handles.AddOns.CellMask(:,:,1,z,t);
                currentMask(currentMask==1 & mask==0)=2;
                currentMaskO = handles.AddOns.CellOutside(:,:,1,z,t);
                currentMaskO(currentMaskO==1 & mask==0)=2;
                handles.AddOns.CellMask(:,:,1,z,t) = currentMask;
                handles.AddOns.CellOutside(:,:,1,z,t) = currentMaskO;
                handles.AddOns.ExternalMask(:,:,1,z,t) = mask;
            else
                % now set the pixel values inside this region in the CellMask == 2
                currentMask = handles.AddOns.CellMask(:,:,1,z,t);
                currentMask(currentMask==1 & mask==1)=2;
                currentMaskO = handles.AddOns.CellOutside(:,:,1,z,t);
                currentMaskO(currentMaskO==1 & mask==1)=2;
                handles.AddOns.CellMask(:,:,1,z,t) = currentMask;
                handles.AddOns.CellOutside(:,:,1,z,t) = currentMaskO;
                handles.AddOns.ExternalMask(:,:,1,z,t) = mask;
            end
        end
        
    end
    waitbar(t/handles.current.numTFrames,h);
    guidata(hObject, handles);
end
close(h)
guidata(hObject, handles);

% --------------------------------------------------------------------
function getGeometry_Callback(hObject, eventdata, handles)
% hObject    handle to getGeometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Variables.SectionArea = zeros(handles.current.numTFrames,1);

h = waitbar(0,'Getting Geometrical Estimates');
for t = 1:handles.current.numTFrames
    set(handles.sliderT,'Value',t)
    handles.current.currentT=t;
    z = handles.current.currentZ;
    furrowplane = handles.AddOns.FurrowPlane{t};
    CellMask = handles.AddOns.CellMask(:,:,1,z,t);
    
    % Jochen's approximation of sphere with spherical cap cut off
    handles.Variables.FurrowWidth(t,1) = furrowplane(z,4)*handles.Variables.px(1);
    sectionArea(1,1) = length(CellMask(CellMask==1))*handles.Variables.px(1)^2;
    sectionArea(1,2) = length(CellMask(CellMask==2))*handles.Variables.px(1)^2;
    handles.Variables.Radius(t,1) = LSCAN_Radius(handles.Variables.FurrowWidth(t,1)./2,sectionArea(1,1));
    handles.Variables.Radius(t,2) = LSCAN_Radius(handles.Variables.FurrowWidth(t,1)./2,sectionArea(1,2));
    Height1=handles.Variables.Radius(t,1).*(1+sqrt(1-(handles.Variables.FurrowWidth(t,1)./2).^2./(handles.Variables.Radius(t,1).^2)));
    Height2=handles.Variables.Radius(t,2).*(1+sqrt(1-(handles.Variables.FurrowWidth(t,1)./2).^2./(handles.Variables.Radius(t,2).^2)));
    handles.Variables.Volume(t,1) = pi/6*Height1.*(3*handles.Variables.FurrowWidth(t,1).^2/4+Height1.^2);
    handles.Variables.Volume(t,2) = pi/6*Height2.*(3*handles.Variables.FurrowWidth(t,1).^2/4+Height2.^2);
    handles.Variables.SurfaceArea(t,1) = 2*pi*(Height1.*handles.Variables.Radius(t,1));
    handles.Variables.SurfaceArea(t,2) = 2*pi*(Height2.*handles.Variables.Radius(t,2));
    
    %     % for dye equilibration
    %     handles.Variables.SectionArea(t) = sectionArea(1,1)+sectionArea(1,2);
    %
    waitbar(t/handles.current.numTFrames,h);
    guidata(hObject, handles);
end
close(h)



figure
plot(1:length(handles.Variables.SurfaceArea),handles.Variables.SurfaceArea)
title('SurfaceArea LHS and RHS vs timestep')


guidata(hObject, handles);









% --------------------------------------------------------------------
function tteModelAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to tteModelAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function measureIntensities_Callback(hObject, eventdata, handles)
% hObject    handle to measureIntensities (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Variables.MeanInside = zeros(size(handles.Variables.MeanFI));
h = waitbar(0,'Measuring Intensities');
for t = 1:handles.current.numTFrames
    set(handles.sliderT,'Value',t)
    handles.current.currentT=t;
    
    cellmask = double(handles.AddOns.CellMask(:,:,1,handles.current.currentZ,t));
    celloutsidemask = double(handles.AddOns.CellOutside(:,:,1,handles.current.currentZ,t));
    fimask = celloutsidemask - cellmask;
    
    for c = 1:handles.current.numChannels
        currentImages = double(handles.File.ImageStack(:,:,c,handles.current.currentZ,t));
        % pole 1
        handles.Variables.MeanFI(t,c,1)=nanmean(currentImages(fimask==1))-handles.Variables.BGFI(handles.current.currentZ,c,t);
        % pole 2
        handles.Variables.MeanFI(t,c,2)=nanmean(currentImages(fimask==2))-handles.Variables.BGFI(handles.current.currentZ,c,t);
        % for dye equilibration
        %         % mean inside
        %         handles.Variables.MeanInside(t,c)=nanmean(currentImages(cellmask>0));
        %         handles.Variables.MeanFI(t,c)=nanmean(currentImages(fimask>0))-handles.Variables.BGFI(handles.current.currentZ,c,t);
    end
    handles.Variables.MeanFI(handles.Variables.MeanFI==0)=nan;
    handles.Variables.MeanFI(handles.Variables.MeanFI==0)=nan;
    waitbar(t/handles.current.numTFrames,h);
    guidata(hObject, handles);
end
close(h)
figure
% subplot(1,2,1)
plot(handles.Variables.dt*(0:length(handles.Variables.MeanFI(:,:,1))-1),handles.Variables.MeanFI(:,:,1))
% hold on
% plot(handles.Variables.dt*(0:length(handles.Variables.MeanInside(:,:,1))-1),handles.Variables.MeanInside(:,:,1))
% title('CAAX and DYE PM in round cells (not bleaching-corrected)','FontSize',18);
% xlabel('time (s)','FontSize',18)
% ylabel('mean intensity','FontSize',18);
% legend('DYE','CAAX');
% subplot(1,2,2)
% plot(handles.Variables.dt*(0:length(handles.Variables.MeanFI(:,:,1))-1),handles.Variables.MeanInside)
% title('CAAX and DYE CYTO in round cells (not bleaching-corrected)','FontSize',18);
% xlabel('time (s)','FontSize',18)
% ylabel('mean intensity','FontSize',18);
% legend('DYE','CAAX');


guidata(hObject, handles);


% --------------------------------------------------------------------
function fitTTEModel_Callback(hObject, eventdata, handles)
% hObject    handle to fitTTEModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dt = handles.Variables.dt;
pix = handles.Variables.px(1); %um per pixel

% find lifeact channel
lact = 0;
for c = 1:handles.current.numChannels
    if strcmp(handles.File.ImageInfo{2}{9+c},'LACT')
        lact = c;
    end
end

if lact>0
    % rewrite some variables for fit function
    time = dt*(1:length(handles.Variables.MeanFI(:,lact,1)));
    LHSactinCortexMeanNorm = handles.Variables.MeanFI(:,lact,1)./((handles.Variables.MeanFI(:,lact,1)+handles.Variables.MeanFI(:,lact,2))/2);
    RHSactinCortexMeanNorm = handles.Variables.MeanFI(:,lact,2)./((handles.Variables.MeanFI(:,lact,1)+handles.Variables.MeanFI(:,lact,2))/2);
    geometry.Radius = handles.Variables.Radius;
    geometry.Volume = handles.Variables.Volume;
    geometry.SurfaceArea = handles.Variables.SurfaceArea;
    geometry.FurrowWidth = handles.Variables.FurrowWidth;
    
    % get automatic start and end values based on volume
    smoothy = 20;
    v = smooth((handles.Variables.Volume(:,1)-handles.Variables.Volume(:,2))./(handles.Variables.Volume(:,1)+handles.Variables.Volume(:,2)),smoothy);
    % find volume extrema of large oscillation
    [vmaxima,vmaxlocs] = findpeaks(v,'MINPEAKHEIGHT',0.25); %0.35
    [vminima,vminlocs] = findpeaks(-v,'MINPEAKHEIGHT',0.25); % 0.35
    vpeaks(:,1)=[vmaxlocs; vminlocs];
    vpeaks(:,2)=[vmaxima; -vminima];
    vpeaks = sortrows(vpeaks,1);
    % exclude everything one quarter period before the first extremum and after the
    % last extremum
    exbefore = vpeaks(1,1)-round((vpeaks(2,1)-vpeaks(1,1))/2); exbefore(exbefore<1)=1;
    exafter = vpeaks(end,1)+round((vpeaks(2,1)-vpeaks(1,1))/2); exafter(exafter>length(v))=length(v); %round((vpeaks(2,1)-vpeaks(1,1))/2)
    included = logical(ones(length(v),1));
    included(1:exbefore,1)=0; included(exafter:end,1)=0;
    excluded = ~included;
    
    
    % figure to show analysed data
    h = figure;
    file_title = 'mechanical_analysis_profiles';
    
    plot(time,LHSactinCortexMeanNorm,'r',time,RHSactinCortexMeanNorm,'m',time,1+v,'k','LineWidth',2)
    legend('Actin 1','Actin 2','v+1');
    hold on
    area(time,2*excluded,'FaceColor',[.9 .9 .9]); alpha(0.5);
    xlim([min(time) max(time)]);ylim([0 2]);
    title('Mean Actin Profiles', 'FontWeight', 'bold','FontSize',16);
    xlabel('time (s)','FontWeight','bold','FontSize',14);
    
    saveas(h,[handles.File.pathName,'figures/',handles.File.fileName(1:end-4),'_',file_title,'.fig'])
    saveas(h,[handles.File.pathName,'figures/',handles.File.fileName(1:end-4),'_',file_title,'.eps'],'psc2')
    
    % fits to TTE model
    i = find(included>0);
    from = i(1); to = i(end);
    handles.Parameters(1).TTE = LSCAN_mech_analyseOsciCells2(handles.File,LHSactinCortexMeanNorm,RHSactinCortexMeanNorm,geometry,from,to,pix,dt);
end
guidata(hObject,handles);














% --------------------------------------------------------------------
function linescanAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to linescanAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function fitCellOutline_Callback(hObject, eventdata, handles)
% hObject    handle to fitCellOutline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for t = 1:handles.current.numTFrames
    handles.current.currentT = t;
    set(handles.sliderT,'Value',t);
    
    % POLE 1
    % 1. fit inner cell circumference in polar coordinates to smoothing spline
    cellmask = handles.AddOns.CellMask(:,:,1,handles.current.currentZ,t);
    cellmask(cellmask==2)=0; % only pole 1
    stats = regionprops(cellmask,'Centroid'); % x, y
    centroid1 = [stats(1).Centroid];
    innerCellCirc = bwboundaries(cellmask); % y, x
    innerCellCircX1 = innerCellCirc{1}(:,2)-centroid1(1);
    innerCellCircY1 = innerCellCirc{1}(:,1)-centroid1(2);
    % get polar coordinates
    [THETA,RHO] = cart2pol(innerCellCircX1,innerCellCircY1);
    [innerCircFit1, gof] = LSCAN_innerCircSmoothingSplineFit([THETA-2*pi THETA THETA+2*pi THETA+4*pi],[RHO RHO RHO RHO]);
    handles.AddOns.innerCircFit{t,1} = innerCircFit1; % theta, rho
    handles.AddOns.innerCircFit{t,3} = centroid1; % x,y
    
    % POLE 2
    % 1. fit inner cell circumference in polar coordinates to smoothing spline
    cellmask = handles.AddOns.CellMask(:,:,1,handles.current.currentZ,t);
    cellmask(cellmask==1)=0; cellmask(cellmask==2)=1; % only pole 2
    stats = regionprops(cellmask,'Centroid'); % x, y
    centroid2 = [stats(1).Centroid];
    innerCellCirc = bwboundaries(cellmask); % y, x
    innerCellCircX2 = innerCellCirc{1}(:,2)-centroid2(1);
    innerCellCircY2 = innerCellCirc{1}(:,1)-centroid2(2);
    % get polar coordinates
    [THETA,RHO] = cart2pol(innerCellCircX2,innerCellCircY2);
    [innerCircFit2, gof] = LSCAN_innerCircSmoothingSplineFit([THETA-2*pi THETA THETA+2*pi THETA+4*pi ],[RHO RHO RHO RHO]);
    handles.AddOns.innerCircFit{t,2} = innerCircFit2; % theta, rho
    handles.AddOns.innerCircFit{t,4} = centroid2; % x,y
    guidata(hObject,handles);
    
    % in prep for linescan
    % furrow centroid
    centroidF = handles.AddOns.FurrowPlane{t}(handles.current.currentZ,:);
    
    % Pole 1
    theta = [0:0.01*pi:2*pi,0];
    [innerFit1X, innerFit1Y] = pol2cart(theta(:),innerCircFit1(theta(:)));
    % Pole 2
    [innerFit2X, innerFit2Y] = pol2cart(theta(:),innerCircFit2(theta(:)));
    
    % use furrow points as starting and end-ANGLES
    furrowX = [centroidF(1)-cos(atan(centroidF(3)))*centroidF(4)/2, centroidF(1)+cos(atan(centroidF(3)))*centroidF(4)/2];
    furrowY = [centroidF(2)-sin(atan(centroidF(3)))*centroidF(4)/2, centroidF(2)+sin(atan(centroidF(3)))*centroidF(4)/2];
    
    % 3. get positions on outline at a fixed distances d from intersection
    % point to intersection point, don't forget to translate back with
    % centroid later!
    % [xlinescan ylinescan] = getLinescanPoints(intersection,fitresult, from,to, distance)
    distance = 0.5; %pixels
    % Pole 1:
    [LinescanX1, LinescanY1] = LSCAN_getLinescanPoints(furrowX-centroid1(1),furrowY-centroid1(2),innerCircFit1,1,2,distance);
    % Pole 2:
    [LinescanX2, LinescanY2] = LSCAN_getLinescanPoints(furrowX-centroid2(1),furrowY-centroid2(2),innerCircFit2,2,1,distance);
    
    % 6. store results
    handles.Variables.Linescan.InnerCircLength(t,1) = length(LinescanX1)*distance;
    handles.Variables.Linescan.InnerCircLength(t,2) = length(LinescanX2)*distance;
    handles.Variables.Linescan.InnerCircPoints{t,1} = [LinescanX1(:)+centroid1(1), LinescanY1(:)+centroid1(2)];
    handles.Variables.Linescan.InnerCircPoints{t,2} = [LinescanX2(:)+centroid2(1), LinescanY2(:)+centroid2(2)];
    
    % SHOW IT
    imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,t),...
        [0 handles.current.maximumIntensity(handles.current.currentC)]);
    title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
    LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
    
    % test:
    % hold on
    % plot(LinescanX1+centroid1(1),LinescanY1+centroid1(2),'go',LinescanX2+centroid2(1),LinescanY2+centroid2(2),'ro')
    
    pause(0.1);
end
guidata(hObject,handles);

% --------------------------------------------------------------------
function getLinescans_Callback(hObject, eventdata, handles)
% hObject    handle to getLinescans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h1 = waitbar(0,'Getting linescans...');


for t=1:handles.current.numTFrames
    handles.current.currentT = t;
    set(handles.sliderT,'Value',t);
    
    % 4. get normal vectors of 40 pixel length from those positions outward
    % [ vectorX vectorY ] = getNormalVectorsForLinescan(fitresult,xPoints,yPoints,L)
    L = handles.AddOns.CortexHeight*15;
    % Pole 1:
    LinescanX1 = handles.Variables.Linescan.InnerCircPoints{t,1}(:,1)-handles.AddOns.innerCircFit{t,3}(1);
    LinescanY1 = handles.Variables.Linescan.InnerCircPoints{t,1}(:,2)-handles.AddOns.innerCircFit{t,3}(2);
    innerFit1 = handles.AddOns.innerCircFit{t,1};
    [ NormvectorX1, NormvectorY1 ] = LSCAN_getNormalVectorsForLinescan(innerFit1,LinescanX1,LinescanY1,L);
    % Pole 2:
    LinescanX2 = handles.Variables.Linescan.InnerCircPoints{t,2}(:,1)-handles.AddOns.innerCircFit{t,4}(1);
    LinescanY2 = handles.Variables.Linescan.InnerCircPoints{t,2}(:,2)-handles.AddOns.innerCircFit{t,4}(2);
    innerFit2 = handles.AddOns.innerCircFit{t,2};
    [ NormvectorX2, NormvectorY2 ] = LSCAN_getNormalVectorsForLinescan(innerFit2,LinescanX2,LinescanY2,L);
    
    
    % 5. translate everything back by centroid
    NormvectorX1 = NormvectorX1 + handles.AddOns.innerCircFit{t,3}(1);
    NormvectorX2 = NormvectorX2 + handles.AddOns.innerCircFit{t,4}(1);
    NormvectorY1 = NormvectorY1 + handles.AddOns.innerCircFit{t,3}(2);
    NormvectorY2 = NormvectorY2 + handles.AddOns.innerCircFit{t,4}(2);
    
    % do all channels
    for channel = 1:handles.current.numChannels
        % 7. do linescan along vector points
        % c = improfile(image,xEndPoints,yEndPoints,numberOfPointsInBetween)
        numPoints = round(L*2); % approximately two per pixel
        handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,1}=nan(numPoints,length(LinescanX1));
        handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,2}=nan(numPoints,length(LinescanX2));
        
        currentImage = handles.File.ImageStack(:,:,channel,handles.current.currentZ,t);
        % Pole 1
        for n=1:length(LinescanX1)
            handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,1}(:,n) = improfile(currentImage,NormvectorX1(:,n),NormvectorY1(:,n),numPoints);
        end
        % Pole 2
        for n=1:length(LinescanX2)
            handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,2}(:,n) = improfile(currentImage,NormvectorX2(:,n),NormvectorY2(:,n),numPoints);
        end
        
        % test:
        % figure, imshow(handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,1},[0 180]);
        % figure, imshow(handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,2},[0 180]);
    end
    
    % update waitbar
    waitbar(t/handles.current.numTFrames,h1);
    %
    %                 imshow(handles.current.currentImageStack(:,:,handles.current.currentC,handles.current.currentZ,t),...
    %                 [0 handles.current.maximumIntensity(handles.current.currentC)]);
    %                 title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
    %                 LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
    %                 % test:
    %                 hold on
    %                 plot(NormvectorX1,NormvectorY1,'g',NormvectorX2,NormvectorY2,'r')
    %                 % waitforbuttonpress;
    %                 pause(0.1);
    guidata(hObject,handles);
    
end
close(h1);
guidata(hObject,handles);



% --------------------------------------------------------------------
function makeLinescanMovies_Callback(hObject, eventdata, handles)
% hObject    handle to makeLinescanMovies (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Movie1(1:(handles.current.numTFrames)*handles.current.numZSlices) = struct('cdata', [],'colormap', []);
Movie2(1:(handles.current.numTFrames)*handles.current.numZSlices) = struct('cdata', [],'colormap', []);

% Pole 1
numCircPoints1 = zeros(1,handles.current.numTFrames);
numCircPoints2 = zeros(1,handles.current.numTFrames);
for t = 1:handles.current.numTFrames
    numCircPoints1(t)=size(handles.Variables.Linescan.LinescanIntensityMatrix{t,1,1},2);
    numCircPoints2(t)=size(handles.Variables.Linescan.LinescanIntensityMatrix{t,1,2},2);
end
maximumCircLength = max([nanmax(numCircPoints1), nanmax(numCircPoints2)]);
if rem(maximumCircLength,2)~=0
    maximumCircLength = maximumCircLength+1;
end
handles.Variables.Linescan.LinescanIntensityMatrixAligned = zeros([size(handles.Variables.Linescan.LinescanIntensityMatrix{1,1,1},1),maximumCircLength,handles.current.numChannels,handles.current.numZSlices,handles.current.numTFrames,2]);

for channel = 1:handles.current.numChannels
    h = figure;
    for t = 1:handles.current.numTFrames
        % POLE 1
        % 1. align linescan rows (peak position)
        M1 = handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,1};
        
        % 2. align linescan cols (center position)
        Mframes1 = zeros(size(M1,1),maximumCircLength);
        Mframes1(:,maximumCircLength/2-numCircPoints1(t)/2+1:maximumCircLength/2+numCircPoints1(t)/2)=flipud(M1);
        handles.Variables.Linescan.LinescanIntensityMatrixAligned(:,:,channel,handles.current.currentZ,t,1)=Mframes1; %{t,channel,1}=Mframes1; %(:,:,channel,handles.current.currentZ,t,1)
        
        % 3. make movies
        imshow(Mframes1,[0 handles.current.maximumIntensity(channel)]);
        title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
        Movie1(t) = getframe(h);
        
        pause(0.2);
        
    end
    close(h);
    h=figure;
    for t = 1:handles.current.numTFrames
        % POLE 2
        % 1. align linescan rows (peak position)
        M2 = handles.Variables.Linescan.LinescanIntensityMatrix{t,channel,2};
        
        % 2. align linescan cols (center position)
        Mframes2 = zeros(size(M2,1),maximumCircLength);
        Mframes2(:,maximumCircLength/2-numCircPoints2(t)/2+1:maximumCircLength/2+numCircPoints2(t)/2)=flipud(M2);
        handles.Variables.Linescan.LinescanIntensityMatrixAligned(:,:,channel,handles.current.currentZ,t,2)=Mframes2; %{t,channel,1}=Mframes1; %(:,:,channel,handles.current.currentZ,t,1)
        
        % 3. make movies
        imshow(Mframes2,[0 handles.current.maximumIntensity(channel)]);
        title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)));
        Movie2(t) = getframe(h);
        
        pause(0.2);
    end
    close(h);
    
    movie2avi(Movie1, strcat(handles.File.pathName,'figures/',handles.File.fileName(1:end-4),'_linescan_',handles.File.ImageInfo{2}{9+channel},'_pole1.avi'), 'compression', 'None');
    movie2avi(Movie2, strcat(handles.File.pathName,'figures/',handles.File.fileName(1:end-4),'_linescan_',handles.File.ImageInfo{2}{9+channel},'_pole2.avi'), 'compression', 'None');
    
end

guidata(hObject,handles);



% --------------------------------------------------------------------
function getLinescanParameters_Callback(hObject, eventdata, handles)
% hObject    handle to getLinescanParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Linescan = handles.Variables.Linescan;
Movie(1:(handles.current.numTFrames)*handles.current.numZSlices) = struct('cdata', [],'colormap', []);

% find caax channel
caax = 0;
for c = 1:handles.current.numChannels
    if strcmp(handles.File.ImageInfo{2}{9+c},'CAAX')
        caax = c;
    end
end

h = figure;
Linescan.MaximumIntensity = nan(handles.current.numTFrames,2,2);
Linescan.Width = nan(handles.current.numTFrames,2);
Linescan.IntegratedIntensity = nan(handles.current.numTFrames,2);
for t = 1:handles.current.numTFrames
    
    % Pole 1
    li = Linescan.LinescanIntensityMatrix{t,caax,1};
    % take only the central 50% for fitting and analysis
    s = size(li); startpl = round(s(2)/4); endpl = s(2)-startpl;
    % fit mean linescan-background
    
    [fitresult, gof] = LSCAN_linescanMeanFit([1:s(1)], nanmean(li(1:s(1),startpl:endpl),2)-handles.Variables.BGFI(handles.current.currentZ,caax,t));
    % get maximum and width
    [Linescan.MaximumIntensity(t,1,1), Linescan.MaximumIntensity(t,1,2)]=max(fitresult(1:round(s(1)/2)));
    % width: position where intensity has decreased to 0.37 of maximum (by 0.63 of maximum)
    [y, Linescan.Width(t,1)]=min(abs(((fitresult(Linescan.MaximumIntensity(t,1,2):s(1)))/Linescan.MaximumIntensity(t,1,1))-0.37));
    Linescan.IntegratedIntensity(t,1)= sum((fitresult(0:0.01:100))*0.01);
    
    % figure
    plot(((0:s(1)-1)-Linescan.MaximumIntensity(t,1,2))*handles.Variables.px(1)/2, (nanmean(li(1:s(1),startpl:endpl),2)-handles.Variables.BGFI(handles.current.currentZ,caax,t)),'r.',...
        ((0:s(1)-1)-Linescan.MaximumIntensity(t,1,2))*handles.Variables.px(1)/2,(fitresult(1:s(1))),'k');
    axis([-20*handles.Variables.px(1)/2 100*handles.Variables.px(1)/2 0 300]);
    xlabel('Position (um)','FontSize',18);
    ylabel('Intensity (a.u.)','FontSize',18)
    title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)),'FontSize',18);
    
    
    % Pole 2
    li = Linescan.LinescanIntensityMatrix{t,caax,2};
    % take only the central 50% for fitting and analysis
    s = size(li); startpl = round(s(2)/4); endpl = s(2)-startpl;
    % fit mean linescan-background
    [fitresult, gof] = LSCAN_linescanMeanFit([1:s(1)], nanmean(li(1:s(1),startpl:endpl),2)-handles.Variables.BGFI(handles.current.currentZ,caax,t));
    % get maximum and width
    [Linescan.MaximumIntensity(t,2,1), Linescan.MaximumIntensity(t,2,2)]=max(fitresult(1:round(s(1)/2)));
    % width: position where intensity has decreased to 0.37 of maximum (by 0.63 of maximum)
    [y, Linescan.Width(t,2)]=min(abs(((fitresult(Linescan.MaximumIntensity(t,2,2):s(1)))/Linescan.MaximumIntensity(t,2,1))-0.37));
    Linescan.IntegratedIntensity(t,2)= sum((fitresult(0:0.01:100))*0.01);
    
    % figure
    hold on
    plot(((0:s(1)-1)-Linescan.MaximumIntensity(t,2,2))*handles.Variables.px(1)/2, (nanmean(li(1:s(1),startpl:endpl),2)-handles.Variables.BGFI(handles.current.currentZ,caax,t)),'m.',...
        ((0:s(1)-1)-Linescan.MaximumIntensity(t,2,2))*handles.Variables.px(1)/2,(fitresult(1:s(1))),'k');
    Movie(t) = getframe(h);
    
    hold off
    pause(0.1);
end
close(h)

movie2avi(Movie, strcat(handles.File.pathName,'figures/',handles.File.fileName(1:end-4),'_caax_linescanProfiles.avi'), 'compression', 'None');

handles.Variables.Linescan = Linescan;
guidata(hObject,handles);

%
% %%%%%%%%%%%%%%%%%%%%%%%
% %%% PEAKS and RATES %%%
% %%%%%%%%%%%%%%%%%%%%%%%


[Parameters] = LSCAN_getMembraneParameters(handles.Variables,strcat(handles.File.pathName,'figures/',handles.File.fileName(1:end-4)));
handles.Parameters(1).Linescan = Parameters.Linescan;
guidata(hObject,handles);












% --------------------------------------------------------------------
function results_Callback(hObject, eventdata, handles)
% hObject    handle to results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function makeMovieAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to makeMovieAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear mov
mov(1:(handles.current.numTFrames)*handles.current.numZSlices) = struct('cdata', [],'colormap', []);% ZSlices
handles.current.showAddons = 1;

for t=1:handles.current.numTFrames
    set(handles.sliderT,'Value',t)
    handles.current.currentT = t;
    for z = handles.current.currentZ %1:handles.current.numZSlices
        set(handles.sliderZ,'Value',z)
        handles.current.currentZ = z;
        
        imshow(handles.current.currentImageStack(:,:,handles.current.currentC,z,t), [0 handles.current.maximumIntensity(handles.current.currentC)]);
        LSCAN_drawAddOns(handles.current.showAddOns,handles.AddOns,handles.current);
        title(strcat(num2str(t),' / ',num2str(handles.current.numTFrames)),'FontSize',18);
        
        mov(handles.current.numZSlices*(t-1)+z) = getframe(gcf);
        pause(0.1);
    end
end
movie2avi(mov, strcat(handles.File.pathName,'figures/',handles.File.fileName(1:end-4),'_movieAnalysis.avi'), 'compression', 'None');
guidata(hObject, handles);


% --------------------------------------------------------------------
function concatenateResults_Callback(hObject, eventdata, handles)
% hObject    handle to concatenateResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% select different folders
% for all folders, load all mat files and concatenate flow / folding data
% then boxplot
% first get the pathname of the folder which contains all the subfolders with the data
dirname = uigetdir();
LSCAN_concatenateResults_LACT( dirname );
%LSCAN_concatenateResults_CAAX( dirname );


% --------------------------------------------------------------------
function saveData_Callback(hObject, eventdata, handles)
% hObject    handle to saveData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
File = handles.File;
AddOns = handles.AddOns;
Variables = handles.Variables;
Parameters = handles.Parameters;
save(strcat(handles.File.pathName,handles.File.fileName(1:end-4),'.mat'), 'File','AddOns','Variables','Parameters','-v7.3');

% for dye equilibration
% Variables=rmfield(Variables,'Bleaching');
% Variables=rmfield(Variables,'Radius');
% Variables=rmfield(Variables,'Volume');
% Variables=rmfield(Variables,'SurfaceArea');
% Variables=rmfield(Variables,'Linescan');
% Variables=rmfield(Variables,'FurrowWidth');
% Variables.MeanFI=Variables.MeanFI(:,1,1);
% Variables.MeanInside=Variables.MeanInside(:,1,1);
% save(strcat(handles.File.pathName,handles.File.fileName(1:end-4),'_Variables.mat'), 'Variables','-v7.3');
