% Source code for paper BUILDING DETECTION USING PROCESSING OF MONOCHROMATIC EARTH OBSERVATION IMAGE
% https://ieeexplore.ieee.org/document/9252585
% https://www.dl.begellhouse.com/journals/0632a9d54950b268,624c7ae26c910972,12beaa1f7e24a83c.html#
%
% Briefly - it finds buildings on space infrared images basing on the shape
% pattern. Assumed any building consist of stright edges and right angles.
% Alex.N.Gorobets@gmail.com

function [] = Demo()
%Just unpack image (Train1/5.tiff) to neighboring folder
% and run This script

%% Predefined consts
IMAGE_IN_FILE_NAME = 'Train1/5.tif';
MIN_BRIGHTNESS = 0;
MAX_BRIGHTNESS = 2^16-1;

% Process mode
HstgrmStretchDepth=0.95;

%% Prepocess
% Load Image
imageIn = imread(IMAGE_IN_FILE_NAME);

% Have effect looks like shifted histogramm - fix it.
fixedImageIn = HstgrmAutoRotate(imageIn, MIN_BRIGHTNESS, MAX_BRIGHTNESS);
% Do Autocontrast
readyImageIn = HstgrmAutoStretch(fixedImageIn, HstgrmStretchDepth, MIN_BRIGHTNESS, MAX_BRIGHTNESS);
% Crop image to save time.
readyImageIn = readyImageIn(:, end*7/8:end);

%% Morphologic process
% Get edges
strightEdges=FindEdgeLines(readyImageIn);
% Mark buildings
FindBldObj(strightEdges,readyImageIn,[]);

end

