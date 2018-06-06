%% function to get a map to explore
function m = getMap(mapName)
%% Pull in image
RGB = imread(mapName);
I = rgb2gray(RGB);
BW = imbinarize(I); % convert to logical
BW = ~BW; % invert so 1's correspond to black pixels

%% Find desired formation position and shift to be centered at origin
[Y,X] = meshgrid(1:size(BW,1),1:size(BW,2)); % create grid of positions
X = X'; Y = Y'; X = flipud(X); Y = flipud(Y); % reorient for coordinates to align
Xeq = X(BW); % take only black pixels
Yeq = Y(BW);
% Xeq = Xeq - mean(Xeq); % shift average to the origin
% Yeq = Yeq - mean(Yeq);
m = [Xeq,Yeq]'; % concatenate and transpose so each xeq_i is a column in m
figure
plot(m(1,:),m(2,:),'k.')
axis square
axis equal
