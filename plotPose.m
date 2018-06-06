function f = plotPose(x,y,theta,f)
if nargin < 4
    f = figure;
end
figure(f);

% Determine orientation
dx = cos(theta);
dy = sin(theta);
x_begin = x - dx;
y_begin = y - dy;

% Plot arrow with tip at x,y oriented by theta
headWidth = 10;
headLength = 15;
ah = annotation('arrow','LineStyle','none','headStyle','cback1',... % style
    'HeadLength',headLength,'HeadWidth',headWidth,... % head size
    'x',[x_begin x],'y',[y_begin y]); % position and orientation
set(ah,'parent',gca);