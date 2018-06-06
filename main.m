% Main function for graph slam

function main

% Add all relevant folders to path
addpath('./')
% addpath('./Test')
% addpath('./Data')
addpath('./Functions')
addpath('./Classes')

dt = 0.001;
tMax = 20;
t = 0:dt:tMax;
Q = 10*dt^2*eye(3);
R = 0.1*eye(2);
R_big = 0.1*eye(20);
P = 0.95; % error ellipse confidence level

numRobots = 1;
anchors = [10;10; pi/5];
reorderRate = 100;
mu0 = anchors;
sig0 = 0.1*eye(3);

controlTraj = [3*ones(1,length(t)-1); 0.2*sin(t(1:(end-1)))];
f = @(x,u) x + dt*[u(1)*cos(x(3)); u(1)*sin(x(3)); u(2)]; 
g = @(x,m,I) getMeasurements(x,m,I);
        
m = getMap('lofts_simple20.jpg');
        
% Initial Conditions
x0  = anchors;

% initialize state and measurement arrays
x = zeros(3,length(t));
x(:,1) = x0; % initial state
y = zeros(20,length(t)-1);
I = zeros(10,length(t)-1); % which features are measured each time

% create trajectories
for i = 1:(length(t)-1)
    % update states and take measurements
    processNoise = sqrtm(Q)*randn(3,1); % convert from std normal
    measureNoise = sqrtm(R_big)*randn(20,1); % convert from std normal
    x(:,i+1) = f(x(:,i),controlTraj(:,i)) + processNoise; % motion model
    I(:,i) = nearestFeatures(x(:,i+1),m); % determine which features are measured
%     I_all = getDetailedI(I(:,i));
    y(:,i) = g(x(:,i+1),m,I(:,i)) + measureNoise; % range measurement model
end

hold on;
plot(x(1,:),x(2,:),'k-')

measurements = y;
measCorrespond = I;
robot = Robot(numRobots,anchors,dt,Q,R,reorderRate,mu0,sig0,controlTraj,measurements,measCorrespond);


fig = figure(gcf);
for i = 1:(length(t)-1)
    robot.runIteration();
    plotRobot(robot,fig);
end
    
%% Function to determine which features are measured
function I = nearestFeatures(x,m)
numMapFeatures = size(m,2);
x_matrix = repmat(x(1:2),1,numMapFeatures);
range = sqrt(sum((m - x_matrix).^2,1))'; % transpose for col vector
[~,sortedIdx] = sort(range,1,'ascend');
I = sortedIdx(1:10,1);

%% function to get measurements
function g = getMeasurements(x,m_matrix,I)
numFeatures = length(I);
numMeas = 2*numFeatures;
g = zeros(numMeas,1);
for i = 1:numFeatures
    rangeRow = 2*i-1;
    bearingRow = 2*i;
    m = m_matrix(:,I(i));
    
    g([rangeRow,bearingRow],1) = [norm(m - x(1:2),2); % range
            (atan2(m(2) - x(2),m(1) - x(1)) - x(3))]; % bearing 
end

