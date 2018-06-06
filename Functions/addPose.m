% function to compute a pose composition (add)
function poseSum = addPose(x1,x2)
% rotation matrix
if size(x2,2) > size(x1,1)
    x1 = repmat(x1,1,size(x2,2));
end

if size(x2,1) == 3

    s = sin(x1(3));
    c = cos(x1(3)); 
    rot = [c -s 0; 
           s  c 0;
           0  0 1];

    % Sum
    poseSum = x1 + rot*x2;

elseif size(x2,1) == 2
    s = sin(x1(3));
    c = cos(x1(3)); 
    rot = [c -s; 
           s  c];

    % Sum
    poseSum = x1(1:2,:) + rot*x2;
        
else
    error('invalid size');
end