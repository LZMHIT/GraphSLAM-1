% function to reverse a pose composition (subtract)
function poseDiff = subtractPose(x1,x2)
% rotation matrix
s = sin(x1(3));
c = cos(x1(3)); 
invRot = [ c  s 0; 
          -s  c 0;
           0  0 1];
   
% Subtract
poseDiff = invRot*(x1 - x2);