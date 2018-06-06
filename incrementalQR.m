% Function to add rows to R from QR and incrementally make it upper
% triangular again, while keeping RHS vector d up to date
function [R,d] = incrementalQR(R,d,newRowsR,newRowsD)
curMeas = size(R,1);
curVars = size(R,2);
newVars = size(newRowsR,2);

% If variables added, resize R
if newVars > curVars
    R(curMeas,newVars) = 0;
end

% Vertically concatenate new rows
R = [R;newRowsR];
d = [d;newRowsD];

% Iteratively zero out nonzero entries below diagonal
for col = 1:newVars
    nz_row = find(R((col+1):end,col));
    for row = (col + nz_row') 
        affectedRows = [col;row];
        affectedCols = (col+1):newVars;

        % Get givens rotation data and apply to positions
        [G,y] = planerot(full(R(affectedRows,col)));
        R(affectedRows,col) = sparse(y);
        R(affectedRows,affectedCols) = G*R(affectedRows,affectedCols);
        d(affectedRows) = G*d(affectedRows); 
    end
end
remainingRows = 1:newVars;
R = R(remainingRows,remainingRows);
d = d(remainingRows);