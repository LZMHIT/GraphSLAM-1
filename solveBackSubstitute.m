% solve system of eqns Rx = d given upper triangular matrix R
function x = solveBackSubstitute(R,d)
[nz_row,nz_col,nz_val] = find(R); % find nonzero entries
numNZ = length(nz_col);
x = d; % initially

% Loop through nonzeros entries
for k = numNZ:-1:1
    % Extract nonzero row, col, and value
    row = nz_row(k);
    col = nz_col(k);
    val = nz_val(k);
    
    if row == col
        x(row) = x(row)/val;
    elseif col > row
        x(row) = x(row) - val*x(col);
    else
        error('R must be upper triangular for back substitution');
    end
end