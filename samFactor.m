% function to factor sqrt SAM matrix
function [R,d] = samFactor(A,b)
R_aug = qr([A,b],0); % Q-less QR factorization. Augment with b to update b without Q
                     % Economy factorization leaves out zero rows so R is
                     % square
numCols = size(A,2);
numRows = min(numCols,size(A,1));
R_col_range = 1:numCols;
R_row_range = 1:numRows;
R = R_aug(R_row_range,R_col_range); % R is all but last row/column of augmented R
d = R_aug(R_row_range,end); % d is last column of augmented R without last row

