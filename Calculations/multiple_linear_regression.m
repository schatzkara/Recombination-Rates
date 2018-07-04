function [out] = mlr(X,Y)
%This function creates a multilinear regression model of the data X
%Given by the matrix X and vector Y. The output is a vector B
%with the least-squares estimators and a scalar E,
the maximum error in the model Y=XB+E.
%   Determine the least-squares estimators
A = X'*X;
K = inv(A);
B = K*X'*Y;
%   Determine the predicted values
M = X*B;
%   Determine the random error and maximum error.
E = Y-M;
maxerr = max(abs(E));
%   Concatenate and return the output.
out = [B;maxerr];
end