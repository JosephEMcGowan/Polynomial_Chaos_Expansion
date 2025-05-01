
%loads data <--- 100x2 matrix (1st col x, 2nd col y)
load('LINE_Students_Linear_Least_Squares.mat')

%in this example we want to find a best fit line:
%       y = b0 + b1*x

%Vector of y values
yVec = DATA(:,2);

%Vector of x values
xVec = DATA(:,1);

%Creates a column vector of 1's, same size as vector x
oneVec = ones(size(xVec));

A = [oneVec xVec];

%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec
