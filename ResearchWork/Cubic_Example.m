%loads data <--- 100x2 matrix (1st col x, 2nd col y)
load('CUBIC_Students_Linear_Least_Squares.mat')

%in this example we want to find a best fit line:
%       y = b0 + b1*x

%Vector of y values
yVec = DATA(:,2);

%Vector of x values
xVec = DATA(:,1);

%Creates a column vector of 1's, same size as vector x
oneVec = ones(size(xVec));

A = [oneVec xVec xVec.^2 xVec.^3];

%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec

%makes a plot
plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on


xData = linspace(-2, 4, 100);
%xData = -2:0.1:4; <- another way to do it

yData = betaVec(1) + betaVec(2) * xData + betaVec(3) * xData.^2 + betaVec(4) * xData.^3;

plot(xData, yData, 'r-', 'LineWidth', 2);