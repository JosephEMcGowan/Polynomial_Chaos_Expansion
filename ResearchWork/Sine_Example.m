%loads data <--- 100x2 matrix (1st col x, 2nd col y)
load('SINE_EXP_Linear_Least_Squares.mat')

%in this example we want to find a best fit line:
%       y = b0 + b1*x

%Vector of y values
yVec = DATA(:,2);

%Vector of x values
xVec = DATA(:,1);

%Creates a column vector of 1's, same size as vector x
oneVec = ones(size(xVec));

A = [oneVec sin(2*xVec) exp(xVec / 4)];

%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec

%makes a plot
plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on

xData = linspace(0, 20, 100);
%xData = 0:0.1:20; <- another way to do it

yData = betaVec(1) + betaVec(2) * sin(2*xData) + betaVec(3) * exp(xData/4);

plot(xData, yData, 'r-', 'LineWidth', 2);