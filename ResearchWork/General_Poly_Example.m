%loads data <--- 100x2 matrix (1st col x, 2nd col y)
load('CUBIC_Students_Linear_Least_Squares.mat')

N = 5;

%Vector of y values
yVec = DATA(:,2);

%Vector of x values
xVec = DATA(:,1);

%Creates a column vector of 1's, same size as vector x
oneVec = ones(size(xVec));

A = [oneVec];
for n=1:N
    A = [A xVec.^n];
end    

%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec

%makes a plot

xData = linspace(-2, 4, 100);

for m = 1:length(xData)
    x = xData(m);
    sum = 0;
    for n = 0:N
        sum = sum + betaVec(n+1) * x.^n;
    end
    yData(m) = sum;
end

plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on

plot(xData, yData, 'r-', 'LineWidth', 2);