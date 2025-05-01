%Need proper sine data points

N = 200;

%Vector of x values
xVec = 2*pi*rand(1000,1);

%Vector of y values
yVec = xVec.^3;



%Creates a column vector of 1's, same size as vector x
zeroVec = zeros(size(xVec));

A = [];
for n=1:N
    A = [A sin((n*xVec)/2)];
end    

%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec;

%makes a plot

xData = linspace(0, 2*pi, 100);

for m = 1:length(xData)
    x = xData(m);
    sum = 0;
    for n = 1:N
        sum = sum + betaVec(n) * sin((n*x)/2);
    end
    yData(m) = sum;
end

plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on

plot(xData, yData, 'r-', 'LineWidth', 2);