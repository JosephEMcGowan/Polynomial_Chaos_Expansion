%Need proper Legendre data points
clear

N = 10;
P = 5;


%Vector of x values
xVec = linspace(0,2*pi,N);

%Vector of y values
yVec = sin(xVec)';

%transforms x from [0, 2pi] to [-1,1]
xVec = (1/pi)*xVec - 1;


%creates our A matrix for the function
for n=1:N
    for m=1:P+1
        A(n,m) = legendreP(m-1,xVec(n));
    end
end    

%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec;



%puts data in order
xData = linspace(-1,1,200);

for m = 1:length(xData)
    x = xData(m);
    sum = 0;
    for n = 1:P+1
        sum = sum + betaVec(n) * legendreP(n-1,x);
    end
    yData(m) = sum;
end

plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on
xlabel('xVec')
ylabel('yVec')

plot(xData, yData, 'r-', 'LineWidth', 2);
legend('training data', 'pred')