%Need proper Legendre data points
clear;
clf;

P = 4;
N = P+1;

%Vector of x values
xVec = give_me_roots(N);
xVec = xVec.';
%Vector of y values
yVec = tan(xVec);


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

    % true function
    yDataReal(m) = tan(xData(m));
end

plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on
xlabel('xVec')
ylabel('yVec')

plot(xData, yDataReal, 'k-', 'LineWidth', 4);
plot(xData, yData, 'r-', 'LineWidth', 2);
legend('training data','real', 'prediction')