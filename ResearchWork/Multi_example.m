clear;
clf;

%sets the amount of legendre polynomials we will be using
P = 3;                                                                                                                                                                                                                                                                                       

%Vector of x values, picks the amount of poly we chose +1 roots
xVec = give_me_roots(P+1);


%creates rootMAT
ct = 0;
for r1 = 1:length(xVec)
    for r2 = 1:length(xVec)
        ct = ct + 1;
        rootMAT(ct,1) = xVec(r1);
        rootMAT(ct,2) = xVec(r2);
    end
end
rootMAT

%Vector of y values, sets the shape of the function we "draw"
yVec = (rootMAT(:,1)).^3 + (rootMAT(:,2)).^3;


%forms alphaMAT matrix
ct = 0;
for a1 = 0:P
    for a2 = 0:P
        if a1+a2 <= P
            ct = ct + 1;
            alphaMAT(ct,1) = a1;
            alphaMAT(ct,2) = a2;
        end
    end
end
alphaMAT


%creates our A matrix for the function
for n=1:length(rootMAT(:,1))
    for m=1:length(alphaMAT(:,1))
        A(n,m) = legendreP(alphaMAT(m,1),rootMAT(n,1))*legendreP(alphaMAT(m,2),rootMAT(n,2));
    end
end    
A


%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec;

betaVec

%-------------------------------------------------------------------------

%puts data in order
xData = linspace(-1,1,10);

%This takes in x and y data and forms the equation to solve for the betaVec
%components
ct = 0;
for m1 = 1:length(xData)
    for m2 = 1:length(xData)
        x1 = xData(m1);
        x2 = xData(m2);
        sum = 0;
        for n = 1:length(alphaMAT(:,1))
            sum = sum + betaVec(n) * legendreP(alphaMAT(n,1),x1) * legendreP(alphaMAT(n,2), x2);
        end
        ct = ct+1;
        yData(ct) = sum;
        % true function
        yDataReal(ct) = (x1).^3 + (x2).^3;
    end
end


%plot(xVec, yVec, '.', 'MarkerSize',20, 'Color', [0.2 0.6 0.7]); hold on
%xlabel('xVec')
%ylabel('yVec')


%This outputs our graphs but it ignores the xdata since there are more
%values then it would see outputted, looks for 10 but there are 100 also it
%is not the 2D graph we want, just normal
plot(yDataReal, 'k-', 'LineWidth', 4); hold on;
plot(yData, 'r-', 'LineWidth', 2);
%legend('training data','real', 'prediction')

