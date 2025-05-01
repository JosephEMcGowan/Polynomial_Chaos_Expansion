clear;
clf;

%sets the amount of legendre polynomials we will be using
P = 5;                                                                                                                                                                                                                                                                                       

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
rootMAT;

%Vector of y values, sets the shape of the function we "draw"
yVec = tan(rootMAT(:,1)) + sin(5*rootMAT(:,2));


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
alphaMAT;


%creates our A matrix for the function
for n=1:length(rootMAT(:,1))
    for m=1:length(alphaMAT(:,1))
        A(n,m) = legendreP(alphaMAT(m,1),rootMAT(n,1))*legendreP(alphaMAT(m,2),rootMAT(n,2));
    end
end    
A;


%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec;

betaVec;

%-------------------------------------------------------------------------

%puts data in order and makes a 2D domain
x1_Data = linspace(-1,1,10);
x2_Data = linspace(-1,1,10);
[X1,X2] = meshgrid(x1_Data, x2_Data);

%This takes in x and y data and forms the equation to solve for the betaVec
%components
ct = 0;
for m1 = 1:length(x1_Data)
    for m2 = 1:length(x2_Data)
        x1 = x1_Data(m1);
        x2 = x2_Data(m2);
        sum = 0;
        for n = 1:length(alphaMAT(:,1))
            ct = ct + 1;
            sum = sum + betaVec(n) * legendreP(alphaMAT(n,1),x1) * legendreP(alphaMAT(n,2), x2);
        end
        gpc(m2,m1) = sum;
        % true function
        gpcReal(m2,m1) = tan(x1) + sin(5*x2);
    end
end

%Graphs the 3D image of our prediction
figure(1)
surf(X1,X2,gpc)
xlabel('x1')
ylabel('x2')
zlabel('gpc Prediction')
colorbar

%Grapgs the 3D image of the error between our preciction and the real
figure(2)
surf(X1,X2,abs(gpc-gpcReal))
xlabel('x1')
ylabel('x2')
zlabel('Error Between')
colorbar