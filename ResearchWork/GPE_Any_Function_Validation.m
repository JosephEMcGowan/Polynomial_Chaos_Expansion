clear;
clf;

tic

%sets the amount of legendre polynomials we will be using
P = 3;                                                                                                                                                                                                                                                                                       

%Vector of x values, picks the amount of poly we chose +1 roots
xVec = give_me_roots(P+1);


%creates rootMAT
ct = 0;
for r1 = 1:length(xVec)
    for r2 = 1:length(xVec)
        for r3 = 1:length(xVec)
        ct = ct + 1;
        rootMAT(ct,1) = xVec(r1);
        rootMAT(ct,2) = xVec(r2);
        rootMAT(ct,3) = xVec(r3);
        end
    end
end

%this for loop takes the columns of the root matrix and takes the value of
%the distance formula using the columns values and puts them into a vector
for r1 = 1:length(rootMAT(:,1))
    distVec(r1) = sqrt((rootMAT(r1,1).^2)+(rootMAT(r1,2).^2)+(rootMAT(r1,3).^2));
end

%this sorts both our distance vector and our root matrix indices, this
%allows us to then see an organized version of our root matrix
[SortedDistVec,rootInd] = sort(distVec);
rootMAT = rootMAT(rootInd,:);

%this simplifies the amount of data we take from our root matrix, as we do
%not need all data from the given roots, just some. This is one of the big
%advantages of polynomial chaos, it allows us to get more from less
M = 3;
%NSamples = (M-1)*(factorial(M+P)/(factorial(M)*factorial(P)));
NSamples = (P+1)^3;
NSamples
rootMAT = rootMAT(1:NSamples, :);

fprintf('finished making the sorted distance vector and sorted root indices vector\n')
fprintf('finished making the sorted and simplified root matrix\n')
SortedDistVec;
rootInd;
rootMAT;



%Vector of y values, sets the shape of the function we "draw"
%polynomial example
%yVec = (1/8)*( ( 3.*( (1/2).*rootMAT(:,1)+(1/2)).^2 + 1 ) .*( 3.*( (1/2).*rootMAT(:,2)+(1/2)).^2 + 1 ) .* ( 3.*((1/2).*rootMAT(:,3)+(1/2)).^2 + 1));

%Ishigami function example where a = 7 and b = 0.1
yVec = sin(pi .* rootMAT(:,1)) + (7*(sin(pi .* rootMAT(:,2)).^2)) + (((0.1).*(pi .* rootMAT(:,3)).^4) .* (sin(pi .* rootMAT(:,1))));

%forms alphaMAT matrix
ct = 0;
for a1 = 0:P
    for a2 = 0:P
        for a3 = 0:P
            if a1+a2+a3 <= P
                ct = ct + 1;
                alphaMAT(ct,1) = a1;
                alphaMAT(ct,2) = a2;
                alphaMAT(ct,3) = a3;
            end
        end      
    end
end
fprintf('finished making the alpha matrix\n')
alphaMAT;

%creates our A matrix for the function
for n=1:length(rootMAT(:,1))
    for m=1:length(alphaMAT(:,1))
        A(n,m) = legendreP(alphaMAT(m,1),rootMAT(n,1))*legendreP(alphaMAT(m,2),rootMAT(n,2))*legendreP(alphaMAT(m,3),rootMAT(n,3));
    end
end  
fprintf('finished making the A matrix\n')
A;


%Get the coefficients for our model function
%   Through pseudo-inverse
betaVec = inv(A'*A)*A'*yVec;

betaVec;

SVar = 0;
for totalVar=2:length(betaVec)
    SVar = SVar + (betaVec(totalVar)).^2 .* 1/(2.*alphaMAT(totalVar,1) + 1)  .* 1/(2.*alphaMAT(totalVar,2) + 1) .* 1/(2.*alphaMAT(totalVar,3) + 1);
end

SVar;

%This creates and stores an indicies vector that finds and stores when the
%second and third columns of the alphaMAT are both 0, this keeps varrying the
%first term but keeps the second and third columns equal to 0. This set of
%code works to find only when the first column varries and the second and
%third are 0.
ct = 0;
for S1=2:length(alphaMAT)
    if alphaMAT(S1,2) + alphaMAT(S1,3) == 0
        ct = ct+1;
        IndsVec1(ct) = S1;
    end
end

IndsVec1;

%This creates and stores an indicies vector that finds and stores when the
%first and third columns of the alphaMAT are both 0, this keeps varrying the
%second term but keeps the first and third columns equal to 0. This set of
%code works to find only when the second column varries and the first and
%third are 0.
ct = 0;
for S2=2:length(alphaMAT)
    if alphaMAT(S2,1) + alphaMAT(S2,3) == 0
        ct = ct+1;
        IndsVec2(ct) = S2;
    end
end

IndsVec2;

%This creates and stores a first order index vector that finds and stores when the
%second and first columns of the alphaMAT are both 0, this keeps varrying the
%third term but keeps the second and first columns equal to 0. This set of
%code works to find only when the third column varries and the second and
%first are 0.
ct = 0;
for S3=2:length(alphaMAT)
    if alphaMAT(S3,1) + alphaMAT(S3,2) == 0
        ct = ct+1;
        IndsVec3(ct) = S3;
    end
end

IndsVec3;


%This is creating the second order index vector where both columns one and two are
%being varried, and column 3 is left at 0. This picks all instances where
%columns 1 and 2 are non zero and column 3 is zero.
ct = 0;
for S12=2:length(alphaMAT)
    if alphaMAT(S12,1) ~= 0 && alphaMAT(S12,2) ~= 0 &&alphaMAT(S12,3) == 0
        ct = ct+1;
        IndsVec12(ct) = S12;
    end
end

IndsVec12;

%This is creating the second order index vector where both columns 1 and 3 are
%being varried, and column 2 is left at 0. This picks all instances where
%columns 1 and 3 are non zero and column 2 is zero.
ct = 0;
for S13=2:length(alphaMAT)
    if alphaMAT(S13,1) ~= 0 && alphaMAT(S13,3) ~= 0 &&alphaMAT(S13,2) == 0
        ct = ct+1;
        IndsVec13(ct) = S13;
    end
end

IndsVec13;

%This is creating the second order index vector where both columns 2 and 3 are
%being varried, and column 1 is left at 0. This picks all instances where
%columns 2 and 3 are non zero and column 1 is zero.
ct = 0;
for S23=2:length(alphaMAT)
    if alphaMAT(S23,2) ~= 0 && alphaMAT(S23,3) ~= 0 &&alphaMAT(S23,1) == 0
        ct = ct+1;
        IndsVec23(ct) = S23;
    end
end

IndsVec23;

%This is creating the third order index vector where both columns one and
%two and three are being varried. This picks all instances where
%columns 1 and 2 and 3 are non zero.
ct = 0;
for S123=2:length(alphaMAT)
    if alphaMAT(S123,1) ~= 0 && alphaMAT(S123,2) ~= 0 && alphaMAT(S123,3) ~= 0
        ct = ct+1;
        IndsVec123(ct) = S123;
    end
end

IndsVec123;


%This finds and makes a vector of all indicies where the first column is not
%equal to 0. Here we only care about the first column, this gives us the
%total order of all cases where column one is changing
ct = 0;
for S1T=2:length(alphaMAT)
    if alphaMAT(S1T,1) ~= 0
        ct = ct+1;
        IndsVec1T(ct) = S1T;
    end
end

IndsVec1T;


%This finds and makes a vector of all indicies where the second column is not
%equal to 0. Here we only care about the second column, this gives us the
%total order of all cases where column two is changing
ct = 0;
for S2T=2:length(alphaMAT)
    if alphaMAT(S2T,2) ~= 0
        ct = ct+1;
        IndsVec2T(ct) = S2T;
    end
end

IndsVec2T;


%This finds and makes a vector of all indicies where the third column is not
%equal to 0. Here we only care about the third column, this gives us the
%total order of all cases where column three is changing
ct = 0;
for S3T=2:length(alphaMAT)
    if alphaMAT(S3T,3) ~= 0
        ct = ct+1;
        IndsVec3T(ct) = S3T;
    end
end

IndsVec3T;

%This computes the value of S1 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S1 is the first
%order partial variance for trait 1 which corresponds to column 1 of the
%alphaMAT. 
S1 = 0;
for ct = 1:length(IndsVec1)
S1 = S1 + betaVec(IndsVec1(ct)).^2*(1/(2*alphaMAT(IndsVec1(ct),1)+1));
end

S1;
SU1 = S1/SVar;
SU1;


%This computes the value of S2 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S2 is the first
%order partial variance for trait 2 which corresponds to column 2 of the
%alphaMAT. 
S2 = 0;
for ct = 1:length(IndsVec2)
S2 = S2 + betaVec(IndsVec2(ct)).^2*(1/(2*alphaMAT(IndsVec2(ct),2)+1));
end

S2;
SU2 = S2/SVar;
SU2;


%This computes the value of S3 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S3 is the first
%order partial variance for trait 3 which corresponds to column 3 of the
%alphaMAT. 
S3 = 0;
for ct = 1:length(IndsVec3)
S3 = S3 + betaVec(IndsVec3(ct)).^2*(1/(2*alphaMAT(IndsVec3(ct),3)+1));
end

S3;
SU3 = S3/SVar;
SU3;



%This computes the value of S12 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S12 is the first
%order partial variance for trait 1 and 2 which corresponds to column 1 and 2 of the
%alphaMAT. 
S12 = 0;
for ct = 1:length(IndsVec12)
S12 = S12 + betaVec(IndsVec12(ct)).^2*(1/(2*alphaMAT(IndsVec12(ct),1)+1))*(1/(2*alphaMAT(IndsVec12(ct),2)+1));
end

SU12 = S12/SVar;
SU12;


%This computes the value of S13 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S13 is the first
%order partial variance for trait 1 and 3 which corresponds to column 1 and 3 of the
%alphaMAT. 
S13 = 0;
for ct = 1:length(IndsVec13)
S13 = S13 + betaVec(IndsVec13(ct)).^2*(1/(2*alphaMAT(IndsVec13(ct),1)+1))*(1/(2*alphaMAT(IndsVec13(ct),3)+1));
end

SU13 = S13/SVar;
SU13;


%This computes the value of S23 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S23 is the first
%order partial variance for trait 2 and 3 which corresponds to column 2 and 3 of the
%alphaMAT. 
S23 = 0;
for ct = 1:length(IndsVec23)
S23 = S23 + betaVec(IndsVec23(ct)).^2*(1/(2*alphaMAT(IndsVec23(ct),2)+1))*(1/(2*alphaMAT(IndsVec23(ct),3)+1));
end

SU23 = S23/SVar;
SU23;


%This computes the value of S123 with the found betaVec based on the given
%yVec, using the rows of alphaMAT marked by the IndsVec. S123 is the first
%order partial variance for trait 1, 2 and 3 which corresponds to column 1, 2 and 3 of the
%alphaMAT. 
S123 = 0;
for ct = 1:length(IndsVec123)
S123 = S123 + betaVec(IndsVec123(ct)).^2*(1/(2*alphaMAT(IndsVec123(ct),1)+1))*(1/(2*alphaMAT(IndsVec123(ct),2)+1))*(1/(2*alphaMAT(IndsVec123(ct),3)+1));
end

SU123 = S123/SVar;
SU123;

%This finds the total order of all indicies where we only look at where 1
%is nonzero, we dont care about 2 or 3
S1T = 0;
for ct = 1:length(IndsVec1T)
S1T = S1T + betaVec(IndsVec1T(ct)).^2*(1/(2*alphaMAT(IndsVec1T(ct),1)+1))*(1/(2*alphaMAT(IndsVec1T(ct),2)+1))*(1/(2*alphaMAT(IndsVec1T(ct),3)+1));
end

S1T;
SU1T = S1T/SVar;
SU1T;

%This finds the total order of all indicies where we only look at where 2
%is nonzero, we dont care about 1 or 3
S2T = 0;
for ct = 1:length(IndsVec2T)
S2T = S2T + betaVec(IndsVec2T(ct)).^2*(1/(2*alphaMAT(IndsVec2T(ct),1)+1))*(1/(2*alphaMAT(IndsVec2T(ct),2)+1))*(1/(2*alphaMAT(IndsVec2T(ct),3)+1));
end

S2T;
SU2T = S2T/SVar;
SU2T;

%This finds the total order of all indicies where we only look at where 3
%is nonzero, we dont care about 1 or 2
S3T = 0;
for ct = 1:length(IndsVec3T)
S3T = S3T + betaVec(IndsVec3T(ct)).^2*(1/(2*alphaMAT(IndsVec3T(ct),1)+1))*(1/(2*alphaMAT(IndsVec3T(ct),2)+1))*(1/(2*alphaMAT(IndsVec3T(ct),3)+1));
end

S3T;
SU3T = S3T/SVar;
SU3T;

%---------------------------------------------------------------------------------

%puts data in order and makes a 2D domain
x1_Data = linspace(-1,1,20);
x2_Data = linspace(-1,1,20);
x3 = 0.26;
[X1,X2] = meshgrid(x1_Data, x2_Data);

%This takes in x and y data and forms the equation to solve for the betaVec
%components
for m1 = 1:length(x1_Data)
    for m2 = 1:length(x2_Data)
        x1 = x1_Data(m1);
        x2 = x2_Data(m2);
        sum = 0;
        for n = 1:length(alphaMAT(:,1))
            sum = sum + betaVec(n) * legendreP(alphaMAT(n,1),x1) * legendreP(alphaMAT(n,2), x2) * legendreP(alphaMAT(n,3), x3);
        end
        gpc(m2,m1) = sum;
        % true function choose between poly and ishigami
        %gpcRealPoly(m2,m1) = (1/8)*( ( 3.*( (1/2).* x1 + (1/2)).^2 + 1 ) .*( 3.*( (1/2).* x2 +(1/2)).^2 + 1 ) .* ( 3.*((1/2).* x3 +(1/2)).^2 + 1));
        gpcRealIshi(m2,m1) = sin(pi .* x1) + (7*(sin(pi .* x2.^2))) + (((0.1).*(pi .* x3).^4) .* (sin(pi .* x1)));
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
%choose a surf pending on if you are using polynomial function or ishigami
figure(2)
%surf(X1,X2,abs(gpc-gpcRealPoly))
surf(X1,X2,abs(gpc-gpcRealIshi))
xlabel('x1')
ylabel('x2')
zlabel('Error Between')
colorbar

toc