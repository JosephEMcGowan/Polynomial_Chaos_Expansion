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
        for r3 = 1:length(xVec)
        ct = ct + 1;
        rootMAT(ct,1) = xVec(r1);
        rootMAT(ct,2) = xVec(r2);
        rootMAT(ct,3) = xVec(r3);
        end
    end
end
fprintf('finished making the root matrix\n')
rootMAT;

%Vector of y values, sets the shape of the function we "draw"
%yVec = (1/8)*( ( 3.*( (1/2).*rootMAT(:,1)+(1/2)).^2 + 1 ) .*( 3.*( (1/2).*rootMAT(:,2)+(1/2)).^2 + 1 ) .* ( 3.*((1/2).*rootMAT(:,3)+(1/2)).^2 + 1));

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
alphaMAT

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

IndsVec123

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

S2 = 0;
for ct = 1:length(IndsVec2)
S2 = S2 + betaVec(IndsVec2(ct)).^2*(1/(2*alphaMAT(IndsVec2(ct),2)+1));
end

S2;
SU2 = S2/SVar;
SU2;

S3 = 0;
for ct = 1:length(IndsVec3)
S3 = S3 + betaVec(IndsVec3(ct)).^2*(1/(2*alphaMAT(IndsVec3(ct),3)+1));
end

S3;
SU3 = S3/SVar;
SU3;


S12 = 0;
for ct = 1:length(IndsVec12)
S12 = S12 + betaVec(IndsVec12(ct)).^2*(1/(2*alphaMAT(IndsVec12(ct),1)+1))*(1/(2*alphaMAT(IndsVec12(ct),2)+1));
end

SU12 = S12/SVar;
SU12;

S12 = 0;
for ct = 1:length(IndsVec12)
S12 = S12 + betaVec(IndsVec12(ct)).^2*(1/(2*alphaMAT(IndsVec12(ct),1)+1))*(1/(2*alphaMAT(IndsVec12(ct),2)+1));
end

SU12 = S12/SVar;
SU12;

S12 = 0;
for ct = 1:length(IndsVec12)
S12 = S12 + betaVec(IndsVec12(ct)).^2*(1/(2*alphaMAT(IndsVec12(ct),1)+1))*(1/(2*alphaMAT(IndsVec12(ct),2)+1));
end

SU12 = S12/SVar;
SU12;



S123 = 0;
for ct = 1:length(IndsVec123)
S123 = S123 + betaVec(IndsVec123(ct)).^2*(1/(2*alphaMAT(IndsVec123(ct),1)+1))*(1/(2*alphaMAT(IndsVec123(ct),2)+1))*(1/(2*alphaMAT(IndsVec123(ct),3)+1));
end

SU123 = S123/SVar;
SU123