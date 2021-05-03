function [eqlast_error] = Err_smooth(alpha)

global nseries  hvec kvec eps kappa kappa_X delta delta_X  nu...
    xtabm xtabp xtabX xtabm_m1 xtabp_m1 xtabX_m1 xm xp xX...
    dh cfull sfull clist slist scfull ssfull...
    anC bnC anS bnS epsbar yX csi kvec_ref sclist sslist kvec_ref...
    yTop zTop eq_option xtabTop xtabTop_m1


eqlast_error = NaN;

if((eq_option == 1)|(eq_option == 2))
    % Miller profile or Double Null equilibrium
    neqs = 6;
elseif(eq_option == 3)
    % Single Null equilibrium
    neqs = 11;
end
% This excludes the last equation, which is used for the definition of the
% error.

% Also define here the separation constants.

eteps = sqrt(1+eps^2)*alpha;

%------------Version with k_n distribution---------
%Here we are changing k_ns, which are assigned elsewhere

kvec = kvec_ref;
kvec = kvec*alpha;
hvec = sqrt((alpha^2-kvec.^2)*(1+eps^2));

%-------------------------------------------
% First part: define the series coefficients

lamtest = alpha*sqrt(epsbar*nu);

% Notice the "+1": I prefer to enter the index of the last term, rather than
% the number of terms.
% Also notice that the values of the Xc and Xs are the same for the sine
% and cosine y functions. Thus, we still only need four of them.
anC = NaN(nseries+1,4);
bnC = anC;

anS = anC;
bnS = bnC;

if(eq_option==3)
    Xcvalue = NaN(4,4); %First index is the point location
else
    Xcvalue = NaN(3,4); %First index is the point location
end

XSvalue = Xcvalue;
Xcprim = Xcvalue;
Xsprim = XSvalue;

astartC = [1,0,0]';
bstartC = [0,0,0]';

astartS = [0,0,0]';
bstartS = [1,0,0]';

for i=1:4 % it may be possible to do this with a vector call, but it's not worth it
    [anC(:,i), bnC(:,i)] = an_bn_fun_new(nseries, lamtest, epsbar, kvec(i), astartC, bstartC);
    [anS(:,i), bnS(:,i)] = an_bn_fun_new(nseries, lamtest, epsbar, kvec(i), astartS, bstartS);
end


%-------------------------------------------
% Second part: fill the arrays of "S" and "C" functions with values

% Point xm

Xcvalue(1,:) = dot(xtabm',anC).*cos(kvec*xm) + ...
    dot(xtabm',bnC).*sin(kvec*xm);

Xsvalue(1,:) = dot(xtabm',anS).*cos(kvec*xm) + ...
    dot(xtabm',bnS).*sin(kvec*xm);

Xcprim(1,:) = dot(xtabm_m1',anC).*cos(kvec*xm) + ...
    dot(xtabm_m1',bnC).*sin(kvec*xm) + ...
    kvec.*(dot(xtabm',bnC).*cos(kvec*xm) - ...
    dot(xtabm',anC).*sin(kvec*xm));

Xsprim(1,:) = dot(xtabm_m1',anS).*cos(kvec*xm) + ...
    dot(xtabm_m1',bnS).*sin(kvec*xm) + ...
    kvec.*(dot(xtabm',bnS).*cos(kvec*xm) - ...
    dot(xtabm',anS).*sin(kvec*xm));


% Point xp

Xcvalue(2,:) = dot(xtabp',anC).*cos(kvec*xp) + ...
    dot(xtabp',bnC).*sin(kvec*xp);

Xsvalue(2,:) = dot(xtabp',anS).*cos(kvec*xp) + ...
    dot(xtabp',bnS).*sin(kvec*xp);

Xcprim(2,:) = dot(xtabp_m1',anC).*cos(kvec*xp) + ...
    dot(xtabp_m1',bnC).*sin(kvec*xp) + ...
    kvec.*(dot(xtabp',bnC).*cos(kvec*xp) - ...
    dot(xtabp',anC).*sin(kvec*xp));

Xsprim(2,:) = dot(xtabp_m1',anS).*cos(kvec*xp) + ...
    dot(xtabp_m1',bnS).*sin(kvec*xp) + ...
    kvec.*(dot(xtabp',bnS).*cos(kvec*xp) - ...
    dot(xtabp',anS).*sin(kvec*xp));

% The first two points are in common for all the shapes. For the
% remaining one or two we need to differentiate depending on the
% equilibrium shape.

if(eq_option==1)
    % Point zTop
    
    Xcvalue(3,:) = dot(xtabTop',anC).*cos(kvec*zTop) + ...
        dot(xtabTop',bnC).*sin(kvec*zTop);
    
    Xsvalue(3,:) = dot(xtabTop',anS).*cos(kvec*zTop) + ...
        dot(xtabTop',bnS).*sin(kvec*zTop);
    
    Xcprim(3,:) = dot(xtabTop_m1',anC).*cos(kvec*zTop) + ...
        dot(xtabTop_m1',bnC).*sin(kvec*zTop) + ...
        kvec.*(dot(xtabTop',bnC).*cos(kvec*zTop) - ...
        dot(xtabTop',anC).*sin(kvec*zTop));
    
    Xsprim(3,:) = dot(xtabTop_m1',anS).*cos(kvec*zTop) + ...
        dot(xtabTop_m1',bnS).*sin(kvec*zTop) + ...
        kvec.*(dot(xtabTop',bnS).*cos(kvec*zTop) - ...
        dot(xtabTop',anS).*sin(kvec*zTop));
    
elseif(eq_option==2)
    % Point zTop <-> X-point
    
    Xcvalue(3,:) = dot(xtabX',anC).*cos(kvec*xX) + ...
        dot(xtabX',bnC).*sin(kvec*xX);
    
    Xsvalue(3,:) = dot(xtabX',anS).*cos(kvec*xX) + ...
        dot(xtabX',bnS).*sin(kvec*xX);
    
    Xcprim(3,:) = dot(xtabX_m1',anC).*cos(kvec*xX) + ...
        dot(xtabX_m1',bnC).*sin(kvec*xX) + ...
        kvec.*(dot(xtabX',bnC).*cos(kvec*xX) - ...
        dot(xtabX',anC).*sin(kvec*xX));
    
    Xsprim(3,:) = dot(xtabX_m1',anS).*cos(kvec*xX) + ...
        dot(xtabX_m1',bnS).*sin(kvec*xX) + ...
        kvec.*(dot(xtabX',bnS).*cos(kvec*xX) - ...
        dot(xtabX',anS).*sin(kvec*xX));
    
elseif(eq_option==3)
    % This is a repetition, but it is worth to keep things clean.
    
    % Point zTop
    
    Xcvalue(3,:) = dot(xtabTop',anC).*cos(kvec*zTop) + ...
        dot(xtabTop',bnC).*sin(kvec*zTop);
    
    Xsvalue(3,:) = dot(xtabTop',anS).*cos(kvec*zTop) + ...
        dot(xtabTop',bnS).*sin(kvec*zTop);
    
    Xcprim(3,:) = dot(xtabTop_m1',anC).*cos(kvec*zTop) + ...
        dot(xtabTop_m1',bnC).*sin(kvec*zTop) + ...
        kvec.*(dot(xtabTop',bnC).*cos(kvec*zTop) - ...
        dot(xtabTop',anC).*sin(kvec*zTop));
    
    Xsprim(3,:) = dot(xtabTop_m1',anS).*cos(kvec*zTop) + ...
        dot(xtabTop_m1',bnS).*sin(kvec*zTop) + ...
        kvec.*(dot(xtabTop',bnS).*cos(kvec*zTop) - ...
        dot(xtabTop',anS).*sin(kvec*zTop));
    
    % Point zBot <-> X-point
    
    Xcvalue(4,:) = dot(xtabX',anC).*cos(kvec*xX) + ...
        dot(xtabX',bnC).*sin(kvec*xX);
    
    Xsvalue(4,:) = dot(xtabX',anS).*cos(kvec*xX) + ...
        dot(xtabX',bnS).*sin(kvec*xX);
    
    Xcprim(4,:) = dot(xtabX_m1',anC).*cos(kvec*xX) + ...
        dot(xtabX_m1',bnC).*sin(kvec*xX) + ...
        kvec.*(dot(xtabX',bnC).*cos(kvec*xX) - ...
        dot(xtabX',anC).*sin(kvec*xX));
    
    Xsprim(4,:) = dot(xtabX_m1',anS).*cos(kvec*xX) + ...
        dot(xtabX_m1',bnS).*sin(kvec*xX) + ...
        kvec.*(dot(xtabX',bnS).*cos(kvec*xX) - ...
        dot(xtabX',anS).*sin(kvec*xX));
    
end

%-------------------------------------------
% Third part: fill the matrix of coefficients and right-hand side
% Very important: the order of variables is {c2,s2,c3,s3,c4,s4,sc1,sc2,sc3,ss2,ss3}
% Note that we are using different names with respect to Jeff's writeup for
% consistency with previous versions and to avoid confusion with the an and
% bn coefficients of the series expansion.

mcoef = NaN(neqs);
rhs = NaN(neqs,1);

indC = [1 3 5];
indS = indC + 1;

if(eq_option==3)
    indsC = [7 8 9];
    indsS = [10 11];
end

eqind = 0;

% Equation 1 (a): psi(-1,0) = 0
% This one is the same for all geometries.

eqind = eqind+1;

mcoef(eqind,indC) = Xcvalue(1,2:4);
mcoef(eqind,indS) = Xsvalue(1,2:4);

if(eq_option==3)
    mcoef(eqind,indsC) = 0;
    mcoef(eqind,indsS) = 0;
end

rhs(eqind) = -Xcvalue(1,1);


% Equation 2 (b): psi(1,0) = 0
% This one is also the same for all geometries.

eqind = eqind+1;

mcoef(eqind,indC) = Xcvalue(2,2:4);
mcoef(eqind,indS) = Xsvalue(2,2:4);

if(eq_option==3)
    mcoef(eqind,indsC) = 0;
    mcoef(eqind,indsS) = 0;
end

rhs(eqind) = -Xcvalue(2,1);


% Equation 3 (c): psi(z_top,y_top) = 0

eqind = eqind+1;

mcoef(eqind,indC) = Xcvalue(3,2:4).*cos(yTop*hvec(2:4));
mcoef(eqind,indS) = Xsvalue(3,2:4).*cos(yTop*hvec(2:4));

if(eq_option==3)
    mcoef(eqind,indsC) = Xcvalue(3,1:3).*sin(yTop*hvec(1:3));
    mcoef(eqind,indsS) = Xsvalue(3,2:3).*sin(yTop*hvec(2:3));
end

rhs(eqind) = -Xcvalue(3,1).*cos(yTop*hvec(1));

if(eq_option==3)
    % Equation 4 (d)
    % The new Eq. (d) is imposed in the lower X point
    
    eqind = eqind+1;
    
    mcoef(eqind,indC) = Xcvalue(4,2:4).*cos(yX*hvec(2:4));
    mcoef(eqind,indS) = Xsvalue(4,2:4).*cos(yX*hvec(2:4));
    mcoef(eqind,indsC) = Xcvalue(4,1:3).*sin(yX*hvec(1:3));
    mcoef(eqind,indsS) = Xsvalue(4,2:3).*sin(yX*hvec(2:3));
    
    rhs(eqind) = -Xcvalue(4,1).*cos(yX*hvec(1));
    
end

if(eq_option==3)
    % Equation 5 (e): 0 derivative at (-1,0)
    % This one is new.
    
    eqind = eqind+1;
    
    mcoef(eqind,indC) = 0;
    mcoef(eqind,indS) = 0;
    mcoef(eqind,indsC) = Xcvalue(1,1:3).*hvec(1:3);
    mcoef(eqind,indsS) = Xsvalue(1,2:3).*hvec(2:3);
    
    rhs(eqind) = 0;
    
end

if(eq_option==3)
    
    % Equation 6 (f): 0 derivative at (1,0)
    % This one is new.
    
    eqind = eqind+1;
    
    mcoef(eqind,indC) = 0;
    mcoef(eqind,indS) = 0;
    mcoef(eqind,indsC) = Xcvalue(2,1:3).*hvec(1:3);
    mcoef(eqind,indsS) = Xsvalue(2,2:3).*hvec(2:3);
    
    rhs(eqind) = 0;
    
end

% Equation 4 (d): psi_z(z_top,y_top) = 0
% This is Equation 7 for the SN case

eqind = eqind + 1;

if(eq_option==1)
    mcoef(eqind,indC) = Xcprim(3,2:4).*cos(yTop*hvec(2:4));
    mcoef(eqind,indS) = Xsprim(3,2:4).*cos(yTop*hvec(2:4));
elseif(eq_option==2)
    mcoef(eqind,indC) = Xcprim(3,2:4).*cos(yX*hvec(2:4));
    mcoef(eqind,indS) = Xsprim(3,2:4).*cos(yX*hvec(2:4));
elseif(eq_option==3)
    mcoef(eqind,indC) = Xcprim(3,2:4).*cos(yTop*hvec(2:4));
    mcoef(eqind,indS) = Xsprim(3,2:4).*cos(yTop*hvec(2:4));
    mcoef(eqind,indsC) = Xcprim(3,1:3).*sin(yTop*hvec(1:3));
    mcoef(eqind,indsS) = Xsprim(3,2:3).*sin(yTop*hvec(2:3));
end

if(eq_option==2)
    rhs(eqind) = -Xcprim(3,1).*cos(yX*hvec(1));
else
    rhs(eqind) = -Xcprim(3,1).*cos(yTop*hvec(1));
end

if(eq_option==3)
    
    % Equation 8 (h): psi_z(xX,yX)=0
    
    
    eqind = eqind+1;
    
    mcoef(eqind,indC) = Xcprim(4,2:4).*cos(yX*hvec(2:4));
    mcoef(eqind,indS) = Xsprim(4,2:4).*cos(yX*hvec(2:4));
    mcoef(eqind,indsC) = Xcprim(4,1:3).*sin(yX*hvec(1:3));
    mcoef(eqind,indsS) = Xsprim(4,2:3).*sin(yX*hvec(2:3));
    
    rhs(eqind) = -Xcprim(4,1).*cos(yX*hvec(1));
    
end

if(eq_option==3)
    % Equation 9 (i): psi_y(xX,yX)=0
    % This is to define the X-point [this is Eq. (g) in the DN version].
    
    eqind = eqind+1;
    % CHECK SIGNS!!!!!
    
    mcoef(eqind,indC) =  -hvec(2:4).*sin(hvec(2:4)*yX).*Xcvalue(4,(2:4));
    mcoef(eqind,indS) = -hvec(2:4).*sin(hvec(2:4)*yX).*Xsvalue(4,2:4);
    mcoef(eqind,indsC) = hvec(1:3).*cos(hvec(1:3)*yX).*Xcvalue(4,(1:3));
    mcoef(eqind,indsS) = hvec(2:3).*cos(hvec(2:3)*yX).*Xsvalue(4,2:3);
    
    rhs(eqind) = hvec(1).*sin(hvec(1)*yX).*Xcvalue(4,(1));
    
end

% Equation 10 (j): assign curvature in (-1,0)
% This is Equation 5 (e) for the symmetrical cases.

eqind = eqind+1;

% Note: we define this with the opposite sign with respect to Jeff's notes.
if(eq_option==1)
    Lambda1 = -(1-eps)*(1-dh)^2/(kappa^2);
elseif(eq_option==2)
    Lambda1 = -(1-eps)*(1-delta_X)*(1+csi)/(kappa_X^2);
elseif(eq_option==3)
    Lambda1 =  -(1-eps)*(1-dh)^2/(kappa^2)-(1-eps)*(1-delta_X)*(1+csi)/(kappa_X^2);
    Lambda1 = Lambda1/2;
end

mcoef(eqind,indC) = Xcvalue(1,2:4).*(hvec(2:4).^2) + ...
    Lambda1*Xcprim(1,2:4);
mcoef(eqind,indS) = Xsvalue(1,2:4).*(hvec(2:4).^2 )+ ...
    Lambda1*Xsprim(1,2:4);

if(eq_option==3)
    mcoef(eqind,indsC) = 0;
    mcoef(eqind,indsS) = 0;
end


rhs(eqind) = -(Xcvalue(1,1)*hvec(1)^2 + Lambda1*Xcprim(1,1));


% Equation 11 (k): assign curvature in (1,0)
% This is Equation 6 (f) for the symmetrical cases.

eqind = eqind+1;

% Here the sign is the same as in Jeff's notes
if(eq_option==1)
    Lambda2 = (1+eps)*(1+dh)^2/(kappa^2);
elseif(eq_option==2)
    Lambda2 = (1+eps)*(1+delta_X)*(1+csi)/(kappa_X^2);
elseif(eq_option==3)
    Lambda2 = (1+eps)*(1+dh)^2/(kappa^2)+(1+eps)*(1+delta_X)*(1+csi)/(kappa_X^2);
    Lambda2 = Lambda2/2;
end

mcoef(eqind,indC) = Xcvalue(2,2:4).*(hvec(2:4).^2) + ...
    Lambda2*Xcprim(2,2:4);
mcoef(eqind,indS) = Xsvalue(2,2:4).*(hvec(2:4).^2 )+ ...
    Lambda2*Xsprim(2,2:4);

if(eq_option==3)
    mcoef(eqind,indsC) = 0;
    mcoef(eqind,indsS) = 0;
end

rhs(eqind) = -(Xcvalue(2,1)*hvec(1)^2 + Lambda2*Xcprim(2,1));

%-------------------------------------------
% Fourth part: solve the linear system and replace values in equation 7 (g)
% This is the main difference from the Miller profile version

lvec = NaN(neqs,1);

solvec = linsolve(mcoef,rhs);

clist = solvec(indC);
slist = solvec(indS)';

if(eq_option==3)
    sclist = solvec(indsC)';
    sslist = solvec(indsS)';
end

%curv = kappa/((1-eps*delta)^2*(1-delta^2));

% Note that Eq. (g) has seven or twelve coefficients. The first one was on the
% right-hand side in the previous equations. To keep things in order, we
% shift all indexes by 1.
% This can be done in a more elegant manner defining an additional vector.

sfull = [0, slist];

cfull = [1; clist]';
indCfull = [1 indC+1];

if(eq_option==3)
    scfull = sclist;
    ssfull = [0, sslist];
end

if((eq_option==1)||(eq_option==3))
    
    Lambda3 = kappa/((1-eps*delta)^2*(1-delta^2));
    
    d2term = -(kvec.^2+lamtest^2*zTop)/(1+epsbar*zTop);
    
    lvec(indCfull) = (cos(yTop*hvec).*d2term + ...
        Lambda3*hvec.*sin(yTop*hvec)).*cfull.*Xcvalue(3,:);
    lvec(indS+1) = (cos(yTop*hvec(2:4)).*d2term(2:4) + ...
        Lambda3*hvec(2:4).*sin(yTop*hvec(2:4))).*slist.*Xsvalue(3,2:4);
    
    if(eq_option==3)
        lvec(indsC+1) = (sin(yTop*hvec(1:3)).*d2term(1:3) - ...
            Lambda3*hvec(1:3).*cos(yTop*hvec(1:3))).*sclist.*Xcvalue(3,1:3);
        lvec(indsS+1) = (sin(yTop*hvec(2:3)).*d2term(2:3) - ...
            Lambda3*hvec(2:3).*cos(yTop*hvec(2:3))).*sslist.*Xsvalue(3,2:3);
    end
    
else
    
    lvec(indCfull) = hvec.*sin(hvec*yX).*cfull.*Xcvalue(3,:);
    lvec(indS+1) = hvec(2:4).*sin(hvec(2:4)*yX).*slist.*Xsvalue(3,2:4);
    
end

eqlast_error = (sum(lvec)/sum(abs(lvec)))^2;


end