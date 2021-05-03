function make_analytic_eqdsk(eq_option_in,R0_in,eps_in,delta_X_in,kappa_X_in,delta_in,kappa_in,nu_in,B0)
% MAKE_ANALYTIC_EQDSK(eq_option_in,R0_in,eps_in,delta_X_in,kappa_X_in,delta_in,kappa_in,nu_in,B0)
% generates an analytic GEQDSK equilibrium file based on the method of
% Guazzotto and Freidberg (JPP 2021).  Input arguments:
% eq_option_in = 1 for elongated limiter equilibrium, 2 for double null, 3
% for single null (lower)
% R0_in = major radius of plasma center, eps_in = inverse aspect ratio,
% delta_X_in = triangularity of plasma near x-point (eq_options 2 and 3)
% kappa_X_in = elongation of plasma near x-point (eq_options 2 and 3)
% delta_in = triangularity of plasma away from x-point (eq_options 1, 2,
% and 3)
% kappa_in = elongation of plasma away from x-point (eq_options 1, 2, and
% 3)
% nu_in = beta_p like parameter (see JOPP 2021 paper for details)
% B0 = magnetic field on axis in Tesla
% The script generates an eqdsk output file with the name 'geqdsk_analytic'
% upon completion.



% This solves all the shapes.

global nseries  hvec kvec eps kappa kappa_X delta_X xtab xtab_m1 nu...
    xtabm xtabp xtabTop xtabX xtabm_m1 xtabp_m1 xtabTop_m1 xtabX_m1 xm xp zTop xX...
    dh powers powers_m1 twoD epsbar yTop yX cfull sfull csi psisign kappa delta...
    kvec_ref eq_option nlevels vlevel f_J f_P psi_surf psi_axis_input psimax...
    pedestal_option alphasol eq_option


pedestal_option = 0;
% eq_option = 3;
% 1 for limiter profile,
% 2 for Double Null
% 3 for Single Null

twoD = 0;
nR = 129;
nZ = 129;

noplots = 0;

eq_option = eq_option_in;
R0 = R0_in;
eps = eps_in;
delta_X = delta_X_in;
kappa_X = kappa_X_in;
delta = delta_in;
kappa = kappa_in;
epsbar = 2*eps/(1+eps^2);
nu = nu_in;




dh = asin(delta);

if((eq_option==2)||(eq_option==3))
    
    sqdX=sqrt(1-delta_X^2);
    csi = sqdX/(kappa_X-sqdX);
    sqcsi = sqrt(1-csi^2);
    kappa_0 = kappa_X/sqcsi;
    
    xx1 = (csi-delta_X)/(1-csi);
    xx2 = (csi+delta_X)/(1-csi);
    
    theta0 = atan(sqcsi/csi);
    
end


a = eps*R0;


nseries = 220;
psisign = 1;


%-------------------------------------------
% Part 0: define the points.
% Observe that once the geometry is defined all terms involving z never
% change, so they can be calculated once and for all outside the error
% function.

% Here we define the plasma reference shape, too. This is needed only for
% plotting.


if (noplots==0)
    
    if(eq_option==1)
        
        theta = linspace(0, 2*pi,300);
        bigR_edge = R0*(1+eps*cos(theta+dh*sin(theta)));
        bigZ_edge = kappa*R0*eps*sin(theta);
        
    elseif(eq_option==2)
        
        thetaL = linspace(pi-theta0,pi+theta0,150);
        thetaR = linspace(-theta0,theta0,150);
        
        bigX_L = xx1+(1+xx1)*cos(thetaL);
        bigY_L = kappa_0*sin(thetaL);
        
        bigX_R = -xx2+(1+xx2)*cos(thetaR);
        bigY_R = kappa_0*sin(thetaR);
        
        bigX=[bigX_L, bigX_R];
        bigY=[bigY_L, bigY_R];
        
        bigR_edge = a*bigX+R0;
        bigZ_edge = a*bigY;
        
    elseif(eq_option==3)
        
        thetaup = linspace(0,pi,150);
        thetaL = linspace(pi,pi+theta0,150);
        thetaR = linspace(2*pi-theta0,2*pi,150);
        
        bigX_up = cos(thetaup+dh*sin(thetaup));
        bigY_up = kappa*sin(thetaup);
        
        bigX_L = xx1+(1+xx1)*cos(thetaL);
        bigY_L = kappa_0*sin(thetaL);
        
        bigX_R = -xx2+(1+xx2)*cos(thetaR);
        bigY_R = kappa_0*sin(thetaR);
        
        bigX=[bigX_up, bigX_L, bigX_R];
        bigY=[bigY_up, bigY_L, bigY_R];
        
        bigR_edge = a*bigX+R0;
        bigZ_edge = a*bigY;
        
    end
    
end

xm = -1;
xp = 1;

if(eq_option==1)
    zTop = -delta-eps/2*(1-delta^2);
    yTop = kappa;
elseif(eq_option==2)
    xX = -delta_X-eps/2*(1-delta_X^2);
    yX = kappa_X;
    zTop = xX;
    yTop = yX;
elseif(eq_option==3)
    zTop = -delta-eps/2*(1-delta^2);
    yTop = kappa;
    xX = -delta_X-eps/2*(1-delta_X^2);
    yX = -kappa_X;
end

powers = 0:nseries;
powers_m1 = [0,0:nseries-1]; %The first 0 is just a place holder.

if(eq_option==1)
    xtab = [xm; xp; zTop].^powers;
    xtab_m1 = powers.*[xm; xp; zTop].^powers_m1;
elseif(eq_option==2)
    xtab = [xm; xp; xX].^powers;
    xtab_m1 = powers.*[xm; xp; xX].^powers_m1;
elseif(eq_option==3)
    xtab = [xm; xp; zTop; xX].^powers;
    xtab_m1 = powers.*[xm; xp; zTop; xX].^powers_m1;
end

% Still deciding the best way: it may work better to have matrices of the
% same points for different values of k.

xtabm = repmat(xtab(1,:),4,1);
xtabp = repmat(xtab(2,:),4,1);

xtabm_m1 = repmat(xtab_m1(1,:),4,1);
xtabp_m1 = repmat(xtab_m1(2,:),4,1);

if(eq_option==1)
    xtabTop = repmat(xtab(3,:),4,1);
    xtabTop_m1 = repmat(xtab_m1(3,:),4,1);
elseif(eq_option==2)
    xtabX = repmat(xtab(3,:),4,1);
    xtabX_m1 = repmat(xtab_m1(3,:),4,1);
elseif(eq_option==3)
    xtabTop = repmat(xtab(3,:),4,1);
    xtabTop_m1 = repmat(xtab_m1(3,:),4,1);
    xtabX = repmat(xtab(4,:),4,1);
    xtabX_m1 = repmat(xtab_m1(4,:),4,1);
end

% This combination works well for most cases. Feel free to experiment!
kvec_ref = [0 1/6 6/7 1];

alphaminplot = 0.9;
alphamaxplot = 2.75;


nplot = 600;

% This is used to get a good guess for the location of the root.

alphaplot = linspace(alphaminplot,alphamaxplot,nplot);
vals=NaN(size(alphaplot));
for i = 1:length(alphaplot)
    vals(i) = Err_smooth(alphaplot(i));
end
% figure(1)
% plot(alphaplot,vals)
dalpha = (alphamaxplot-alphaminplot)/(nplot-1);



% Now find alpha. First find local minima in vals, then use them as
% starting points for the minimum search. Stop when the value of the
% function is less than 10^-8;

locs = islocalmin(vals);
alphastarts=alphaplot(locs);

options = optimset('Display','iter','TolX',1e-8,'TolFun',1e-8);
% options = optimset('TolX',1e-8,'TolFun',1e-8);

alphaerr=1;
ind = 0;

while(alphaerr>10^-8)
    ind = ind+1;
    alphaL = alphastarts(ind)-dalpha;
    alphaR = alphastarts(ind)+dalpha;
    %    alphasol = fminsearch(@Err_general,alphastarts(ind),options)
    alphasol = fminbnd(@Err_smooth,alphaL,alphaR,options);
    alphaerr = Err_smooth(alphasol);
end

% end





Err_smooth(alphasol);

% In this case we need to check the sign of psi_axis.

psimax = 1;
psiref = psi_any_shape(0,0);

if(psiref<0)
    psisign = -1;
end

psiref = psisign*psiref;

if(noplots==0)
    
    % Let's plot the solution, just to visualize what we have done.
    bigR = linspace(R0*(1-eps)-0.05,R0*(1+eps)+0.05,nR);
    if(eq_option==1)
        bigZ = linspace(-a*kappa-0.05,a*kappa+0.05,nZ);
    elseif(eq_option==2)
        bigZ = linspace(-a*kappa_X-0.05,a*kappa_X+0.05,nZ);
    elseif(eq_option==3)
        bigZ = linspace(-a*kappa_X-0.05,a*kappa+0.05,nZ);
    end
    
    z_for_plot = (bigR.^2/R0^2-1-eps^2)/(2*eps);
    y_for_plot = bigZ./a;
    
    twoD = 1;
    psi2D = psi_any_shape(z_for_plot,y_for_plot);
    
    % Keep in mind that psi_any_shape is normalized to 1 in the center. We can
    % change this if we like.
    nlevels = 10;
    vlevel = linspace(0,psiref,nlevels);
    
    figure(1)
    contour(bigR,bigZ,psi2D,50)
    hold on
    contour(bigR,bigZ,psi2D,[0 0],'r','LineWidth',2)
%     hold off
    grid on
    axis equal
    axis tight
    
end

psirz = -psi2D/2/pi;
lcfs = contour(bigR,bigZ,psirz,[-0.005 -0.005],'LineWidth',2);
hold off

%% Begin post-processing

pick_q0 = 1;

safe_small = 5*10.^-3;

twoD = 0;
psimax = 1;

fun_for_max = @(xx) -psi_any_shape(xx(1),xx(2));

xxstart = [0 0];
xxmax = fminsearch(fun_for_max,xxstart);
psimax = abs(psi_any_shape(xxmax(1),xxmax(2)));

zaxis = xxmax(1);
yaxis = xxmax(2);
R_axis = sqrt(2*eps*zaxis + 1 + eps^2) * R0;
Z_axis = a*yaxis;

% Define remaining input parameters (by default, we have already defined R0
% elsewhere).

if(pick_q0==0)
    
    beta0 = 0.0605; % Toroidal beta on the magnetic axis
    
elseif(pick_q0==1)
    
    q0 = 1;
    psi_xx_0 = psi_xx_any_shape(zaxis,yaxis);
    psi_yy_0 = psi_yy_any_shape(zaxis,yaxis);
    xloc = 1+eps^2+2*eps*zaxis;
    
    psiloc = 1;
    beta0 = nu*eps^2*alphasol^2/(q0^2*xloc^2*psi_xx_0*psi_yy_0-...
        (1+eps^2)*(1-nu)*eps^2*alphasol^2);
    
end

mu0  = 4*pi*10^-7;

zmin = -1; %min(bigZ);
zmax = 1; %max(bigZ);

if((eq_option==1)||(eq_option==2))
    % Miller shape and Double Null
    ymin = -yTop;
    ymax = yTop;
elseif(eq_option==3)
    % Singe Null
    ymin = yX; % yX<0 by definition
    ymax = yTop;
end


% Calculate algebraic quantities.

p0 = beta0 * B0^2/(2*mu0); % Plasma pressure on axis

dB_ov_B = beta0*(1+eps^2)/2*((1-nu)/nu); % Plasma diamagnetism delta B / B0

psi_axis = eps*B0*R0^2/alphasol*sqrt(beta0/nu); % Magnetic flux on axis

% Plot normalized p and J_phi

zofRfun = @(R) (R.^2./R0^2 - 1- eps^2)./(2*eps);
yofZfun = @(Z) Z/a;
Jphifunb = @(R,Z) (2*mu0*p0*R.^2 + 2*R0^2*B0^2*dB_ov_B).* ...
    psi_any_shape(zofRfun(R),yofZfun(Z))./psi_axis./R;

bigR1D = linspace(R0*(1-eps),R0*(1+eps),nR);
z_for_plot_1D = (bigR1D.^2/R0^2-1-eps^2)/(2*eps);
psi1D = psi_any_shape(z_for_plot_1D,zeros(size(bigR1D)));

p1D = psi1D.^2;
p1D = p1D./max(p1D);
Jphi1D = Jphifunb(bigR1D,zeros(size(bigR1D)))./mu0;
Jphi1D = abs(Jphi1D)./max(abs(Jphi1D));

% figure
% plot(bigR1D,p1D,'LineWidth',1.5)
% hold on
% plot(bigR1D,Jphi1D,'LineWidth',1.5)
% legend({'Pressure','J_\phi'},'FontSize',14)
% legend('Location','northwest')
% axis tight





% Now calculate the various integrals that enter into all other quantities.

% First define any needed functions and function handles

% This is used for the area.
fun_Jac2D = @(z) a^2./(sqrt(1+eps^2+2*eps*z));

% We are going to need the plasma shape for various reasons, so we might as
% well determine it here and use it to calculate the surface and line integrals.

ntheta0 = 500; % Number of points on the surface
ntemp = 200;

% Define the input plasma shape.
% Use the desired plasma shape as a starting guess for the outermost
% surface.
if(eq_option==1)
    
    zsfun = @(theta) cos(theta+dh*sin(theta)) - ...
        eps/2 * sin(theta+dh*sin(theta)).^2;
    ysfun = @(theta) kappa*sin(theta);
    
    % Create two y(z) functions that we'll need later.
    % First, positive y.
    
    thetatemp = linspace(0, pi, ntemp);
    ztop_list = zsfun(thetatemp);
    ytop_list = ysfun(thetatemp);
    % Negative y is very similar.
    thetatemp = linspace(pi, 2*pi, ntemp);
    zbot_list = zsfun(thetatemp);
    ybot_list = ysfun(thetatemp);
    
elseif(eq_option==2)
    
    thetaL = linspace(pi-theta0,pi+theta0,ntemp);
    thetaR = linspace(-theta0,theta0,ntemp);
    
    bigX_L = xx1+(1+xx1)*cos(thetaL);
    bigY_L = kappa_0*sin(thetaL);
    
    bigX_R = -xx2+(1+xx2)*cos(thetaR);
    bigY_R = kappa_0*sin(thetaR);
    
    bigX=[bigX_L, bigX_R];
    bigY=[bigY_L, bigY_R];
    
    bigR_edge = a*bigX+R0;
    bigZ_edge = a*bigY;
    
    z_edge = (bigR_edge.^2/R0^2-1-eps^2)/(2*eps);
    y_edge = bigY;
    
    [zboundary_max izmax]=max(z_edge);
    [zboundary_min izmin]=min(z_edge);
    
    zbot_list_step = z_edge(izmin+1:izmax);
    ybot_list_step = bigY(izmin+1:izmax);
    
    [zbot_list isort0] = unique(zbot_list_step);
    ybot_list = ybot_list_step(isort0);
    
    ztop_list_step = [z_edge(1:izmin-1), z_edge(izmax+1:end-1)];
    ytop_list_step = [y_edge(1:izmin-1), y_edge(izmax+1:end-1)];
    
    [ztop_list_step2 isort] = sort(ztop_list_step);
    ytop_list_step2 = ytop_list_step(isort);
    
    [ztop_list isort2] = unique(ztop_list_step2);
    ytop_list = ytop_list_step2(isort2);
    
    
elseif(eq_option==3)
    
    thetaup = linspace(0,pi,ntemp);
    thetaL = linspace(pi,pi+theta0,ntemp);
    thetaR = linspace(2*pi-theta0,2*pi,ntemp);
    
    bigX_up = cos(thetaup+dh*sin(thetaup));
    bigY_up = kappa*sin(thetaup);
    
    bigX_L = xx1+(1+xx1)*cos(thetaL);
    bigY_L = kappa_0*sin(thetaL);
    
    bigX_R = -xx2+(1+xx2)*cos(thetaR);
    bigY_R = kappa_0*sin(thetaR);
    
    bigX=[bigX_up, bigX_L, bigX_R];
    bigY=[bigY_up, bigY_L, bigY_R];
    
    bigR_edge = a*bigX+R0;
    bigZ_edge = a*bigY;
    
    z_edge = (bigR_edge.^2/R0^2-1-eps^2)/(2*eps);
    y_edge = bigY;
    
    ztop_list = z_edge(1:ntemp);
    ytop_list = bigY(1:ntemp);
    
    zbot_list = z_edge(ntemp+1:end);
    ybot_list = bigY(ntemp+1:end);
    
end

fysofzT = spline(ztop_list,ytop_list);
funysofzT = @(z) ppval(fysofzT,z);

fysofzB = spline(zbot_list,ybot_list);
funysofzB = @(z) ppval(fysofzB,z);

% Create the arrays.
npoints = nR;
zlist = NaN(npoints,1);
yuplist = zlist;
ydownlist = zlist;

% Set up known edge points.
zlist = linspace(-1,1,npoints);
yuplist(1) = 0;
yuplist(end) = 0;
ydownlist(1) = 0;
ydownlist(end) = 0;

% Work on the top shape first.
for i = 2:npoints-1
    
    zloc = zlist(i);
    ffory = @(y) psi_any_shape(zloc,y);
    zerostartb = yuplist(i-1);
    zerostart = [0 zerostartb];
    while(ffory(zerostart(1,1))*ffory(zerostart(1,2))>0)
        temp = zerostart(1,2) + 0.002;
        zerostart(1,2) = temp;
        
    end
    yloc = fzero(ffory,zerostart);
    yuplist(i) = yloc;
    
end

% Repeat the calculation for the bottom shape.
for i = 2:npoints-1
    
    zloc = zlist(i);
    ffory = @(y) psi_any_shape(zloc,y);
    zerostartb = ydownlist(i-1);
    zerostart = [0 zerostartb];
    while(ffory(zerostart(1,1))*ffory(zerostart(1,2))>0)
        temp = zerostart(1,2) - 0.002;
        zerostart(1,2) = temp;
        
    end
    yloc = fzero(ffory,zerostart);
    ydownlist(i) = yloc;
    
end

zlist_for_FLOW = [flipud(zlist'); zlist(2:end)'];
ylist_for_FLOW = [flipud(yuplist); ydownlist(2:end)];

R_for_FLOW = R0 * sqrt(1 + eps^2 + 2*eps*zlist_for_FLOW);
Z_for_FLOW = a*ylist_for_FLOW;

save_for_FLOW = [R_for_FLOW Z_for_FLOW];

save('R_Z_shape.dat','save_for_FLOW','-ascii');

fysofzT = spline(zlist,yuplist);
funysofzT = @(z) ppval(fysofzT,z);

fysofzB = spline(zlist,ydownlist);
funysofzB = @(z) ppval(fysofzB,z);

ytoplast = @(z) funysofzT(z);
ybotlast = @(z) funysofzB(z);

% Notice that we calculate the same integrals twice using different
% routines and methods for debugging purposes. This is overkill for this
% problem. The older version can be commented out.

fun_A_new = @(z,y) fun_Jac2D(z);
fun_A_new = @(z,y) 1+z-z;

fun_psi2_new = @(z,y) psi_any_shape(z,y).^2;
fun_psi2_new = @(z,y) psi_for_integral(z,y).^2;

fun_t1_new = @(z,y) (1+nu*epsbar*z)./(1+epsbar*z).* psi_for_integral(z,y);

fun_t2_new = @(z,y) (1+nu*epsbar*z)./(1+epsbar*z).* fun_psi2_new(z,y);

% Now get the integrals

int_A_new = integral2(fun_A_new,zmin,zmax,ybotlast,ytoplast);

int_psi2_new = integral2(fun_psi2_new,zmin,zmax,ybotlast,ytoplast);

int_t2_new = integral2(fun_t2_new,zmin,zmax,ybotlast,ytoplast);

int_t1_new = integral2(fun_t1_new,zmin,zmax,ybotlast,ytoplast);


fun_A = @(z,y) fun_Jac2D(z).*psi_for_integral(z,y)./(psi_for_integral(z,y)+10^-16);

fun_psi2 = @(z,y) psi_for_integral(z,y).^2;

fun_t1 = @(z,y) (1+nu*epsbar*z)./(1+epsbar*z).* psi_for_integral(z,y);

fun_t2 = @(z,y) (1+nu*epsbar*z)./(1+epsbar*z).* fun_psi2(z,y);

% Now get the integrals

int_A = integral2(fun_A,zmin,zmax,ymin,ymax,'AbsTol',1e-2,'RelTol',1e-4);

int_psi2 = integral2(fun_psi2,zmin,zmax,ymin,ymax);

int_t2 = integral2(fun_t2,zmin,zmax,ymin,ymax);

int_t1 = integral2(fun_t1,zmin,zmax,ymin,ymax);

% Now we can calculate the desired integral quantities.

% betator = beta0 * int_psi2/int_A;
betator_new = beta0 * int_psi2_new/int_A_new;

% betapol = nu * int_psi2/int_t2;
betapol_new = nu * int_psi2_new/int_t2_new;

I_plasma = eps*B0*R0*alphasol/mu0 * sqrt(beta0/nu) * int_t1;
I_plasma_new = eps*B0*R0*alphasol/mu0 * sqrt(beta0/nu) * int_t1_new;

l_inductance = 4*pi/alphasol^2 * int_t2/int_t1^2;
l_inductance_new = 4*pi/alphasol^2 * int_t2_new/int_t1_new^2;


% We have moved q_star to the end because it requires kappa_95 for the
% cases with X points

% We need to explicitly define F(psi) for the FLOW version.
bigF = @(psi) R0*B0 *sqrt(1 + 2*dB_ov_B *psi.^2);

% q_star_2 = 2*int_A*bigF(0)/(mu0*R0^2*I_plasma);
% q_star_2_new = 2*int_A_new*bigF(0)/(mu0*R0^2*I_plasma);


% Now we move to line integrals. These only matter for q(psi).
% It is easier to work on one contour at a time, so add here any other line
% integrals that may come up.

% We use Jeff's simplification of the integrand to
% ymax(psi)*sin(theta)/psi_z(z(theta),y(theta))

% Set a psi value for each surface. Use the psi value to determine the z
% values for the y=0 axis first, then set up a grid in z and determine y.

enq = 10;
ntheta = 401*ones(enq,1); % Edit into a non-constant array if needed. Any odd number will work.
qtab_smooth = NaN(enq+1,2); % Includes the axis

psi_for_q_list = NaN(enq,1);
psi_for_q_list = linspace(0,enq/(enq+1),enq);
%psi_for_q_list(end) = 0.999; % Uncomment this line to confirm q0.

% Make sure that one value is 0.05, so that we can give q95
[temppsi, i95] = min(abs(psi_for_q_list(2:end)-0.05));
psi_for_q_list(i95+1) = 0.05;

zlastL = -1;
zlastR = 1;

% The safety factor goes to infinity on the LCFS if there is an X point.
% For that reason, we start from the second value of psi.
% If so desired, one can run the integrals starting from the psi=0 surface
% for all shapes. The code will return a "large" value if there is an
% X-point.
if(eq_option==1)
    iqstart = 1;
else
    iqstart = 2;
end


for i = iqstart:enq    
    % Start from the outside and move inwards.
    
    psiloc = psi_for_q_list(i);
    nthetaloc = ntheta(i);
    nthetamid = (nthetaloc-1)/2;
    
    % Note that this overwrites the previous lists.
    zlist = NaN(nthetaloc,1);
    ylist = NaN(nthetaloc,1);
    % Find the left and right crossings of the horizontal axis for this
    % value of psi.
    
    if(i==1)
        % We set up the arrays in a different way following earlier versions.
        % The easiest thing to do is to redefine the arrays in the same way as
        % for all other surfaces.
        zcrossL = -1;
        zcrossR = 1;
        
        zlastL = zcrossL;
        zlastR = zcrossR;
        
        % Fill the z grid
        zlist(1:nthetamid) = linspace(zcrossR,zcrossL,nthetamid);
        zlist(nthetamid:end) = linspace(zcrossL,zcrossR,nthetamid+2);
        
        % Fill the y grid using the functions created earlier.
        % (Repeated point is just to agree with the rest, it does not create any
        % issues)
        ylist(1:nthetamid) = ytoplast(zlist(1:nthetamid));
        ylist(nthetamid:end) = ybotlast(zlist(nthetamid:end));
        
    else
        
        fforcross = @(z) psi_any_shape(z,0)-psiloc;
        crossstart = [zlastL zaxis];
        zcrossL = fzero(fforcross,crossstart);
        crossstart = [ zaxis zlastR];
        zcrossR = fzero(fforcross,crossstart);
        
        zlastL = zcrossL;
        zlastR = zcrossR;
        
        % Fill the z grid
        zlist(1:nthetamid) = linspace(zcrossR,zcrossL,nthetamid);
        zlist(nthetamid:end) = linspace(zcrossL,zcrossR,nthetamid+2);
        
        % Now fill the y grid. We need to do this by steps.
        
        ylist(1) = 0;
        
        for j = 2:nthetamid-1
            
            zloc = zlist(j);
            ffory = @(y) psi_any_shape(zloc,y)-psiloc;
            zerostartb = ytoplast(zloc);
            zerostart = [0 zerostartb];
            jj = 0;
            while(ffory(zerostart(1,1))*ffory(zerostart(1,2))>0)
                temp = zerostart(1,2) + 0.002;
                zerostart(1,2) = temp;
            end
            yloc = fzero(ffory,zerostart);
            ylist(j) = yloc;
            
        end
        
        ylist(nthetamid) = 0;
        
        for j = nthetamid+1:nthetaloc-1
            
            zloc = zlist(j);
            ffory = @(y) psi_any_shape(zloc,y)-psiloc;
            zerostart = [ybotlast(zloc) 0];
            while(ffory(zerostart(1))*ffory(zerostart(2))>0)
                zerostart(1) = zerostart(1)-0.002;
            end
            yloc = fzero(ffory,zerostart);
            ylist(j) = yloc;
            
        end
        
        ylist(nthetaloc) = 0;
        
        % Now redefine the y(z) functions and use them to find ymax(s) for this
        % surface.
        
        fystepofzB = spline(zlist(1:nthetamid),ylist(1:nthetamid));
        ytoplast = @(z) ppval(fystepofzB,z);
        
        fystepofzT = spline(zlist(nthetamid:end),ylist(nthetamid:end));
        ybotlast = @(z) ppval(fystepofzT,z);
        
    end
    yuplist = ylist(1:nthetamid);
    ydownlist = ylist(nthetamid+1:end);
    
    [zbot, ybot] = fminbnd(ybotlast,zlastL,zlastR);
    ymaxB = -ybot;
    
    ytemp = @(z) -ytoplast(z);
    [ztop, ytop] = fminbnd(ytemp,zlastL,zlastR);
    ymaxT = -ytop;
    
    if(psiloc==0.05)
        kappa95 = (ymaxB+ymaxT)/2;
    end
    
    % Now determine theta
    
    thetatempT = real(asin(yuplist./ymaxT));
    thetatempB = real(asin(ydownlist./ymaxB));
    
    [thmaxT iTop]=max(thetatempT);
    thetatempT(iTop+1:end) = pi - thetatempT(iTop+1:end);
    
    [thminB iBot]=min(thetatempB);
    thetatempB(1:iBot) = -thetatempB(1:iBot)+pi;
    thetatempB(iBot+1:end) = 2*pi + thetatempB(iBot+1:end);
    
    theta_list = [thetatempT; thetatempB];
    
    
    spline_z_level = spline(theta_list,zlist);
    zlvlfun = @(x) ppval(spline_z_level,x);
    
    % We need two integrands because the plasma shape may not be symmetric.
    fq_smoothT = @(th) -ymaxT*cos(th)./((1+eps^2+2*eps*zlvlfun(th)) .* ...
        psi_x_any_shape(zlvlfun(th),ymaxT*sin(th)));
    fq_smoothB = @(th) -ymaxB*cos(th)./((1+eps^2+2*eps*zlvlfun(th)) .* ...
        psi_x_any_shape(zlvlfun(th),ymaxB*sin(th)));
    
    % Now, the issue is that this function is analytically smooth, but
    % numerically it can blow up near the zeros of the denominator. As of
    % now, it is unclear why this happens, but we'll deal with it as we can
    % for a start.
    
    % First, find the zeros.
    
    fden = @(th) psi_x_any_shape(zlvlfun(th),ymaxT*sin(th));
    zero1 = fzero(fden,pi/2);
    fden = @(th) psi_x_any_shape(zlvlfun(th),ymaxB*sin(th));
    zero2 = fzero(fden,3*pi/2);
    
    % Then break the integral in four parts, just "dropping" the trouble
    % spots.
    
    % Part 1
    
    th1_end = zero1-safe_small;
    int_part1 = integral(fq_smoothT,0,th1_end);
    th_interval_1 = th1_end;
    th_interval_1_right = zero1;
    
    int_part1 = int_part1/th_interval_1*th_interval_1_right;
    
    % Part 2
    th2_start = zero1+safe_small;
    th2_end = pi;
    int_part2 = integral(fq_smoothT,th2_start,th2_end);
    th_interval_2 = th2_end-th2_start;
    th_interval_2_right = pi-zero1;
    
    int_part2 = int_part2/th_interval_2*th_interval_2_right;
    
    % Part 3
    th3_start = pi;
    th3_end = zero2-safe_small;
    int_part3 = integral(fq_smoothB,th3_start,th3_end);
    th_interval_3 = th3_end-th3_start;
    th_interval_3_right = th_interval_3+safe_small;
    
    int_part3 = int_part3/th_interval_3*th_interval_3_right;
    
    
    % Part 4
    th4_start = zero2+safe_small;
    th4_end = 2*pi;
    int_part4 = integral(fq_smoothB,th4_start,th4_end);
    th_interval_4 = th4_end-th4_start;
    th_interval_4_right = th_interval_4+safe_small;
    
    int_part4 = int_part4/th_interval_4*th_interval_4_right;
    
    inte_smooth = int_part1+int_part2+int_part3+int_part4;
    
    qloc_smooth = inte_smooth * bigF(psiloc)/(R0*B0)*sqrt(nu/beta0)* ...
        eps*alphasol/(2*pi);
    
    qtab_smooth(i,:) = [psiloc qloc_smooth];
    if(psiloc==0)
        qa = qloc_smooth;
    elseif(psiloc==0.05)
        q95 = qloc_smooth;
    end
    
end



% q_axis calculation

if(pick_q0==0)
    
    psi_xx_0 = psi_xx_any_shape(zaxis,yaxis);
    psi_yy_0 = psi_yy_any_shape(zaxis,yaxis);
    xloc = 1+eps^2+2*eps*zaxis;
    
    psiloc = 1;
    q0 = bigF(psiloc)/(R0*B0) * sqrt(nu/beta0) * eps * alphasol / ...
        (xloc * sqrt(psi_xx_0*psi_yy_0));
    
else
    
    psiloc=1;
    q0;
    
end

qtab_smooth(end,:) = [psiloc q0];

% End of q calculation
% Keep in mind that we are only using 50 points for the level curves, but
% one is free to use as many as desired. Since the solution is analytical,
% there is no intrinsic limitation on the resolution!

% Now also calculate q_star
if(eq_option==1)
    q_star = pi*eps/alphasol * sqrt(nu/beta0) * (1+kappa^2)/int_t1;
    q_star_new = pi*eps/alphasol * sqrt(nu/beta0) * (1+kappa^2)/int_t1_new;
elseif(eq_option==2)
    q_star = pi*eps/alphasol * sqrt(nu/beta0) * (1+kappa95^2)/int_t1;
    q_star_new = pi*eps/alphasol * sqrt(nu/beta0) * (1+kappa95^2)/int_t1_new;
elseif(eq_option==3)
    q_star = pi*eps/alphasol * sqrt(nu/beta0) * (1+kappa95^2)/int_t1;
    q_star_new = pi*eps/alphasol * sqrt(nu/beta0) * (1+kappa95^2)/int_t1_new;
end

% Save some more FLOW input.

psilist = linspace(0,1,nR);
plist = beta0*B0^2/(2*mu0)*psilist.^2;
Flist = bigF(psilist);

pofpsi_out = [psilist; plist]';
Fofpsi_out = [psilist; Flist]';

% save('p_iso.dat','pofpsi_out','-ascii');
% save('b0.dat','Fofpsi_out','-ascii');
% 
% save('R_Z_shape.dat','save_for_FLOW','-ascii');
% 
% if(noplots==0)
%     figure
%     plot(qtab_smooth(iqstart:end,1),qtab_smooth(iqstart:end,2),'-o')
%     hold on
% end

[zmag_ind rmag_ind] = find(psirz==min(min(psirz)));
rmag = bigR(rmag_ind);
zmag = bigZ(zmag_ind);

%% Write results to eqdsk file

filename = ['geqdsk_analytic'];

nw = nR;
nh = nZ;
rdim = abs(bigR(end)-bigR(1));
zdim = abs(bigZ(end)-bigZ(1));
rleft = bigR(1);
rcentr = rleft+rdim/2;
zmid = 0;
rmaxis = rmag;
zmaxis = zmag;
simag = min(min(psirz));
sibry = 0;
bcentr = B0;
current = I_plasma_new;
xdum = 0;

fpol = fliplr(Flist);
pres = fliplr(plist);
ffprim =fpol.*diff(interp1(fpol,1:(length(fpol)+1)))*length(fpol);
pprime =diff(interp1(pres,1:(length(pres)+1)))*length(pres);
qpsi = interp1(-qtab_smooth(2:end,1)'/pi,qtab_smooth(2:end,2)',linspace(simag,sibry,nR),'spline');
nbbbs = ceil((size(lcfs,2)-1)/2);
limitr = [5];
bbbs(1:2:2*nbbbs) = lcfs(1,2:2:end);
bbbs(2:2:2*nbbbs) = lcfs(2,2:2:end);
lim = [min(bigR) max(bigZ) max(bigR) max(bigZ) max(bigR) min(bigZ) min(bigR) min(bigZ) min(bigR) max(bigZ)];


dlmwrite(filename,[1 nw nh],'delimiter',' ','coffset',51)
dlmwrite(filename,[rdim zdim rcentr rleft zmid],'-append','precision','%16.9e','delimiter','')
dlmwrite(filename,[rmaxis zmaxis simag sibry bcentr],'-append','precision','%16.9e','delimiter','')
dlmwrite(filename,[current simag xdum rmaxis xdum],'-append','precision','%16.9e','delimiter','')
dlmwrite(filename,[zmaxis xdum sibry xdum xdum],'-append','precision','%16.9e','delimiter','')
for ind=1:floor(nw/5)
    dlmwrite(filename,fpol(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,fpol(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
for ind=1:floor(nw/5)
    dlmwrite(filename,pres(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,pres(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
for ind=1:floor(nw/5)
    dlmwrite(filename,ffprim(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,ffprim(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
for ind=1:floor(nw/5)
    dlmwrite(filename,pprime(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,pprime(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
psirz_reshape=reshape(psirz',1,nw*nh);
for ind=1:floor(nw*nh/5)
    dlmwrite(filename,psirz_reshape(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,psirz_reshape(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
for ind=1:floor(nw/5)
    dlmwrite(filename,qpsi(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,qpsi(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
dlmwrite(filename,[nbbbs limitr],'-append','precision','  %i','delimiter','')
for ind=1:floor(nbbbs*2/5)
    dlmwrite(filename,bbbs(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,bbbs(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')
for ind=1:floor(limitr*2/5)
    dlmwrite(filename,lim(((ind-1)*5+1):((ind-1)*5+5)),'-append','precision','%16.9e','delimiter','')
end
dlmwrite(filename,lim(((ind)*5+1):end),'-append','precision','%16.9e','delimiter','')

disp([filename ' saved'])
