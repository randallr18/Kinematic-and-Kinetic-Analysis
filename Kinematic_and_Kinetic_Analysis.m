%Plots are commented out; Please uncomment one at a time to analyze if interested;

% first 4-bar mechanism - driver
pto2x = 0.0;
pto2y = 0.0;

rocker = 8.139;
coupler = 20.048;
crank = 8.058;
ground = 20.081;

a = crank;
b = coupler;
c = rocker;
d = ground;

K1 = d/a;
K2 = d/c;
K3 = (a*a-b*b+c*c+d*d)/(2*a*c);
K4 = d/b;
K5 = (c*c-d*d-a*a-b*b)/(2*a*b);

% second 4-bar mechansim - triple-rocker
pto2x1 = 0.0;
pto2y1 = 0.0;

a1 = 28.912;  % crank1;
b2 = 12;  % coupler2;
c3 = 42.747;  % rocker3;
d4 = -21.77;  % ground4;

K11 = d4/a1;
K21 = d4/c3;
K31 = (a1*a1-b2*b2+c3*c3+d4*d4)/(2*a1*c3);
K41 = d4/b2;
K51 = (c3*c3-d4*d4-a1*a1-b2*b2)/(2*a1*b2);

omega2 = 2.0; %rad/sec
omega2not = omega2;
alpha2 = 0.0; %rad/sec^2
phi2 = 0.0; %rad/sec^3
T4 = 1.0; %Newton meters? torque on link 1

%Masses for each link
m1 = .10; %slugs for all mass
m2 = .15;
m3 = .18;
m4 = .1;
m5 = .15;
msnow = .310;
g = 32.2; %ft/s^2


IO1 = m1*a^2/3; %slugs*ft^2 for all I's
I2 = m2*b^2/12;
IO2 = m3*c^2/3;
I4 = m4*b2^2/12;
IO3 = m5*c3^2/3;

%coefficients of O1x O1y Ax Ay Bx By O2x O2y Cx Cy Dx Dy O3x O3y Tinput

coef=zeros(9,9);
coef(1,1) = 1;
coef(1,3) = 1;
coef(2,2) = 1;
coef(2,4) = 1;
coef(3,15) = 1;
coef(4,3) = -1;
coef(4,5) = 1;
coef(5,4) = -1;
coef(5,6) = 1;
coef(7,7) = 1;
coef(7,5) = -1;
coef(7,9) = 1;
coef(8,8) = 1;
coef(8,6) = -1;
coef(8,10) = 1;
coef(10,9) = -1;
coef(10,11) = 1;
coef(11,10) = -1;
coef(11,12) = 1;
coef(13,13) = 1;
coef(13,11) = -1;
coef(14,14) = 1;
coef(14,12) = -1;


figure
hold on
for i = 1:1:628.8

    % Iteration are factored up by a scale of 10 and start at 1 due to the
    % export formulas used at the end of the for loop. i is scaled back down to the
    % proper values needed on line 92.

    %kinematic - position, velocity and acceleration analysis for the
    %4-bar driver
    o2 = i / 10

    A = cos(o2)-K1-K2*cos(o2)+K3;
    B = -2*sin(o2);
    C = K1 - (K2+1)*cos(o2)+K3;
    D = cos(o2)-K1+K4*cos(o2)+K5;
    E = B;
    F = K1+(K4-1)*cos(o2) + K5;

    o3pos = 2*atan(((-E+sqrt(E*E-4*D*F))/(2*D)));
    o3neg = 2*atan(((-E-sqrt(E*E-4*D*F))/(2*D)));

    o4pos =  2*atan(((-B+sqrt(B*B-4*A*C))/(2*A)));
    o4neg =  2*atan(((-B-sqrt(B*B-4*A*C))/(2*A)));

    otrans = abs(o3pos-o4pos);

    ptax = pto2x+a*cos(o2);
    ptay = pto2y+a*sin(o2);

    pto4x = d;
    pto4y = 0.0;

    ptbaxp = ptax+b*cos(o3pos);
    ptbayp = ptay+b*sin(o3pos);
    ptbaxn = ptax+b*cos(o3neg);
    ptbayn = ptay+b*sin(o3neg);

    ptbo4xp = pto4x + c*cos(o4pos);
    ptbo4yp = pto4y + c*sin(o4pos);
    ptbo4xn = pto4x + c*cos(o4neg);
    ptbo4yn = pto4y + c*sin(o4neg);

    x = [pto2x ptax  ptbaxn ptbo4xn pto4x];
    y = [pto2y ptay  ptbayn ptbo4yn pto4y];

    omega3p = a*omega2/b*sin(o4pos-o2)/sin(o3pos-o4pos);
    omega3m = a*omega2/b*sin(o4neg-o2)/sin(o3neg-o4neg);
    omega4p = a*omega2/c*sin(o2-o3pos)/sin(o4pos-o3pos);
    omega4m = a*omega2/c*sin(o2-o3neg)/sin(o4neg-o3neg);

    %plot(o4, omega3m, '.');

    mvp = omega4p/omega2;
    mvm = omega4m/omega2;

    Ap = c*sin(o4pos);
    Am = c*sin(o4neg);
    Bp = b*sin(o3pos);
    Bm = b*sin(o3neg);
    Cp = a*alpha2*sin(o2)+a*omega2^2*cos(o2)+b*omega3p^2*cos(o3pos)+c*omega4p^2*sin(o4pos);
    Cm = a*alpha2*sin(o2)+a*omega2^2*cos(o2)+b*omega3m^2*cos(o3neg)+c*omega4m^2*sin(o4neg);
    Dp = c*cos(o4pos);
    Dm = c*cos(o4neg);
    Ep = b*cos(o3pos);
    Em = b*cos(o3neg);
    Fp = a*alpha2*cos(o2)-a*omega2^2*sin(o2)-b*omega3p^2*sin(o3pos)+c*omega4p^2*sin(o4pos);
    Fm = a*alpha2*cos(o2)-a*omega2^2*sin(o2)-b*omega3m^2*sin(o3neg)+c*omega4m^2*sin(o4neg);

    alpha3p = (Cp*Dp-Ap*Fp)/(Ap*Ep-Bp*Dp);
    alpha3m = (Cm*Dm-Am*Fm)/(Am*Em-Bm*Dm);
    alpha4p = (Cp*Ep-Bp*Fp)/(Ap*Ep-Bp*Dp);
    alpha4m = (Cm*Em-Bm*Fm)/(Am*Em-Bm*Dm);

    %plot(o4, alpha3m, '.');

    %kinematic - position, velocity and acceleration analysis for the
    %4-bar driver

    %Established relationship between theta of the crank for both four bar
    %mechanisms through trigonometry.
    o4 = 1.454 - o4neg

    A1 = cos(o4)-K11-K21*cos(o4)+K31;
    B1 = -2*sin(o4);
    C1 = K11 - (K21+1)*cos(o4)+K31;
    D1 = cos(o4)-K11+K41*cos(o4)+K51;
    E1 = B1;
    F1 = K11+(K41-1)*cos(o4) + K51;

    o3pos1 = 2*atan(((-E1+sqrt(E1*E1-4*D1*F1))/(2*D1)));
    o3neg1 = 2*atan(((-E1-sqrt(E1*E1-4*D1*F1))/(2*D1)));

    o4pos1 =  2*atan(((-B1+sqrt(B1*B1-4*A1*C1))/(2*A1)));
    o4neg1 =  2*atan(((-B1-sqrt(B1*B1-4*A1*C1))/(2*A1)));

    otrans1 = abs(o3pos1-o4pos1);

    ptax1 = pto2x1+a1*cos(o4);
    ptay1 = pto2y1+a1*sin(o4);

    pto4x1 = d4;
    pto4y1 = 0.0;

    ptbaxp1 = ptax1+b2*cos(o3pos1);
    ptbayp1 = ptay1+b2*sin(o3pos1);
    ptbaxn1 = ptax1+b2*cos(o3neg1);
    ptbayn1 = ptay1+b2*sin(o3neg1);

    ptbo4xp1 = pto4x1 + c3*cos(o4pos1);
    ptbo4yp1 = pto4y1 + c3*sin(o4pos1);
    ptbo4xn1 = pto4x1 + c3*cos(o4neg1);
    ptbo4yn1 = pto4y1 + c3*sin(o4neg1);

    x1 = [pto2x1 ptax1  ptbaxp1 ptbo4xp1 pto4x1];
    y1 = [pto2y1 ptay1  ptbayp1 ptbo4yp1 pto4y1];


    omega3p1 = a1*omega2/b2*sin(o4pos1-o4)/sin(o3pos1-o4pos1);
    omega3m1 = a1*omega2/b2*sin(o4neg1-o4)/sin(o3neg1-o4neg1);
    omega4p1 = a1*omega2/c3*sin(o4-o3pos1)/sin(o4pos1-o3pos1);
    omega4m1 = a1*omega2/c3*sin(o4-o3neg1)/sin(o4neg1-o3neg1);

    %plot(o4, omega4p1, '.');

    mvp = omega4p/omega2;
    mvm = omega4m/omega2;

    Ap1 = c3*sin(o4pos1);
    Am1 = c3*sin(o4neg1);
    Bp1 = b2*sin(o3pos1);
    Bm1 = b2*sin(o3neg1);
    Cp1 = a1*alpha2*sin(o4)+a1*omega2^2*cos(o4)+b2*omega3p1^2*cos(o3pos1)+c3*omega4p1^2*sin(o4pos1);
    Cm1 = a1*alpha2*sin(o4)+a1*omega2^2*cos(o4)+b2*omega3m1^2*cos(o3neg1)+c3*omega4m1^2*sin(o4neg1);
    Dp1 = c3*cos(o4pos1);
    Dm1 = c3*cos(o4neg1);
    Ep1 = b2*cos(o3pos1);
    Em1 = b2*cos(o3neg1);
    Fp1 = a1*alpha2*cos(o4)-a1*omega2^2*sin(o4)-b2*omega3p1^2*sin(o3pos1)+c3*omega4p1^2*sin(o4pos1);
    Fm1 = a1*alpha2*cos(o4)-a1*omega2^2*sin(o4)-b2*omega3m1^2*sin(o3neg1)+c3*omega4m1^2*sin(o4neg1);

    alpha3p1 = (Cp1*Dp1-Ap1*Fp1)/(Ap1*Ep1-Bp1*Dp1);
    alpha3m1 = (Cm1*Dm1-Am1*Fm1)/(Am1*Em1-Bm1*Dm1);
    alpha4p1 = (Cp1*Ep1-Bp1*Fp1)/(Ap1*Ep1-Bp1*Dp1);
    alpha4m1 = (Cm1*Em1-Bm1*Fm1)/(Am1*Em1-Bm1*Dm1);

    %plot(o4, alpha4p1, '.');

    %Rotated cordinate system. See notes for details.
    betaDyad = -3.71/180*pi;
    beta2 = o2 + betaDyad;
    beta3m = o3neg + betaDyad;
    beta4m = o4neg + betaDyad;
    betaFour = 57.59/180*pi;
    beta5 = o4 + betaFour;
    beta6p = o3pos1 + betaFour;
    beta7p = o4pos1 + betaFour;

    %coefficients of O1x O1y Ax Ay Bx By O2x O2y Cx Cy Dx Dy O3x O3y Tinput
    coef(3,3) = -a*sin(beta2);
    coef(3,4) = a*cos(beta2);
    coef(6,3) = b/2.0*sin(beta3m);
    coef(6,4) = b/2.0*cos(beta3m);
    coef(6,5) = b/2.0*sin(beta3m);
    coef(6,6) = b/2.0*cos(beta3m);
    coef(9,10) = c*cos(beta5);
    coef(9,9) = c*sin(beta5);
    coef(9,6) = a1*cos(beta4m);
    coef(9,5) = a1*sin(beta4m);
    coef(12,10) = b2/2.0*sin(beta6p);
    coef(12,11) = b2/2.0*cos(beta6p);
    coef(12,12) = b2/2.0*sin(beta6p);
    coef(12,13) = b2/2.0*cos(beta6p);
    coef(15,11) = c3*sin(beta7p);
    coef(15,12) = c3*cos(beta7p);


    %accelerations at the center of mass for each link. Necessary for
    %kinetic analysis. Calculations found through kinematic acceleration
    %analysis.
    a2x = -alpha2/2.0*a*sin(beta2)-omega2^2/2.0*a*cos(beta2);
    a2y = alpha2/2.0*a*cos(beta2)-omega2^2/2.0*a*sin(beta2);

    aAx = -alpha2*2.0*a*sin(beta2)-omega2^2*a*cos(beta2);
    aAy = -alpha2*2.0*a*cos(beta2)-omega2^2*a*sin(beta2);

    a3x = aAx-b/2.0*sin(beta3m)*alpha3m-b/2.0*cos(beta3m)*omega3m^2;
    a3y = aAy+b/2.0*cos(beta3m)*alpha3m-b/2.0*sin(beta3m)*omega3m^2;

    a4x = -a1/2.0*sin(beta5)*alpha4m-a1/2.0*cos(beta5)*omega4p1^2;
    a4y = a1/2.0*cos(beta5)*alpha4m-a1/2.0*sin(beta5)*omega4p1^2;

    aCx = -alpha2*2.0*a1*sin(beta5)-omega2^2*a1*cos(beta5);
    aCy = -alpha2*2.0*a1*cos(beta5)-omega2^2*a1*sin(beta5);

    a5x = aCx-b2/2.0*sin(beta6p)*alpha3p1-b2/2.0*cos(beta6p)*omega3p1^2;
    a5y = aCy+b2/2.0*cos(beta6p)*alpha3p1-b2/2.0*sin(beta6p)*omega3p1^2;

    a6x = -c3/2.0*sin(beta7p)*alpha4p1-c3/2.0*cos(beta7p)*omega4p1^2;
    a6y = c3/2.0*cos(beta7p)*alpha4p1-c3/2.0*sin(beta7p)*omega4p1^2;


    %right hand side of matrix

    rhs = zeros(15,1);
    rhs(1,1) = m1*a2x;
    rhs(2,1) = m1*a2y+m1*g;
    rhs(3,1) = IO1*alpha2+m1*g*a/2.0*cos(beta2);
    rhs(4,1) = m2*a3x;
    rhs(5,1) = m2*a3y+m2*g;
    rhs(6,1) = I2*alpha3m;
    rhs(7,1) = m3*a4x;
    rhs(8,1) = m3*a4y+m3*g;
    rhs(9,1) = IO2*alpha4m+m3*g*a1/2.0*cos(beta5);
    rhs(10,1) = m4*a5x;
    rhs(11,1) = m4*a5y+m4*g+msnow*g;
    rhs(12,1) = I4*alpha3p1+msnow*g*b2/4.0*cos(beta6p);
    rhs(13,1) = m5*a6x;
    rhs(14,1) = m5*a6y+m5*g;
    rhs(15,1) = IO3*alpha4p1+m3*g/2.0*cos(beta7p);

    X = linsolve(coef, rhs);

    %plot of required input torque with set input omega and mass of snow
    plot(o2, X(15), '.');
    xlim([0 20]);
    ylim([-100 100]);


    %Implict forumula used to export data.
    k(i) = pto2x;
    l(i) = ptax;
    m11(i) = ptbaxn;
    n(i) = ptbo4xn;
    o(i) = pto4x;

    p(i) = pto2y;
    q(i) = ptay;
    r(i) = ptbayn;
    s(i) = ptbo4yn;
    t(i) = pto4y;

    k1(i) = o4neg;
    l1(i) = o3neg;

    m22(i) = o2;


    %plots of the position analysis of both mechanims. Comment out lines 86
    %and 87 to view.

    %plot(x,y,'o-','linewidth',3);
    %xlim([-40 40]);
    %ylim([-40 40]);
    %linkdata on

    %plot(x1,y1,'o-','linewidth',3);
    %xlim([-40 40]);
    %ylim([-40 40]);
    %linkdata on

    %omega2 = sqrt(2.0*alpha2*o2+omega2not^2);
end
