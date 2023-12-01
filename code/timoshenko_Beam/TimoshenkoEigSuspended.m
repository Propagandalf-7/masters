function [u,p,Eig] = TimoshenkoEigSuspended()
format long g;
syms A;
syms B;
syms C;
syms D;
syms x;
syms m;
syms o;
syms lam;
syms k;
syms a;
syms t;

%ell = 0.04004;
%d = 0.01004;
%b = 0.01009;
%mass = 0.01143;
%rho = mass/(ell*d*b);
%rho = 2817.9111;
%E = 72666110000;
%G = 27174810000;
%nu = E/(2*G)-1;
%kap = 5*(1+nu)/(6+5*nu);
%A = ell*d;
%I = ell*d^3/12;

%k_i = 0;
%gamma = (G*kap^2)/E;
%alpha = 12/((d/ell))^2;
%T = ell*sqrt(rho/(G*kap^2));

%syms w;
%a = rho/E*(A/I-rho*w^2/(kap^2*G));
%b = (rho/E*(1+E/(kap^2*G)))^(1/2);
%lam1 = ((b^4*w^4/4+a*w^2)^(1/2)-b^2*w^2/2)^(1/2);
%lam2 = ((b^4*w^4/4+a*w^2)^(1/2)+b^2*w^2/2)^(1/2);
%a1 = (rho*w^2/(kap^2*G)+lam1^2);
%a2 = (rho*w^2/(kap^2*G)-lam2^2);

delt = 4*gamma/(1+gamma)^2*alpha/lam + (1-gamma)^2/(1+gamma)^2;
omega2 = 1/2*lam*(1+gamma)*(delt^(1/2)+1);
mu2 = simplify(1/2*lam*(1+gamma)*(delt^(1/2)-1));

u = A*sinh(m*x) + B*cosh(m*x) + C*sin(o*x) + D*cos(o*x);
p = A*((lam+m^2)/m*cosh(m*x)) + B*((lam+m^2)/m*sinh(m*x)) + C*(-(lam-o^2)/o*cos(o*x)) + D*((lam-o^2)/o*sin(o*x));

%u = B + C*sin(o*x) + D*cos(o*x);
%p = A + B*a*x + C*(-(l-o^2)/o)*cos(o*x) + D*((l-o^2)/o)*sin(o*x);

%u = A*sin(t*x) + B*cos(t*x) + C*sin(o*x) + D*cos(o*x);
%p = A*(-(l-t^2)/t)*cos(t*x) + B*(l-t^2)/t*sin(t*x) + C*(-(l-o^2)/o)*cos(o*x) + D*(l-o^2)/o*sin(o*x);

subs(diff(p),x,0);
(subs(diff(u) - p - k*u,x,0));

u = subs(u,[D,C],[(m^2+lam)/(o^2-lam)*B,o/m*A + k*(B + (m^2+lam)/(o^2-lam)*B)]);
p = subs(p,[D,C],[(m^2+lam)/(o^2-lam)*B,o/m*A + k*(B + (m^2+lam)/(o^2-lam)*B)]);

%u = subs(u,[B,C],[-D*((l-o^2)/a),o/l*(A + k*(-D*((l-o^2)/a) +D))]);
%p = subs(p,[B,C],[-D*((l-o^2)/a),o/l*(A + k*(-D*((l-o^2)/a) +D))]);

%u = subs(u,[B,C],[-D*((l-o^2)/(l-t^2)),o/l*(-A*l/t + k*(-D*((l-o^2)/(l-t^2)) +D))]);
%p = subs(p,[B,C],[-D*((l-o^2)/(l-t^2)),o/l*(-A*l/t + k*(-D*((l-o^2)/(l-t^2)) +D))]);

M1  = (subs(diff(p),x,1));
M2 = (subs(diff(u) - p +k*u,x,1));
M = [subs(M1,[A,B],[1,0]) subs(M1,[A,B],[0,1]);subs(M2,[A,B],[1,0]) subs(M2,[A,B],[0,1])];
 
L = subs(M,k,k_i);
L = subs(L,o,(omega2)^(1/2));
L = subs(L,m,(mu2)^(1/2));

Y = simplify(det(L));
Y = subs(Y,lam,x);
Y = simplify(Y);

%syms w
%A = d*ell;
%I = ell*d^3/12;
%b2 = rho/E*(1+E/(kap^2*G));
%a = rho/E*(A/I-rho*w^2/(kap^2*G));
%lam1 = (b2^2*w^4/4+a*w^2) - b2*w^2/2;
%lam2 = (b2^2*w^4/4+a*w^2) + b2*w^2/2;
%a1 = (rho*w^2/(kap^2*G)+lam1^2);
%a2 = (rho*w^2/(kap^2*G)-lam1^2);

%lam1 = m;
%lam2 = o;
%a1 = (lam+m^2);
%a2 = (lam-o^2);
%ell =1;
Y2 = 2+((lam1*a1)/(lam2*a2)-(lam2*a2)/(lam1*a1))*sinh(lam1*ell)*sin(lam2*ell)-2*cosh(lam1*ell)*cos(lam2*ell);
%Y2 = subs(Y2,o,(omega2)^(1/2));
%Y2 = subs(Y2,m,(mu2)^(1/2));
%Y2 = subs(Y2,lam,x)
Y2 = subs(Y2,w,x)

clf
%ezplot(Y,[0,200])
hold on
ezplot(Y2,[0,1000000])
%ezplot(Y/T,[0,200])
grid on

R = FindRoots(Y2,0.0001,1000000,1000)

%R = FindRoots(Y,100,200,0.1)
RF = R;
for i = 1:size(RF,2)
    RF(i) = RF(i)/(2*pi);
end
RF
%RF = zeros(1,size(R,2));
%for i = 1:size(R,2)
%    if(R(i)-0.1>0)
%        RF(i) = FindRoots(Y,R(i)-0.1,R(i)+0.1,0.0001)
%    else
%        RF(i) = FindRoots(Y,0.0001,R(i)+0.1,0.0001)
%    end
%end

RFreq = zeros(1,size(RF,2));
for i = 1:size(RF,2)
    RFreq(i) = (sqrt(RF(i))/(2*pi))/T;
end
Eig = RFreq'

REAL = [27417.5;
        60851.1;
        97796.0];
ERROR = (REAL-Eig(1:3))/REAL*100
%ModeNum = 1;

%{

L = subs(L,lam,RF(ModeNum));
[a,L1] = gauss(L,[0;0]);
A1 = -L1(1,2)/L1(1,1);
%ANum = [-0.954152205207778 -0.930891281851854 -0.872827099855040];
u = subs(u,[o,m],[(omega2)^(1/2),(mu2)^(1/2)]);
p = subs(p,[o,m],[(omega2)^(1/2),(mu2)^(1/2)]);
u = simplify(subs(u,[lam,A,B,k],[RF(ModeNum),A1,1,sqrt(5/6)]));
p = simplify(subs(p,[lam,A,B,k],[RF(ModeNum),A1,1,sqrt(5/6)]));
%u2 = simplify(subs(u,[l,A,D,k],[Roots(ModeNum),ANum(ModeNum),1,0.0000135]));
%u3 = simplify(subs(u,[l,A,D,k],[Roots(ModeNum),ANum(ModeNum),1,0.0000135]));
clf
normMax = norm(u,Inf);
%u = 1/normMax*u;
ezplot(u,[0,1])
%axis([0 1.2 -1.5 1.5])
%}
return;

function I = IntervalDivision(a,b,TOL)
if abs(b-a)>=TOL
   m = abs(b-a)/2;
   I = [IntervalDivision(a,a+m,TOL); IntervalDivision(a+m,b,TOL)];
   return;
else
    I = [a b];
end
return

function R = FindRoots(Y,a,b,TOL)
syms x;
I = IntervalDivision(a,b,TOL);
n = size(I,1);
SubsI = zeros(size(I));
for i = 1:n
    SubsI(i,1) = subs(Y,x,I(i,1));
    SubsI(i,2) = subs(Y,x,I(i,2));
end

icount = 1;
for i = 1:n
   if SubsI(i,1) == 0
       R(icount) = I(i,1);
       icount = icount+1;
   elseif SubsI(i,1) == 0
       R(icount) = I(i,1);
       icount = icount+1;
   elseif SubsI(i,1)*SubsI(i,2) < 0
       R(icount) = (I(i,2)+I(i,1))/2;
       icount = icount+1;
   end
end
return
