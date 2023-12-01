function [u,p,Eig] = TimoshenkoEigFree(alpha)
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
format long;

gamma = 0.2829767118;

delt = 4*gamma/(1+gamma)^2*alpha/lam + (1-gamma)^2/(1+gamma)^2;
omega2 = 1/2*lam*(1+gamma)*(delt^(1/2)+1);
mu2 = simplifyFraction(1/2*lam*(1+gamma)*(delt^(1/2)-1));
theta2 = 1/2*lam*(1+gamma)*(1-delt^(1/2));


%u = A*sinh(m*x) + B*cosh(m*x) + C*sin(o*x) + D*cos(o*x);
%p = A*((lam+m^2)/m*cosh(m*x)) + B*((lam+m^2)/m*sinh(m*x)) + C*(-(lam-o^2)/o*cos(o*x)) + D*((lam-o^2)/o*sin(o*x));

%u = B + C*sin(o*x) + D*cos(o*x);
%p = A + B*a*x + C*(-(lam-o^2)/o)*cos(o*x) + D*((lam-o^2)/o)*sin(o*x);

u = A*sin(t*x) + B*cos(t*x) + C*sin(o*x) + D*cos(o*x);
p = A*(-(lam-t^2)/t)*cos(t*x) + B*(lam-t^2)/t*sin(t*x) + C*(-(lam-o^2)/o)*cos(o*x) + D*(lam-o^2)/o*sin(o*x);


simplify(subs(diff(p),x,0));
simplify(subs(diff(u) - p,x,0));

%u = subs(u,[D,C],[-((m^2+lam)/(-o^2+lam))*B,o/m*A]);
%p = subs(p,[D,C],[-((m^2+lam)/(-o^2+lam))*B,o/m*A]);

%u = subs(u,[B,C],[-D*((lam-o^2)/a),o/lam*(A + k*(-D*((lam-o^2)/a) +D))]);
%p = subs(p,[B,C],[-D*((lam-o^2)/a),o/lam*(A + k*(-D*((lam-o^2)/a) +D))]);

u = subs(u,[D,C],[-B*((t^2-lam)/(o^2-lam)),-A*o/t]);
p = subs(p,[D,C],[-B*((t^2-lam)/(o^2-lam)),-A*o/t]);

M1 = (subs(diff(p),x,1));
M2 = (subs(diff(u) - p,x,1));
M = [subs(M1,[A,B],[1,0]) subs(M1,[A,B],[0,1]);subs(M2,[A,B],[1,0]) subs(M2,[A,B],[0,1])];

latex(simplify(det(M)));

L = subs(M,k,sqrt(5/6));
L = subs(L,o,(omega2)^(1/2));
L = subs(L,m,(mu2)^(1/2));
L = subs(L,t,(theta2)^(1/2));

Y = det(L);
Y = simplify(subs(Y,lam,x));

R = 0;
R = FindRoots(Y,190.001,600,0.1)
%R = FindRoots(Y,100,200,0.1)
RF = zeros(1,size(R,2));
RF2 = zeros(1,size(R,2));
RF3 = zeros(1,size(R,2));
for i = 1:size(R,2)
    if(R(i)-0.1>0)
        RF(i) = FindRoots(Y,R(i)-0.1,R(i)+0.1,0.0001)
    else
        RF(i) = FindRoots(Y,0.0001,R(i)+0.1,0.0001)
    end
end
RF
for i = 1:size(R,2)
    if(RF(i)-0.0001>0)
        RF2(i) = FindRoots(Y,RF(i)-0.0001,RF(i)+0.0001,0.00001)
    else
        RF2(i) = FindRoots(Y,0.0001,RF(i)+0.0001,0.00001)
    end
end
RF2
for i = 1:size(R,2)
    if(RF2(i)-0.00001>0)
        RF3(i) = FindRoots(Y,RF2(i)-0.00001,RF2(i)+0.00001,0.000001)
    else
        RF3(i) = FindRoots(Y,0.00001,RF2(i)+0.00001,0.000001)
    end
end

Eig = RF3';
%ModeNum = 1;
Eig
Nat = sqrt(Eig)/(2*pi)/1.48223276*10^-5

%%{
for i = 1:size(Eig,1)
    f = figure(i+10);
    RF(i)
    LS = subs(L,lam,RF(i));
    [a,L1] = gauss(LS,[0;0]);
    B1 = double(-L1(1,1)/L1(1,2))
    %ANum = [-0.954152205207778 -0.930891281851854 -0.872827099855040];
    us = subs(u,[o,m],[(omega2)^(1/2),(mu2)^(1/2)]);
    ps = subs(p,[o,m],[(omega2)^(1/2),(mu2)^(1/2)]);
    us = simplify(subs(us,[lam,A,B,k],[RF(i),1,B1,sqrt(5/6)]));
    ps = simplify(subs(ps,[lam,A,B,k],[RF(i),1,B1,sqrt(5/6)]));
    %u2 = simplify(subs(u,[l,A,D,k],[Roots(ModeNum),ANum(ModeNum),1,0.0000135]));
    %u3 = simplify(subs(u,[l,A,D,k],[Roots(ModeNum),ANum(ModeNum),1,0.0000135]));
    clf
    %normMax = norm(u,Inf);
    %u = 1/normMax*u;
    xd = 0:0.01:1;
    uss = subs(us,x,xd);
    max = norm(uss,Inf);
    us = us/max;
    ezplot(us,[0,1])
    
    %pss = subs(ps,x,xd);
    %maxp = norm(pss,Inf);
    %ps = ps/maxp;
    %ezplot(ps,[0,1])
    %axis([0 1.2 -1.5 1.5])
end
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