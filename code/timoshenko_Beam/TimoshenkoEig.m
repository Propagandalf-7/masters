function [u,p,Eig] = TimoshenkoEig(alpha)
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

%gamma = 0.25;
nu = 0.3;
gamma = 1/(2*(1+nu))*5/6;

delt = 4*gamma/(1+gamma)^2*alpha/lam + (1-gamma)^2/(1+gamma)^2;
omega2 = 1/2*lam*(1+gamma)*(delt^(1/2)+1);
mu2 = 1/2*lam*(1+gamma)*(delt^(1/2)-1);
theta2 = 1/2*lam*(1+gamma)*(1-delt^(1/2));


u = A*sinh(m*x) + B*cosh(m*x) + C*sin(o*x) + D*cos(o*x);
p = A*((lam+m^2)/m*cosh(m*x)) + B*((lam+m^2)/m*sinh(m*x)) + C*(-(lam-o^2)/o*cos(o*x)) + D*((lam-o^2)/o*sin(o*x));

%u = B + C*sin(o*x) + D*cos(o*x);
%p = A + B*a*x + C*(-(lam-o^2)/o)*cos(o*x) + D*((lam-o^2)/o)*sin(o*x);

%u = A*sin(t*x) + B*cos(t*x) + C*sin(o*x) + D*cos(o*x);
%p = A*(-(lam-t^2)/t)*cos(t*x) + B*(lam-t^2)/t*sin(t*x) + C*(-(lam-o^2)/o)*cos(o*x) + D*(lam-o^2)/o*sin(o*x);


subs((u),x,0);
subs((p),x,0);

u = subs(u,[D,C],[-B,A*(lam+m^2)/m*o/(lam-o^2)]);
p = subs(p,[D,C],[-B,A*(lam+m^2)/m*o/(lam-o^2)]);

%u = subs(u,[B,C],[-D*((lam-o^2)/a),o/lam*(A + k*(-D*((lam-o^2)/a) +D))]);
%p = subs(p,[B,C],[-D*((lam-o^2)/a),o/lam*(A + k*(-D*((lam-o^2)/a) +D))]);

%u = subs(u,[D,C],[-B,-o/(o^2-lam)*(t^2-lam)/t*A]);
%p = subs(p,[D,C],[-B,-o/(o^2-lam)*(t^2-lam)/t*A]);

M1 = (subs(diff(p),x,1));
M2 = (subs(diff(u) - p,x,1));
M = [subs(M1,[A,B],[1,0]) subs(M1,[A,B],[0,1]);subs(M2,[A,B],[1,0]) subs(M2,[A,B],[0,1])];

%latex(simplify(det(M)))

L = subs(M,k,sqrt(5/6));
L = subs(L,o,(omega2)^(1/2));
L = subs(L,m,(mu2)^(1/2));
L = subs(L,t,(theta2)^(1/2));

Y = det(L);
Y = simplify(subs(Y,lam,x));

%ezplot(Y,[0,300])
%grid on

R = 0;
R = FindRoots(Y,0.001,500,0.1)
%R = FindRoots(Y,100,200,0.1)
RF = zeros(1,size(R,2));
RF2 = zeros(1,size(R,2));
RF3 = zeros(1,size(R,2));
RF4 = zeros(1,size(R,2));
for i = 1:size(R,2)
    if(R(i)-0.1>0)
        RF(i) = FindRoots(Y,R(i)-0.1,R(i)+0.1,0.0001);
    else
        RF(i) = FindRoots(Y,0.0001,R(i)+0.1,0.0001);
    end
end
for i = 1:size(R,2)
    if(RF(i)-0.0001>0)
        RF2(i) = FindRoots(Y,RF(i)-0.0001,RF(i)+0.0001,0.00001);
    else
        RF2(i) = FindRoots(Y,0.0001,RF(i)+0.0001,0.00001);
    end
end
for i = 1:size(R,2)
    if(RF2(i)-0.00001>0)
        RF3(i) = FindRoots(Y,RF2(i)-0.00001,RF2(i)+0.00001,0.000001);
    else
        RF3(i) = FindRoots(Y,0.00001,RF2(i)+0.00001,0.000001);
    end
end
%for i = 1:size(R,2)
 %   if(RF3(i)-0.00001>0)
 %       RF4(i) = FindRoots(Y,RF3(i)-0.000001,RF3(i)+0.000001,0.0000001);
 %   else
%        RF4(i) = FindRoots(Y,0.000001,RF3(i)+0.000001,0.0000001);
%    end
%end

Eig = RF3';
%ModeNum = 1;
Eig
%{
imageDir = fullfile(cd, 'images');
if ~exist(imageDir, 'dir')
   mkdir(imageDir);
end

for i = 1:size(Eig,1)
    % Get values
    LS = subs(L,lam,RF(i));
    [a,L1] = gauss(LS,[0;0]);
    B1 = double(-L1(1,1)/L1(1,2))
    us = subs(u,[o,m],[(omega2)^(1/2),(mu2)^(1/2)]);
    ps = subs(p,[o,m],[(omega2)^(1/2),(mu2)^(1/2)]);
    us = simplify(subs(us,[lam,A,B,k],[RF(i),1,B1,sqrt(5/6)]));
    ps = simplify(subs(ps,[lam,A,B,k],[RF(i),1,B1,sqrt(5/6)]));
    
    xd = 0:0.01:1;
    uss = subs(us,x,xd);
    max = norm(uss,Inf);
    us = us/max;
    
    pss = subs(ps,x,xd);
    maxp = norm(pss,Inf);
    ps = ps/maxp;
    
    % Displacement
    f1 = figure('Name', ['Mode ' num2str(i) ' Displacement']);
    clf(f1)
    ezplot(us,[0,1])
    title(['Mode ' num2str(i) ' Displacement'])
    xlabel('x (Position)')
    ylabel('Displacement (Normalized)')
    legend('Displacement', 'Location', 'best')
    
    % Stress
    f2 = figure('Name', ['Mode ' num2str(i) ' Stress']);
    clf(f2)
    ezplot(ps,[0,1])
    title(['Mode ' num2str(i) ' Stress Distribution'])
    xlabel('x (Position)')
    ylabel('Stress (Normalized)')
    legend('Stress Distribution', 'Location', 'best')
    
    % Both
    f3 = figure('Name', ['Mode ' num2str(i) ' Displacement and Stress']);
    clf(f3)
    hold on
    ezplot(us,[0,1])
    ezplot(ps,[0,1])
    title(['Mode ' num2str(i) ' Displacement and Stress'])
    xlabel('x (Position)')
    legend('Displacement', 'Stress Distribution', 'Location', 'best')
    hold off
    
    saveas(f1, fullfile(imageDir, ['Mode_' num2str(i) '_Displacement.png']));
    saveas(f2, fullfile(imageDir, ['Mode_' num2str(i) '_Stress.png']));
    saveas(f3, fullfile(imageDir, ['Mode_' num2str(i) '_Displacement_and_Stress.png']));
end
writeToExcel(Eig, imageDir);
%}
return;

function writeToExcel(Eig, imageDir)
    % Define the name of the Excel file
    excelFileName = 'TimoshenkoResults.xlsx';

    % Initialize COM server
    Excel = actxserver('Excel.Application');
    Excel.Workbooks.Add;

    % Get active sheet
    WorkSheets = Excel.ActiveWorkBook.Sheets;
    sheet1 = WorkSheets.get('Item', 1);
    sheet1.Activate;

    % Start writing data to Excel
    sheet1.Range('A1').Value = 'Mode Number';
    sheet1.Range('B1').Value = 'Eigen Value';

    for i = 1:size(Eig, 1)
        sheet1.Range(['A' num2str(i + 1)]).Value = i; % Mode Number
        sheet1.Range(['B' num2str(i + 1)]).Value = Eig(i); % Eigen Value
        
        % Insert images
        pic_path = fullfile(imageDir, ['Mode_' num2str(i) '_Displacement.png']);
        disp(['Image path: ', pic_path]);
        Excel.ActiveSheet.Shapes.AddPicture(pic_path, 0, 1, 100, i*100, 200, 200);
        pic_path = fullfile(imageDir, ['Mode_' num2str(i) '_Stress.png']);
        disp(['Image path: ', pic_path]);
        Excel.ActiveSheet.Shapes.AddPicture(pic_path, 0, 1, 100, i*100, 200, 200);
        pic_path = fullfile(imageDir, ['Mode_' num2str(i) '_Displacement_and_Stress.png']);
        disp(['Image path: ', pic_path]);
        Excel.ActiveSheet.Shapes.AddPicture(pic_path, 0, 1, 100, i*100, 200, 200);
    end

    % Save and close the Excel file
    Excel.ActiveWorkBook.SaveAs(excelFileName);
    pause(1); % waits for 1 second

    Excel.ActiveWorkbook.Close;
    Excel.Quit;
    Excel.delete;
    clear Excel;

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