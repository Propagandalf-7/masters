%function [E,wP,xP,yP,size_c] = PlateCantileverCubic(d,n,h,inum,numEig)
function [E,n,m] = PlateCantileverCubic(d,n,h,inum,numEig)
format long g
m = ceil(n*d);
a = 0;
b = 1;
c = 0;

%h = sqrt(12/alpha);
%d = 1;
deltx = (b-a)/n;
delty = (d-c)/m;

size_c = (n+1)*(m+1);

nu = 0.3;

kappa_b = (5/6);
kappa_p = 0.9554;%0.29738*nu + 0.763932;

I = (h^3)/12;
beta = kappa_b/((2*(1+nu))*I);%*alpha;%0.3846*kappa_p/I
A = 1/(beta*(1-nu^2));
B = 1/(2*beta*(1+nu));

[MM,Kxx,Kxy,Kyy,Lx,Ly,Edge] = CalMatrix(n,m,deltx,delty);
LxT = Lx';
LyT = Ly';
Kyx = Kxy';%CHECKED
O = sparse(size(MM,1),size(MM,2));
Mu = [MM O O; O I*MM O; O O I*MM];
Ku = [Kxx+Kyy LxT LyT; h*Lx A*Kxx+B*Kyy+h*MM A*nu*Kyx+B*Kxy; h*Ly A*nu*Kxy+B*Kyx A*Kyy+B*Kxx+h*MM];%The correct one!
%Ku = [Kxx+Kyy Lx Ly; h*LxT A*Kxx+B*Kyy+h*MM A*nu*Kxy+B*Kyx; h*LyT A*nu*Kyx+B*Kxy A*Kyy+B*Kxx+h*MM];%The not correct one!
%Ku = [h*(Kxx+Kyy) h*Lx h*Ly; h*LxT A*Kxx+B*Kyy+h*MM A*nu*Kxy+B*Kyx; h*LyT A*nu*Kyx+B*Kxy A*Kyy+B*Kxx+h*MM];
Mq = [MM O O; O O O; O O O];
F = zeros(size(Mu,1),1);
F(1:(m+1),1) = 0.01;
x = [];

for i = [0 2 4 6 8 10]
    x = [x; Edge+(i)*(m+1)*(n+1)];
end
Mu(x,:) = [];
Mu(:,x) = [];
Ku(x,:) = [];
Ku(:,x) = [];
Mq(x,:) = [];

[R,p,s] = chol(Mu,'vector');
[V,DE,flag] = eigs(Ku,R,numEig,'smallestabs','IsCholesky',true,'CholeskyPermutation',s,'Tolerance',1e-4);
E = diag(DE);

%[V,D] = eigs(Ku,Mu,numEig,'sm');
%E = diag(D);
%V = V(:,E>=0);
%E = E(E>=0);

%tic
%Kug = gpuArray(Ku);
%bg = gpuArray(Mq*(F));
%u = gmres(Kug,bg,30,1e-4,30);
%ueq = gather(u);

%toc
%wP = [ueq(1:(m+1)*(n+1)-(m+1),1); zeros(m+1,1)];

%inum = 1
wP = zeros(inum,(m+1)*(n+1));
%wP(1,:) = [ueq(1:(m+1)*(n+1)-(m+1),1); zeros(m+1,1)];
for i = inum:-1:1
w = V(:,i);
tic
Kug = gpuArray(Ku);
bg = gpuArray(Mu*(w));
u = gmres(Kug,bg,30,1e-4,30);
ueq = gather(u);

%bg = Mu*(w);
%ueq = gmres(Ku,bg,30,1e-4,30);

toc
wP(i,:) = [ueq(1:(m+1)*(n+1)-(m+1),1); zeros(m+1,1)];
end

icx = b;
icy = d;
icount = 1;
xP = zeros(1,(n+1)*(m+1));
yP = zeros(1,(n+1)*(m+1));
for i = 1:n+1
   for j = 1:m+1
      xP(1,icount) = icx;
      yP(1,icount) = icy;
      icy = icy - delty;
      icount = icount + 1;
   end
   icx = icx - deltx;
   icy = d;
end

%for i = inum:-1:1
%figure();
%scatter3(xP(1,:),yP(1,:),wP(i,:));
%end


%}
%w = V(:,8);

%Kyx = Kxy';%CHECKED
%All = (n+1)*(m+1);

%K1 = Kxx + (1-nu)/2*Kyy; 
%K2 = nu*Kyx + (1-nu)/2*Kxy;
%K3 = nu*Kxy + (1-nu)/2*Kyx;
%K4 = Kyy + (1-nu)/2*Kxx;
%O = sparse(size(MM,1),size(MM,2));%CHECKED
%MMu = [MM O;%CHECKED
%       O MM];%CHECKED
%Mf = MMu;
%K = 1/(gamma*(1-nu^2))*[K1 K2; K3 K4];%CHECKED
%x = [7*All:-1:7*All-(m+1)+1 5*All:-1:5*All-(m+1)+1 3*All:-1:3*All-(m+1)+1 1*All:-1:1*All-(m+1)+1];
%K(x,:) = [];
%K(:,x) = [];
%MMu(x,:) = [];
%MMu(:,x) = [];
%Mf(x,:) = [];
%CHECKED
%eig(Mu,K)
%[V,D] = eigs(K,MMu,10,'sm');
%E = diag(D);
%w = V(:,8);




%ueq = Ku\Mu*(-w)

%b =Mf*(-F);
%tic
%ueq = K\Mf*(-F);
%toc

%{
tic
tol = 0.0001;
maxit = 300000;
alpha1 = max(sum(abs(K),2)./diag(K))-2;
L = ichol(K,struct('type','ict','droptol',1e-3,'diagcomp',alpha1));
%ueq = pcg(K,b,tol,maxit,L,L');
toc

Ep = Positions(m,n,deltx,delty);
%ux = 0;
%uy = 0;
ux = [ueq(1:(n+1)*(m+1)-(m+1),1);zeros(m+1,1)] + Ep(:,1);
dxux = [ueq((n+1)*(m+1)-(m+1)+1:2*(n+1)*(m+1)-(m+1),1)];
dyux = [ueq(2*(n+1)*(m+1)-(m+1)+1:3*(n+1)*(m+1)-2*(m+1),1);zeros(m+1,1)];
dxyux = [ueq(3*(n+1)*(m+1)-2*(m+1)+1:4*(n+1)*(m+1)-2*(m+1),1)];
uy = [ueq(4*(n+1)*(m+1)-2*(m+1)+1:5*(n+1)*(m+1)-3*(m+1),1);zeros(m+1,1)] + Ep(:,2);
dxuy = [ueq(5*(n+1)*(m+1)-3*(m+1)+1:6*(n+1)*(m+1)-3*(m+1),1)];
dyuy = [ueq(6*(n+1)*(m+1)-3*(m+1)+1:7*(n+1)*(m+1)-4*(m+1),1);zeros(m+1,1)];
dxyuy = [ueq(7*(n+1)*(m+1)-4*(m+1)+1:8*(n+1)*(m+1)-4*(m+1),1)];

ux = flip(ux);
dxux = flip(dxux);
dyux = flip(dyux);
dxyux = flip(dxyux);
uy = flip(uy);
dxuy = flip(dxuy);
dyuy = flip(dyuy);
dxyuy = flip(dxyuy);

stress = ceil((n+1)/2);
figure();
scatter(ux,uy);
sigma11 = 1/(gamma*(1-nu^2))*(dxux + nu*dyuy);
sigma22 = 1/(gamma*(1-nu^2))*(dyuy + nu*dxux);
sigma12 = 1/(2*gamma*(1+nu))*(dyux + dyux);

T = [sigma11(stress) sigma12(stress);sigma12(stress) sigma22(stress)];
%}
return;

function [Mq,Kxxq,Kxyq,Kyyq,Lxq,Lyq] = matrix(deltx,delty)%CHECKED
syms x;
syms y;

Q = [1 x x^2 x^3 y x*y x^2*y x^3*y y^2 x*y^2 x^2*y^2 x^3*y^2 y^3 x*y^3 x^2*y^3 x^3*y^3];
    
%size_num = size(Q,2)/4;
T = MATRIX_T(Q);

Mq   = zeros(size(Q,2))*x*y;
Kxxq = zeros(size(Q,2))*x*y;
Kxyq = zeros(size(Q,2))*x*y;
Kyyq = zeros(size(Q,2))*x*y;
Lxq  = zeros(size(Q,2))*x*y;
Lyq  = zeros(size(Q,2))*x*y;

for i = 1:size(Q,2)
    for j = 1:size(Q,2)
        Mq(j,i)   = Q(j)*Q(i);
        Kxxq(j,i) = diff(Q(j),x)*diff(Q(i),x);
        Kxyq(j,i) = diff(Q(j),y)*diff(Q(i),x);
        Kyyq(j,i) = diff(Q(j),y)*diff(Q(i),y);
        Lxq(j,i)  = Q(j)*diff(Q(i),x);
        Lyq(j,i)  = Q(j)*diff(Q(i),y);
    end
end
Mq = int(int(Mq,x,0,1),y,0,1);
Kxxq = int(int(Kxxq,x,0,1),y,0,1);
Kxyq = int(int(Kxyq,x,0,1),y,0,1);
Kyyq = int(int(Kyyq,x,0,1),y,0,1);
Lxq = int(int(Lxq,x,0,1),y,0,1);
Lyq = int(int(Lyq,x,0,1),y,0,1);

IT = inv(T);

Mq = (IT)'*Mq*IT;
Kxxq = (IT)'*Kxxq*IT;
Kxyq = (IT)'*Kxyq*IT;
Kyyq = (IT)'*Kyyq*IT;
Lxq = (IT)'*Lxq*IT;
Lyq = (IT)'*Lyq*IT;

Mq = double(Mq*deltx*delty);
Kxxq = double(Kxxq*delty/deltx);
Kyyq = double(Kyyq*deltx/delty);
Kxyq = double(Kxyq);
Lxq = double(Lxq*delty);
Lyq = double(Lyq*deltx);
return;

function [Adj,Type,Edge] = Domain(n,m)
D = zeros(m+1,n+1);
icount = 1;
for i = n+1:-1:1
   for j = 1:m+1
       D(j,i) = icount;
       icount = icount + 1;
   end
end
%D
Edge1 = [];
Edge2 = [];
Edge3 = [];
Edge4 = [];
Edge1 = (D(:,1));
%Edge2 = (D(m+1,:)');
%Edge3 = (D(:,n+1));
%Edge4 = (D(1,:)');

Edge = sort(unique([Edge1; Edge2; Edge3; Edge4]));
%Edge

D0 = [zeros(1,n+3);zeros(m+1,1) D zeros(m+1,1);zeros(1,n+3)];

icount = 1;
Adj = zeros((n+1)*(m+1),9);
for i = n+2:-1:2
   for j = 2:m+2
   Adj(icount,1) = D0(j,i); %middel
   Adj(icount,2) = D0(j-1,i); %bo
   Adj(icount,3) = D0(j-1,i+1); %regsbo
   Adj(icount,4) = D0(j,i+1); %regs
   Adj(icount,5) = D0(j+1,i+1); %regs onder
   Adj(icount,6) = D0(j+1,i); %onder
   Adj(icount,7) = D0(j+1,i-1); %links onder
   Adj(icount,8) = D0(j,i-1); %links
   Adj(icount,9) = D0(j-1,i-1); %linksbo
   icount = icount +1;
   end
end
Type= zeros((n+1)*(m+1),2);
T = [1     0     0     0     0     2     5     4     0;
     2     1     0     0     0     3     6     5     4;
     3     2     0     0     0     0     0     6     5;
     4     0     0     1     2     5     8     7     0;
     5     4     1     2     3     6     9     8     7;
     6     5     2     3     0     0     0     9     8;
     7     0     0     4     5     8     0     0     0;
     8     7     4     5     6     9     0     0     0;
     9     8     5     6     0     0     0     0     0];
 nnz(T);
 
for i = 1:(n+1)*(m+1)
   Type(i,1) = Adj(i,1);
   for j = 1:9
       bflag = true;
       for k = 1:9
           if(any(T(j,k)) ~= any(Adj(i,k)))
               bflag = false;
           end
       end
       if(bflag == true)
          Type(i,2) = j; 
       end
   end
end
return

function [M,Kxx,Kxy,Kyy,Lx,Ly,Edge] = CalMatrix(n,m,deltx,delty)
[Adj, Type, Edge] = Domain(n,m);
[Mq,Kxxq,Kxyq,Kyyq,Lxq,Lyq] = matrix(deltx,delty);

ns = nnz(Adj);
Ms = zeros(16*ns,1);
Kxxs = zeros(16*ns,1);
Kxys = zeros(16*ns,1);
Kyys = zeros(16*ns,1);
Lxs = zeros(16*ns,1);
Lys = zeros(16*ns,1);

x = [1:(n+1)*(m+1)]';
c = sum(Adj~=0,2);

ix = repelem(x,c);
x = [ix;ix+(n+1)*(m+1);ix+2*(n+1)*(m+1);ix+3*(n+1)*(m+1);
    ix;ix+(n+1)*(m+1);ix+2*(n+1)*(m+1);ix+3*(n+1)*(m+1);
    ix;ix+(n+1)*(m+1);ix+2*(n+1)*(m+1);ix+3*(n+1)*(m+1);
    ix;ix+(n+1)*(m+1);ix+2*(n+1)*(m+1);ix+3*(n+1)*(m+1);];

iy = nonzeros(Adj');
y = [iy;iy;iy;iy;
    iy+(n+1)*(m+1);iy+(n+1)*(m+1);iy+(n+1)*(m+1);iy+(n+1)*(m+1);
    iy+2*(n+1)*(m+1);iy+2*(n+1)*(m+1);iy+2*(n+1)*(m+1);iy+2*(n+1)*(m+1);
    iy+3*(n+1)*(m+1);iy+3*(n+1)*(m+1);iy+3*(n+1)*(m+1);iy+3*(n+1)*(m+1)];
B = BMatrix();

Adj(Adj == 0) = nan;
NanAdj = ~isnan(Adj);
NanAdj = NanAdj';

a = 1:9;
b = repelem(a,size(Adj,1),1);
b = b';
iz = b(NanAdj);

for i = 1:ns
  k = 1;
  while(B(Type(ix(i),2),iz(i),k) ~= 0) 
    Ms(i)   = Ms(i) + Mq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kxxs(i) = Kxxs(i) + Kxxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kxys(i) = Kxys(i) + Kxyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kyys(i) = Kyys(i) + Kyyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Lxs(i) = Lxs(i) + Lxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Lys(i) = Lys(i) + Lyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    
    Ms(ns+i)   = Ms(ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1));
    Kxxs(ns+i) = Kxxs(ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1));
    Kxys(ns+i) = Kxys(ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1));
    Kyys(ns+i) = Kyys(ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1));
    Lxs(ns+i) = Lxs(ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1));
    Lys(ns+i) = Lys(ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1));
    
    Ms(2*ns+i)   = Ms(2*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1));
    Kxxs(2*ns+i) = Kxxs(2*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1));
    Kxys(2*ns+i) = Kxys(2*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1));
    Kyys(2*ns+i) = Kyys(2*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1));
    Lxs(2*ns+i) = Lxs(2*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1));
    Lys(2*ns+i) = Lys(2*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1));
    
    Ms(3*ns+i)   = Ms(3*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1));
    Kxxs(3*ns+i) = Kxxs(3*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1));
    Kxys(3*ns+i) = Kxys(3*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1));
    Kyys(3*ns+i) = Kyys(3*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1));
    Lxs(3*ns+i) = Lxs(3*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1));
    Lys(3*ns+i) = Lys(3*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1));
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    Ms(4*ns+i)   = Ms(4*ns+i) + Mq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+4);
    Kxxs(4*ns+i) = Kxxs(4*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+4);
    Kxys(4*ns+i) = Kxys(4*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+4);
    Kyys(4*ns+i) = Kyys(4*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+4);
    Lxs(4*ns+i) = Lxs(4*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+4);
    Lys(4*ns+i) = Lys(4*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+4);
    
    Ms(5*ns+i)   = Ms(5*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+4);
    Kxxs(5*ns+i) = Kxxs(5*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+4);
    Kxys(5*ns+i) = Kxys(5*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+4);
    Kyys(5*ns+i) = Kyys(5*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+4);
    Lxs(5*ns+i) = Lxs(5*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+4);
    Lys(5*ns+i) = Lys(5*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+4);
    
    Ms(6*ns+i)   = Ms(6*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+4);
    Kxxs(6*ns+i) = Kxxs(6*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+4);
    Kxys(6*ns+i) = Kxys(6*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+4);
    Kyys(6*ns+i) = Kyys(6*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+4);
    Lxs(6*ns+i) = Lxs(6*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+4);
    Lys(6*ns+i) = Lys(6*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+4);
    
    Ms(7*ns+i)   = Ms(7*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+4);
    Kxxs(7*ns+i) = Kxxs(7*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+4);
    Kxys(7*ns+i) = Kxys(7*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+4);
    Kyys(7*ns+i) = Kyys(7*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+4);
    Lxs(7*ns+i) = Lxs(7*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+4);
    Lys(7*ns+i) = Lys(7*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+4);
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    Ms(8*ns+i)   = Ms(8*ns+i) + Mq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+8);
    Kxxs(8*ns+i) = Kxxs(8*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+8);
    Kxys(8*ns+i) = Kxys(8*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+8);
    Kyys(8*ns+i) = Kyys(8*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+8);
    Lxs(8*ns+i) = Lxs(8*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+8);
    Lys(8*ns+i) = Lys(8*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+8);
    
    Ms(9*ns+i)   = Ms(9*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+8);
    Kxxs(9*ns+i) = Kxxs(9*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+8);
    Kxys(9*ns+i) = Kxys(9*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+8);
    Kyys(9*ns+i) = Kyys(9*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+8);
    Lxs(9*ns+i) = Lxs(9*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+8);
    Lys(9*ns+i) = Lys(9*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+8);
    
    Ms(10*ns+i)   = Ms(10*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+8);
    Kxxs(10*ns+i) = Kxxs(10*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+8);
    Kxys(10*ns+i) = Kxys(10*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+8);
    Kyys(10*ns+i) = Kyys(10*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+8);
    Lxs(10*ns+i) = Lxs(10*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+8);
    Lys(10*ns+i) = Lys(10*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+8);
    
    Ms(11*ns+i)   = Ms(11*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+8);
    Kxxs(11*ns+i) = Kxxs(11*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+8);
    Kxys(11*ns+i) = Kxys(11*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+8);
    Kyys(11*ns+i) = Kyys(11*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+8);
    Lxs(11*ns+i) = Lxs(11*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+8);
    Lys(11*ns+i) = Lys(11*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+8);
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    Ms(12*ns+i)   = Ms(12*ns+i) + Mq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+12);
    Kxxs(12*ns+i) = Kxxs(12*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+12);
    Kxys(12*ns+i) = Kxys(12*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+12);
    Kyys(12*ns+i) = Kyys(12*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+12);
    Lxs(12*ns+i) = Lxs(12*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+12);
    Lys(12*ns+i) = Lys(12*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1)+12);
    
    Ms(13*ns+i)   = Ms(13*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+12);
    Kxxs(13*ns+i) = Kxxs(13*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+12);
    Kxys(13*ns+i) = Kxys(13*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+12);
    Kyys(13*ns+i) = Kyys(13*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+12);
    Lxs(13*ns+i) = Lxs(13*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+12);
    Lys(13*ns+i) = Lys(13*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+4,B(Type(ix(i),2),iz(i),k+1)+12);
    
    Ms(14*ns+i)   = Ms(14*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+12);
    Kxxs(14*ns+i) = Kxxs(14*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+12);
    Kxys(14*ns+i) = Kxys(14*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+12);
    Kyys(14*ns+i) = Kyys(14*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+12);
    Lxs(14*ns+i) = Lxs(14*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+12);
    Lys(14*ns+i) = Lys(14*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+8,B(Type(ix(i),2),iz(i),k+1)+12);
    
    Ms(15*ns+i)   = Ms(15*ns+i) + Mq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+12);
    Kxxs(15*ns+i) = Kxxs(15*ns+i) + Kxxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+12);
    Kxys(15*ns+i) = Kxys(15*ns+i) + Kxyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+12);
    Kyys(15*ns+i) = Kyys(15*ns+i) + Kyyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+12);
    Lxs(15*ns+i) = Lxs(15*ns+i) + Lxq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+12);
    Lys(15*ns+i) = Lys(15*ns+i) + Lyq(B(Type(ix(i),2),iz(i),k)+12,B(Type(ix(i),2),iz(i),k+1)+12);
    
    k = k + 2;
    if(k > 8)
       break; 
    end
  end
end
M   = sparse(x,y,Ms,4*(n+1)*(m+1),4*(n+1)*(m+1));
Kxx = sparse(x,y,Kxxs,4*(n+1)*(m+1),4*(n+1)*(m+1));
Kxy = sparse(x,y,Kxys,4*(n+1)*(m+1),4*(n+1)*(m+1));
Kyy = sparse(x,y,Kyys,4*(n+1)*(m+1),4*(n+1)*(m+1));
Lx = sparse(x,y,Lxs,4*(n+1)*(m+1),4*(n+1)*(m+1));
Ly = sparse(x,y,Lys,4*(n+1)*(m+1),4*(n+1)*(m+1));
return;

function B = BMatrix()
B = zeros(9,9,8);
B(3,1,:) = [2 2 0 0 0 0 0 0];
B(3,2,:) = [2 3 0 0 0 0 0 0];
B(3,3,:) = [0 0 0 0 0 0 0 0];
B(3,4,:) = [0 0 0 0 0 0 0 0];
B(3,5,:) = [0 0 0 0 0 0 0 0];
B(3,6,:) = [0 0 0 0 0 0 0 0];
B(3,7,:) = [0 0 0 0 0 0 0 0];
B(3,8,:) = [2 1 0 0 0 0 0 0];
B(3,9,:) = [2 4 0 0 0 0 0 0];

B(2,1,:) = [2 2 3 3 0 0 0 0];
B(2,2,:) = [2 3 0 0 0 0 0 0];
B(2,3,:) = [0 0 0 0 0 0 0 0];
B(2,4,:) = [0 0 0 0 0 0 0 0];
B(2,5,:) = [0 0 0 0 0 0 0 0];
B(2,6,:) = [3 2 0 0 0 0 0 0];
B(2,7,:) = [3 1 0 0 0 0 0 0];
B(2,8,:) = [2 1 3 4 0 0 0 0];
B(2,9,:) = [2 4 0 0 0 0 0 0];

B(1,1,:) = [3 3 0 0 0 0 0 0];
B(1,2,:) = [0 0 0 0 0 0 0 0];
B(1,3,:) = [0 0 0 0 0 0 0 0];
B(1,4,:) = [0 0 0 0 0 0 0 0];
B(1,5,:) = [0 0 0 0 0 0 0 0];
B(1,6,:) = [3 2 0 0 0 0 0 0];
B(1,7,:) = [3 1 0 0 0 0 0 0];
B(1,8,:) = [3 4 0 0 0 0 0 0];
B(1,9,:) = [0 0 0 0 0 0 0 0];

B(6,1,:) = [1 1 2 2 0 0 0 0];
B(6,2,:) = [1 4 2 3 0 0 0 0];
B(6,3,:) = [1 3 0 0 0 0 0 0];
B(6,4,:) = [1 2 0 0 0 0 0 0];
B(6,5,:) = [0 0 0 0 0 0 0 0];
B(6,6,:) = [0 0 0 0 0 0 0 0];
B(6,7,:) = [0 0 0 0 0 0 0 0];
B(6,8,:) = [2 1 0 0 0 0 0 0];
B(6,9,:) = [2 4 0 0 0 0 0 0];

B(5,1,:) = [1 1 2 2 3 3 4 4];
B(5,2,:) = [1 4 2 3 0 0 0 0];
B(5,3,:) = [1 3 0 0 0 0 0 0];
B(5,4,:) = [1 2 4 3 0 0 0 0];
B(5,5,:) = [4 2 0 0 0 0 0 0];
B(5,6,:) = [4 1 3 2 0 0 0 0];
B(5,7,:) = [3 1 0 0 0 0 0 0];
B(5,8,:) = [3 4 2 1 0 0 0 0];
B(5,9,:) = [2 4 0 0 0 0 0 0];

B(4,1,:) = [3 3 4 4 0 0 0 0];
B(4,2,:) = [0 0 0 0 0 0 0 0];
B(4,3,:) = [0 0 0 0 0 0 0 0];
B(4,4,:) = [4 3 0 0 0 0 0 0];
B(4,5,:) = [4 2 0 0 0 0 0 0];
B(4,6,:) = [4 1 3 2 0 0 0 0];
B(4,7,:) = [3 1 0 0 0 0 0 0];
B(4,8,:) = [3 4 0 0 0 0 0 0];
B(4,9,:) = [0 0 0 0 0 0 0 0];

B(9,1,:) = [1 1 0 0 0 0 0 0];
B(9,2,:) = [1 4 0 0 0 0 0 0];
B(9,3,:) = [1 3 0 0 0 0 0 0];
B(9,4,:) = [1 2 0 0 0 0 0 0];
B(9,5,:) = [0 0 0 0 0 0 0 0];
B(9,6,:) = [0 0 0 0 0 0 0 0];
B(9,7,:) = [0 0 0 0 0 0 0 0];
B(9,8,:) = [0 0 0 0 0 0 0 0];
B(9,9,:) = [0 0 0 0 0 0 0 0];

B(8,1,:) = [1 1 4 4 0 0 0 0];
B(8,2,:) = [1 4 0 0 0 0 0 0];
B(8,3,:) = [1 3 0 0 0 0 0 0];
B(8,4,:) = [1 2 4 3 0 0 0 0];
B(8,5,:) = [4 2 0 0 0 0 0 0];
B(8,6,:) = [4 1 0 0 0 0 0 0];
B(8,7,:) = [0 0 0 0 0 0 0 0];
B(8,8,:) = [0 0 0 0 0 0 0 0];
B(8,9,:) = [0 0 0 0 0 0 0 0];

B(7,1,:) = [4 4 0 0 0 0 0 0];
B(7,2,:) = [0 0 0 0 0 0 0 0];
B(7,3,:) = [0 0 0 0 0 0 0 0];
B(7,4,:) = [4 3 0 0 0 0 0 0];
B(7,5,:) = [4 2 0 0 0 0 0 0];
B(7,6,:) = [4 1 0 0 0 0 0 0];
B(7,7,:) = [0 0 0 0 0 0 0 0];
B(7,8,:) = [0 0 0 0 0 0 0 0];
B(7,9,:) = [0 0 0 0 0 0 0 0];
return;

function T = MATRIX_T(Q)
syms x;
syms y;

n = size(Q,2);

T = zeros(n);
for j = 1:n
   T(j,1) = subs(Q(j),[x,y],[0,0]);
   T(j,2) = subs(Q(j),[x,y],[1,0]);
   T(j,3) = subs(Q(j),[x,y],[1,1]);
   T(j,4) = subs(Q(j),[x,y],[0,1]);
   if(n > 4)
   T(j,5) = subs(diff(Q(j),x),[x,y],[0,0]);
   T(j,6) = subs(diff(Q(j),x),[x,y],[1,0]);
   T(j,7) = subs(diff(Q(j),x),[x,y],[1,1]);
   T(j,8) = subs(diff(Q(j),x),[x,y],[0,1]);
   end
   if(n > 8)
   T(j,9) = subs(diff(Q(j),y),[x,y],[0,0]);
   T(j,10) = subs(diff(Q(j),y),[x,y],[1,0]);
   T(j,11) = subs(diff(Q(j),y),[x,y],[1,1]);
   T(j,12) = subs(diff(Q(j),y),[x,y],[0,1]);
   end
   if(n > 12)
   T(j,13) = subs(diff(diff(Q(j),y),x),[x,y],[0,0]);
   T(j,14) = subs(diff(diff(Q(j),y),x),[x,y],[1,0]);
   T(j,15) = subs(diff(diff(Q(j),y),x),[x,y],[1,1]);
   T(j,16) = subs(diff(diff(Q(j),y),x),[x,y],[0,1]);
   end
end
T = T';
return

function E = Positions(m,n,dx,dy)
    E = zeros((n+1)*(m+1),2);
    ix = n+1;
    iy = m+1;
    for i = 1:(n+1)*(m+1)
        E(i,:) = [dx*(ix-1),dy*(iy-1)];
        
        iy = iy-1;
        if(iy == 0)
            iy = m+1;
            ix = ix -1;
        end
    end
    %[Cubes,CubeNumbers] = CreateCubes(E,N);
    %Plot(E,N,Cubes)
return