function [E,wP,xP,yP,size_l] = PlateCantilever(n,h,inum,numEig,d)
format long g
m = ceil(n*d);

a = 0;
b = 1;
c = 0;
%alpha = 300;
deltx = (b-a)/n;
delty = (d-c)/m;

size_l = (n+1)*(m+1);

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
%Ku = [h*(Kxx+Kyy) h*Lx h*Ly; h*LxT A*Kxx+B*Kyy+h*MM A*nu*Kxy+B*Kyx; h*LyT A*nu*Kyx+B*Kxy A*Kyy+B*Kxx+h*MM];
Mq = [MM O O; O O O; O O O];

F = zeros(size(Mu,1),1);
F(1:(m+1),1) = -100;
x = [];
for i = [0 1 2]
    x = [x; Edge+(i)*(m+1)*(n+1)];
end

Mu(x,:) = [];
Mu(:,x) = [];
Ku(x,:) = [];
Ku(:,x) = [];
Mq(x,:) = [];

[V,D] = eigs(Ku,Mu,numEig,'sm');
E = diag(D);
%V = V(:,E>=0);
%E = E(E>=0);

wP = zeros(inum,(m+1)*(n+1));
%wP(1,:) = [ueq(1:(m+1)*(n+1)-(m+1),1); zeros(m+1,1)];
for i = inum:-1:1
w = V(:,i);
tic
%Kug = gpuArray(Ku);
%bg = gpuArray(Mu*(w));
%u = gmres(Kug,bg,30,1e-4,30);
%ueq = gather(u);

bg = Mu*(w);
ueq = gmres(Ku,bg,30,1e-4,30);

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
%}

%for i = inum:-1:1
%figure();
%scatter3(xP(1,:),yP(1,:),wP(i,:),5,wP(i,:));
%end
return;

function [Mq,Kxxq,Kxyq,Kyyq,Lyq,Lxq] = matrix(deltx,delty)%CHECKED
syms x;
syms y;

Q = [1 x y x*y];%CHECKED

Mq = zeros(size(Q,2))*x*y;%CHECKED
Kxxq = zeros(size(Q,2))*x*y;%CHECKED
Kxyq = zeros(size(Q,2))*x*y;%CHECKED
Kyyq = zeros(size(Q,2))*x*y;%CHECKED
Lxq = zeros(size(Q,2))*x*y;%CHECKED
Lyq = zeros(size(Q,2))*x*y;%CHECKED

for i = 1:size(Q,2)%CHECKED
    for j = 1:size(Q,2)%CHECKED
        Mq(j,i) = Q(j)*Q(i);%CHECKED
        Kxxq(j,i) = diff(Q(j),x)*diff(Q(i),x);%CHECKED
        Kxyq(j,i) = diff(Q(j),y)*diff(Q(i),x);%CHECKED
        Kyyq(j,i) = diff(Q(j),y)*diff(Q(i),y);%CHECKED
        Lxq(j,i) = Q(j)*diff(Q(i),x);%CHECKED
        Lyq(j,i) = Q(j)*diff(Q(i),y);%CHECKED
    end
end

Mq = int(int(Mq,x,0,1),y,0,1);%CHECKED
Kxxq = int(int(Kxxq,x,0,1),y,0,1);%CHECKED
Kxyq = int(int(Kxyq,x,0,1),y,0,1);%CHECKED
Kyyq = int(int(Kyyq,x,0,1),y,0,1);%CHECKED
Lxq = int(int(Lxq,x,0,1),y,0,1);%CHECKED
Lyq = int(int(Lyq,x,0,1),y,0,1);%CHECKED

%T = inv(MATRIX_T(Q));
%IT = inv([1 0 0 0; 1 1 0 0; 1 1 1 1; 1 0 1 0])
IT = [1 0 0 0; -1 1 0 0; -1 0 0 1; 1 -1 1 -1];%CHECKED

Mq = double((IT)'*Mq*IT);%CHECKED
Kxxq = double((IT)'*Kxxq*IT);%CHECKED
Kxyq = double((IT)'*Kxyq*IT);%CHECKED
Kyyq = double((IT)'*Kyyq*IT);%CHECKED
Lxq = double((IT)'*Lxq*IT);%CHECKED
Lyq = double((IT)'*Lyq*IT);%CHECKED

Mq = double(Mq*deltx*delty);%CHECKED
Kxxq = double(Kxxq*delty/deltx);%CHECKED
Kyyq = double(Kyyq*deltx/delty);%CHECKED
Lxq = double(Lxq*delty);%CHECKED
Lyq = double(Lyq*deltx);%CHECKED
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

function [Adj, Type, Edge, D] = Domain(n,m)
D = zeros(m+1,n+1);
icount = 1;
for i = n+1:-1:1
   for j = 1:m+1
       D(j,i) = icount;
       icount = icount + 1;
   end
end
Edge1 = [];
Edge2 = [];
Edge3 = [];
Edge4 = [];

Edge1 = (D(:,1)); 
%Edge2 = (D(m+1,:)');
%Edge3 = (D(:,n+1));
%Edge4 = (D(1,:)');

Edge = sort(unique([Edge1; Edge2; Edge3; Edge4]));

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
[Adj, Type ,Edge, D] = Domain(n,m);
[Mq,Kxxq,Kxyq,Kyyq,Lyq,Lxq] = matrix(deltx,delty);

ns = nnz(Adj);
Ms = zeros(ns,1);
Kxxs = zeros(ns,1);
Kxys = zeros(ns,1);
Kyys = zeros(ns,1);
Lxs = zeros(ns,1);
Lys = zeros(ns,1);

x = reshape(flip(D')',[],1);
c = sum(Adj~=0,2);

ix = repelem(x,c);
iy = nonzeros(Adj');
B = BMatrix();

Adj(Adj == 0) = nan;
NanAdj = ~isnan(Adj);
NanAdj = NanAdj';

a = 1:9;
b = repelem(a,size(Adj,1),1);
b = b';
iz = b(NanAdj);

[ix iy iz];

for i = 1:ns
  k = 1;
  while(B(Type(ix(i),2),iz(i),k) ~= 0) 
    Ms(i)   = Ms(i) + Mq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kxxs(i) = Kxxs(i) + Kxxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kxys(i) = Kxys(i) + Kxyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kyys(i) = Kyys(i) + Kyyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Lxs(i) = Lxs(i) + Lxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Lys(i) = Lys(i) + Lyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    k = k + 2;
    if(k > 8)
       break; 
    end
  end
end
M   = sparse(ix,iy,Ms,(n+1)*(m+1),(n+1)*(m+1));
Kxx = sparse(ix,iy,Kxxs,(n+1)*(m+1),(n+1)*(m+1));
Kxy = sparse(ix,iy,Kxys,(n+1)*(m+1),(n+1)*(m+1));
Kyy = sparse(ix,iy,Kyys,(n+1)*(m+1),(n+1)*(m+1));
Lx = sparse(ix,iy,Lxs,(n+1)*(m+1),(n+1)*(m+1));
Ly = sparse(ix,iy,Lys,(n+1)*(m+1),(n+1)*(m+1));
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