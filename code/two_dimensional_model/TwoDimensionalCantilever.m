function E = TwoDimensionalCantilever(n,alpha, showGraphs)
format long g
%gpuDevice(2)
beta = 1;
%alpha = 300;
gamma = 0.3205;
nu = 0.3;
iA = 1/(1-nu^2);
iB = 1/(2*gamma*(1+nu));

ex = sqrt(12/alpha);
m = ceil(n*ex);

if(m <= 1)
    m = 2;
end

n
m

a = 0;
b = 1;
c = 0;
d = sqrt(12/alpha);
deltx = (b-a)/n;
delty = (d-c)/m;

[MM,Kxx,Kxy,Kyy,D0] = CalMatrix(n,m,deltx,delty);
Kyx = Kxy';%CHECKED

Admin = (n+1)*(m+1) - (m+1);

K1 = Kxx(1:Admin,1:Admin) + (1-nu)/2*Kyy(1:Admin,1:Admin); 
K2 = nu*Kyx(1:Admin,1:Admin) + (1-nu)/2*Kxy(1:Admin,1:Admin);
K3 = nu*Kxy(1:Admin,1:Admin) + (1-nu)/2*Kyx(1:Admin,1:Admin);
K4 = Kyy(1:Admin,1:Admin) + (1-nu)/2*Kxx(1:Admin,1:Admin);
MMf = MM(1:Admin,:); 
MMu = MM(1:Admin,1:Admin);
Of = sparse(size(MMf,1),size(MMf,2));
Ou = sparse(size(MMu,1),size(MMu,2));

K = 1/(gamma*(1-nu^2))*[K1 K2; K3 K4];
Mf = [MMf Of; Of MMf];
Mu = [MMu Ou; Ou MMu];
%eig(Mu,K)
[V,D] = eigs(K,Mu,20,'sm');
E = diag(D);
%{
f = zeros(size(Mf,2));
f((n+1)*(m+1)+1:(n+1)*(m+1)+(m+1)) = 1;

%f((n+1)*(m+1)+1+ceil((n+1)*(m+1)):(n+1)*(m+1)+(m+1)+ceil((n+1)*(m+1))) = -1;

for i = 20:-1:1
    w = V(:,i);

    %f = -1/200;
    %F1 = zeros((n+1)*(m+1),1);
    %F1(ceil((1+(m+1))/2)) = f;
    %F = zeros(8*(n+1)*(m+1),1);
    %F(4*(n+1)*(m+1)+1:5*(n+1)*(m+1)) = F1;


    %ueq = K\MMu*(-w);
    tic
    %Kg = gpuArray(K);
    %Mfg = gpuArray(Mf);
    %Fg = gpuArray(F);
    toc

    %b =Mf*(-F);
    %tic
    %ueq = K\Mf*(-F);
    %toc
    b = Mu*(-w);
    %tic
    %tol = 0.00001;
    %maxit = 30000;
    %alpha1 = max(sum(abs(K),2)./diag(K))-2;
    %L = ichol(K,struct('type','ict','droptol',1e-4,'diagcomp',alpha1));
    %ueq = pcg(K,b,tol,maxit,L,L');
    ueq = K\b;
    %toc

    %b = Mf*(-f);
    %ueq = K\b;
    
    Ep = Positions(m,n,deltx,delty);
    Ep(:,1) = Ep(:,1)*5;
    %ux = 0;
    %uy = 0;
    size(ueq)
    ux = [ueq(1:(n+1)*(m+1)-(m+1),1);zeros(m+1,1)] + Ep(:,1);
    uy = [ueq((n+1)*(m+1)-(m+1)+1:2*(n+1)*(m+1)-2*(m+1),1);zeros(m+1,1)] + Ep(:,2);
    
    ux = flip(ux);
    uy = flip(uy);
    
    top = D0(floor(m/2)+1,:);
    bot = D0(floor(m/2)-1,:);
    
    %top = D0(m+1,:);
    %bot = D0(1,:);
    
    
    uxmid = ux(D0(ceil((m+1)/2),:));
    uymid = uy(D0(ceil((m+1)/2),:));
    
    
    
    maxsx = norm(ux,Inf);
    maxs = norm(uy,Inf);
    uy = uy/maxs;
    ux = ux/maxsx;
    h = figure(100+i);
    movegui(h,'east');
    
    
    scatter(ux,uy,'linewidth',2);
    hold on;
    %plot(uxmid,uymid);
    
    %for j = 1:size(D0,2)
    %   plot(ux(D0(:,j)),uy(D0(:,j)),'r-','linewidth',3); 
    %   plot([ux(D0(1,j)) ux(D0(m+1,j))],[uy(D0(1,j)) uy(D0(m+1,j))],'k--','linewidth',1)
    %end
    %legend('Line 1','Line 2','Line 3','Line 4')
    ux1 = ux(D0(1,:));
    uy1 = uy(D0(1,:));
    
    maxs2 =  norm(uy1,Inf);
    %scatter(ux1,uy1/maxs2);
    %plot(ux1,uy1/maxs2);
    xx = 0:deltx:1-deltx;
    top(1) = [];
    bot(1) = [];
    TOPx = ux(top);
    TOPy = uy(top);
    BOTx = ux(bot);
    BOTy = uy(bot);
    phi = zeros(1,size(top,2));
    for j = 1:size(top,2)
        u = [TOPx(j) TOPy(j)];
        v = [BOTx(j) BOTy(j)];
        
        vert = u-v;
        hor = [1 0];
        phi(j)= angle(vert,hor);
    end
   
    phi = -flip(phi);
    phi = phi - phi(1);
    phi = phi/norm(phi,Inf)*0.95;
    %plot(xx,phi,'o');
    %scatter(TOPx,TOPy)
    %scatter(BOTx,BOTy)
end

%}
return;

function a = angle(u,v)
a = 0;
dotuv = dot(u,v);
cosalpha = dotuv/(norm(u)*norm(v));
a = acos(cosalpha);
return;

function [Mq,Kxxq,Kxyq,Kyyq] = matrix(deltx,delty)%CHECKED
syms x;
syms y;

Q = [1 x y x*y];%CHECKED

Mq = zeros(size(Q,2))*x*y;%CHECKED
Kxxq = zeros(size(Q,2))*x*y;%CHECKED
Kxyq = zeros(size(Q,2))*x*y;%CHECKED
Kyyq = zeros(size(Q,2))*x*y;%CHECKED

for i = 1:size(Q,2)%CHECKED
    for j = 1:size(Q,2)%CHECKED
        Mq(j,i) = Q(j)*Q(i);%CHECKED
        Kxxq(j,i) = diff(Q(j),x)*diff(Q(i),x);%CHECKED
        Kxyq(j,i) = diff(Q(j),y)*diff(Q(i),x);%CHECKED
        Kyyq(j,i) = diff(Q(j),y)*diff(Q(i),y);%CHECKED
    end
end

Mq = int(int(Mq,x,0,1),y,0,1);%CHECKED
Kxxq = int(int(Kxxq,x,0,1),y,0,1);%CHECKED
Kxyq = int(int(Kxyq,x,0,1),y,0,1);%CHECKED
Kyyq = int(int(Kyyq,x,0,1),y,0,1);%CHECKED


%T = [1 0 0 0; 1 1 0 0; 1 1 1 1; 1 0 1 0];
IT = [1 0 0 0; -1 1 0 0; -1 0 0 1; 1 -1 1 -1];%CHECKED

Mq = double((IT)'*Mq*IT);%CHECKED
Kxxq = double((IT)'*Kxxq*IT);%CHECKED
Kxyq = double((IT)'*Kxyq*IT);%CHECKED
Kyyq = double((IT)'*Kyyq*IT);%CHECKED

Mq = Mq*deltx*delty;%CHECKED
Kxxq = Kxxq*delty/deltx;%CHECKED
Kyyq = Kyyq*deltx/delty;%CHECKED
return;

function [Adj, Type, D] = Domain(n,m)
D = zeros(m+1,n+1);
icount = 1;
for i = n+1:-1:1
   for j = 1:m+1
       D(j,i) = icount;
       icount = icount + 1;
   end
end
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

function [M,Kxx,Kxy,Kyy,D0] = CalMatrix(n,m,deltx,delty)
[Adj, Type, D0] = Domain(n,m);
[Mq,Kxxq,Kxyq,Kyyq] = matrix(deltx,delty);

ns = nnz(Adj);
Ms = zeros(ns,1);
Kxxs = zeros(ns,1);
Kxys = zeros(ns,1);
Kyys = zeros(ns,1);

x = [1:(n+1)*(m+1)]';
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

for i = 1:ns
  k = 1;
  while(B(Type(ix(i),2),iz(i),k) ~= 0) 
    Ms(i)   = Ms(i) + Mq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kxxs(i) = Kxxs(i) + Kxxq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kxys(i) = Kxys(i) + Kxyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
    Kyys(i) = Kyys(i) + Kyyq(B(Type(ix(i),2),iz(i),k),B(Type(ix(i),2),iz(i),k+1));
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