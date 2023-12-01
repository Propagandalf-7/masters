 function [Eig,uxB,uyB,uzB,sz] = CubePC(s,n1,h,inum,numEig)
format long g
warning off;
%mkdir(strcat('\Plots\',sprintf('%.6f',s)));
method = 2;
%gpuDevice(1);
    mwb = MultiWaitBar(3, 1, '3-Dimensional Beam Eigenvalue Calculator', 'g');
    mwb.Update(1, 1, 0, 'Total Progress - Setting parameters');
    mwb.Update(2, 1, 0, ['Matrix Creation ' num2str(0) '%']);
    %mwb.Update(3, 1, 0, 'Plot');
    %alpha = 1200;
    %d2 = sqrt(s/alpha);
    %d1 = sqrt((12*s^2)/(alpha*(1+s^2)));
    d1 = h; 
    d2 = s; 
    %d2 = 1; 
    
    n2 = ceil(n1*h);
    
    if(n2 <= 1)
        n2 = 2;
    end
    if(n2 >= 15)
        n2 = 15;
    end
    n3 = ceil(n2*d2);
    if(n3 <= 1)
        n3 = 2;
    end
    if(n3 >= 15)
        n3 = 15;
    end
    %[n1 n2 n3]
    %n3 = 30;
    %n2 = 2;
    
    sz = n1*n2*n3;
    
    S = [0 1 0 d1 0 d2]; %Set size of the beam
    N = [n1 n2 n3]; %Number of elements
    Delta = [(S(2)-S(1))/N(1) (S(4)-S(3))/N(2) (S(6)-S(5))/N(3)]; %space step size
    
    gamma = 0.3205;
    nu = 0.3;
    f = 0.001;
    %A = 1/(gamma*(1+nu)*(1-2*nu));
    %B = 1/(2*gamma*(1+nu));
    
    %A = 1/(1-nu^2);
    %B = 1/(2*gamma*(1+nu));
    mwb.Update(1, 1, 0.1, 'Total Progress - Creating Matrices');
    [K11,K12,K13,K22,K23,K33,M0,Dom,E] = Matrices(Delta,N,mwb,method);
    mwb.Update(1, 1, 0.3, 'Total Progress - Admissible Basis functions');
    %Om = Omega(N,Dom);
    clear functions
    F = Initial(N,f);
    Mf = M0;
   
   K11(1:(N(2)+1)*(N(3)+1),:) = [];
   K11(:,1:(N(2)+1)*(N(3)+1)) = [];
   K12(1:(N(2)+1)*(N(3)+1),:) = [];
   K12(:,1:(N(2)+1)*(N(3)+1)) = [];
   K13(1:(N(2)+1)*(N(3)+1),:) = [];
   K13(:,1:(N(2)+1)*(N(3)+1)) = [];
   K22(1:(N(2)+1)*(N(3)+1),:) = [];
   K22(:,1:(N(2)+1)*(N(3)+1)) = [];
   K23(1:(N(2)+1)*(N(3)+1),:) = [];
   K23(:,1:(N(2)+1)*(N(3)+1)) = [];
   K33(1:(N(2)+1)*(N(3)+1),:) = [];
   K33(:,1:(N(2)+1)*(N(3)+1)) = [];

   M0(1:(N(2)+1)*(N(3)+1),:) = [];
   M0(:,1:(N(2)+1)*(N(3)+1)) = [];
   Mf(1:(N(2)+1)*(N(3)+1),:) = []; 
    mwb.Update(1, 1, 0.4, 'Total Progress - Concatinating matrices');
    Of = sparse(size(Mf,1),size(Mf,2));
    MF = [Mf Of Of; Of Mf Of; Of Of Mf];
    O = sparse(size(M0,1),size(M0,2));
    M = [M0 O O; O M0 O; O O M0];
    
    %M = [M0 O; O M0];
    K21 = K12';
    K31 = K13';
    K32 = K23';
    
    
    a1 = 1/(gamma*(1+nu));
    a2 = nu/(gamma*(1+nu)*(1-2*nu));
    a3 = 1/(2*gamma*(1+nu));
    
    K1 = a1*K11 + a2*K11 + a3*K22 + a3*K33;
    K2 = a3*K12 + a2*K21;
    K3 = a3*K13 + a2*K31;
    K4 = a2*K12 + a3*K21;
    K5 = a1*K22 + a2*K22 + a3*K11 + a3*K33;
    K6 = a3*K23 + a2*K32;
    K7 = a2*K13 + a3*K31;
    K8 = a2*K23 + a3*K32;
    K9 = a1*K33 + a2*K33 + a3*K11 + a3*K22;
    
    K = [K1 K2 K3; K4 K5 K6; K7 K8 K9];
   % K = [1/(gamma*(1+nu))*K11+nu/(gamma*(1+nu)*(1-2*nu))*(K11+K22+K33)  1/(gamma*(1+nu))*K12 1/(gamma*(1+nu))*K13;
  %       1/(gamma*(1+nu))*K12 1/(gamma*(1+nu))*K22+nu/(gamma*(1+nu)*(1-2*nu))*(K11+K22+K33) 1/(gamma*(1+nu))*K23;
   %      1/(gamma*(1+nu))*K13 1/(gamma*(1+nu))*K23 1/(gamma*(1+nu))*K33+nu/(gamma*(1+nu)*(1-2*nu))*(K11+K22+K33)];
    
   % K33 = (-nu/((1-2*nu)+nu))*(K11+K22);
   % K = [1/(gamma*(1+nu))*K11+nu/(gamma*(1+nu)*(1-2*nu))*(K11+K22+K33)  1/(gamma*(1+nu))*K12;
   %      1/(gamma*(1+nu))*K12 1/(gamma*(1+nu))*K22+nu/(gamma*(1+nu)*(1-2*nu))*(K11+K22+K33)]; 
     
     
    %K = [2*(1-nu)*K11+(1-2*nu)*K22+(1-2*nu)*K33 2*nu*K12+(1-2*nu)*K21 2*nu*K13+(1-2*nu)*K31; 
    %    2*nu*K21+(1-2*nu)*K12 (1-2*nu)*K11+2*(1-nu)*K22+(1-2*nu)*K33 2*nu*K23+(1-2*nu)*K32; 
    %    2*nu*K31+(1-2*nu)*K13 2*nu*K32+(1-2*nu)*K23 (1-2*nu)*K11+(1-2*nu)*K22+2*(1-nu)*K33];
    %K = 1/(2*gamma*(1+nu)*(1-2*nu))*K;
    
    %K = [K11+(1-nu)/2*K22+(1-nu)/2*K33 (1-nu)/2*K12+nu*K21 (1-nu)/2*K13+nu*K31; 
    %    (1-nu)/2*K21+nu*K12 (1-nu)/2*K11+K22+(1-nu)/2*K33 (1-nu)/2*K23+nu*K32; 
    %    (1-nu)/2*K31+nu*K13 (1-nu)/2*K32+nu*K23 (1-nu)/2*K11+(1-nu)/2*K22+K33];
    % K = 1/(gamma*(1-nu^2))*K;
     whos k
   

%{
alpha = max(sum(abs(K),2)./diag(K))-2;
L = ichol(K,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
n = size(K,1);
[V,D] = eigs(@(x)pcg(K,x,1e-3,200,L,L'),n,M,numEig,'sm');
%}
 mwb.Update(1, 1, 0.5, 'Total Progress - Cholsky Decomposition');
[R,p,s] = chol(M,'vector');
 mwb.Update(1, 1, 0.55, 'Total Progress - Eigs');

%Rand = sprand(K);
%[v, lambda] = lobpcg(Rand, K, M, 1e-5, 20,0) 
 
[V,D,flag] = eigs(K,R,numEig,'smallestabs','IsCholesky',true,'CholeskyPermutation',s,'Tolerance',1e-4);
mwb.Update(1, 1, 0.6, 'Total Progress');
Eig = diag(D);
%%Mg = gpuArray(M);
%%Kg = gpuArray(K);

%{
sV = size(Eig,1);
R = zeros(sV,sV);
for i = 1:sV
   for j = 1:sV
    X = K*V(:,i) - M*V(:,i)*D(j);
    NORMX = norm(X,Inf);
    R(j,i) = NORMX;
   end
end
R = K*V-M*V*D;
xlswrite('CompareEigenValues.xlsx',R)
%}

u1p = 0;
u2p = 0;
u3p = 0;
u1s = 0;
u2s = 0;
u3s = 0;
uplx= 0;
uply = 0;
Psize = 0;
%%{
%[KM, KMPat] = sparseinv(K);
%KM = KM*M;
%smallest = 0;
%bestEig = 0;
uxB = zeros(inum,N(1)+1,N(2)+1,1);
uyB = zeros(inum,N(1)+1,N(2)+1,1);
uzB = zeros(inum,N(1)+1,N(2)+1,1);
for i = inum:-1:1
  %smallest = 100;  
  %plot_fig = figure('NumberTitle', 'off', 'Name', strcat('Eigenvalue: ',int2str(i),' - ',num2str(s)));
  w = V(:,i);
  
%  for j = 1:numEig
%    X = K*w - M*w*Eig(j);
%    if norm(X) < smallest
%        smallest = norm(X);
%        bestEig = Eig(j);
%    end
%  end
   %plot_fig.suptitle(strcat(int2str(i),' - ',num2str(s)));
%mwb.Update(3, 1, 0.1, 'Plot');
    %wg = gpuArray(w);
    alpha = max(sum(abs(K),2)./diag(K))-2;
    L = ichol(K,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    u = pcg(K,M*w,1e-1,2000000,L,L');
   
   %u = gmres(K,M*(-w));
    %normalize = norm(w,Inf);
    %u = normalize*u;
   % u = (K)\MF*F;
   
    %alpha = max(sum(abs(K),2)./diag(K))-2;
    %L = ichol(K,struct('type','ict','droptol',1e-3,'diagcomp',alpha));
    %u = pcg(K,MF*F,1e-3,200000,L,L');
   % mwb.Update(3, 1, 0.3, 'Plot');
    u1 = [zeros((N(2)+1)*(N(3)+1),1); u(1:size(u,1)/3)] + E(:,1)*50;
    u2 = [zeros((N(2)+1)*(N(3)+1),1); u(size(u,1)/3+1:2*size(u,1)/3)] + 20*E(:,2);
    u3 = [zeros((N(2)+1)*(N(3)+1),1); u(2*size(u,1)/3+1:3*size(u,1)/3)]+ E(:,3); 
    
    
    u1 = 1/norm(u1,Inf)*u1;
    %u2 = 1/norm(u2,Inf)*u2;
    u3 = 1/norm(u3,Inf)*u3;
    
    [D,E] = Domain(N,Delta);
    %plane = D(:,ceil((N(2)+1)/2),:);
    plane = D(:,:,:);
    plane = D(:,:,ceil((N(3)+1)/2));
    plane = plane(:);
    
   % w1 = w(1:size(w,1)/3);
   % w2 = w(size(w,1)/3+1:2*size(w,1)/3);
   % w3 = w(2*size(w,1)/3+1:3*size(w,1)/3);

  % mwb.Update(3, 1, 0.4, 'Plot');
    %{
    ix = [];
    if(N(1)+1 > 200)
       % iy = [1 (N(2)+2)/2 (N(2)+1) (N(2)+1)*(N(3)+1)-(N(2)+1)+1 (N(2)+1)*(N(3)+1)];
        ih = ((N(2)+1)-1)/2;
        iv = 1+((N(3)+1)-1)/2*(N(2)+1);
        iy = iv+ih; %[1 1+ih 1+2*ih iv iv+ih iv+2*ih 2*iv-1 2*iv+ih-1 2*iv+2*ih-1];%[iv+ih];%
        icount = 1;
        div = floor((N(1)+1)/200);
        for k = 1:div:(N(1)+1)
            for j = 1:size(iy,2)
                ix(icount) = iy(j)+ (k-1)*(N(2)+1)*(N(3)+1);
                icount = icount +1;
            end
        end
        u1p = u1(ix);
        u2p = u2(ix);
        u3p = u3(ix);
       % w1p = w1(ix);
       % w2p = w2(ix);
       % w3p = w3(ix);

        E1p = E(ix,1);
        E2p = E(ix,2);
        E3p = E(ix,3);
        Psize = size(iy,2);
    else
        u1p = u1;
        u2p = u2;
        u3p = u3;
      
      %  w1p = w1;
      %  w2p = w2;
      %  w3p = w3;

        E1p = E(:,1);
        E2p = E(:,2);
        E3p = E(:,3);
        
        Psize = (N(2)+1)*(N(3)+1);
    end
  %}
    %mwb.Update(3, 1, 0.7, 'Plot');
    
    %Psize = (N(3)+1);
    
    u1p = 1/norm(u1,Inf)*u1;
    u2p = 1/norm(u2,Inf)*u2;
    u3p = 1/norm(u3,Inf)*u3;
    
    u1p = u1;
    %u2p = u2;
    u3p = u3;
    
    u1p = u1p(plane,1);
    u2p = u2p(plane,1);
    u3p = u3p(plane,1);
    
    %uxB(i,:) = u1p;
    %uyB(i,:) = u2p;
    %uzB(i,:) = u3p;
    
    %umx = u1p(size(u1p,1)/2
    f = figure(i);
    movegui(f,'west')
    scatter3(u1p,-u2p,u3p);
    title(Eig(i));
    %hold on
    %scatter3(E1p,E2p,E3p,0.1);
    %hold on
    %scatter3(w1p,w2p,w3p);
    %mwb.Update(3, 1, 1, 'Plot');
    %u1s = size(u1p);
    %u2s = size(u2p);
    %u3s = size(u3p);
    %for k = 1:Psize
    %    hold on
    %    plot3(u1p(k:Psize:k+u1s-2*Psize),u2p(k:Psize:k+u2s-2*Psize),u3p(k:Psize:k+u3s-2*Psize),'-');
    %end
%temp_png = strcat('\Plots\',sprintf('%.6f',n1),'\PNG\Plot',sprintf('%.6f',i),'.png');
%temp_fig = strcat('\Plots\',sprintf('%.6f',n1),'\Fig\Plot',sprintf('%.6f',i),'.fig');
%view([0 0 90])
%legend(['Eigenvalue: ' num2str(Eig(i))]);
%saveas(plot_fig,strcat(pwd,temp_png))
%savefig(plot_fig,strcat(pwd,temp_fig))
%close(plot_fig)
end
clear u

%}
mwb.Update(1, 1, 1, 'Total Progress');
   mwb.Close();
return;

function F = Initial(N,f)
F = zeros(3*(N(1)+1)*(N(2)+1)*(N(3)+1),1);
for i = 3*((N(1)+1)*(N(2)+1)*(N(3)+1))-((N(1)+1)*(N(2)+1))+1:3*((N(1)+1)*(N(2)+1)*(N(3)+1))
    F(i) = f;
end
return

function plotCube(C,i,N)
    hold on
    colour = (i/(N(1)*N(2)*N(3))*255)/255*[1 0 0];
    plot3([C(1) C(2)],[C(3) C(3)],[C(5) C(5)],'color',colour);
    plot3([C(1) C(2)],[C(4) C(4)],[C(5) C(5)],'color',colour);
    plot3([C(1) C(2)],[C(3) C(3)],[C(6) C(6)],'color',colour);
    plot3([C(1) C(2)],[C(4) C(4)],[C(6) C(6)],'color',colour);

    plot3([C(1) C(1)],[C(3) C(4)],[C(5) C(5)],'color',colour);
    plot3([C(2) C(2)],[C(3) C(4)],[C(5) C(5)],'color',colour);
    plot3([C(1) C(1)],[C(3) C(4)],[C(6) C(6)],'color',colour);
    plot3([C(2) C(2)],[C(3) C(4)],[C(6) C(6)],'color',colour);

    plot3([C(1) C(1)],[C(3) C(3)],[C(5) C(6)],'color',colour);
    plot3([C(1) C(1)],[C(4) C(4)],[C(5) C(6)],'color',colour);
    plot3([C(2) C(2)],[C(3) C(3)],[C(5) C(6)],'color',colour);
    plot3([C(2) C(2)],[C(4) C(4)],[C(5) C(6)],'color',colour);
return

function [D,E] = Domain(N,Delta)
    D = zeros(N(1)+1,N(2)+1,N(3)+1);
    icount = 1;
    
    for i = 1:N(1)+1
       for k = 1:N(3)+1
          for j = 1:N(2)+1
             D(i,j,k) = icount;
             icount = icount + 1;
          end
       end
    end
    E = zeros((N(1)+1)*(N(2)+1)*(N(3)+1),3);
    ix = 1;
    iy = 1;
    iz = 1;
    ixt = 0;
    for i = 1:(N(1)+1)*(N(2)+1)*(N(3)+1)
        E(i,:) = [Delta(1)*(ix-1),Delta(2)*(iy-1),Delta(3)*(iz-1)];
        
        iy = iy+1;
        
        if(ix == N(1)+2)
            ix = 1;
        end
        if(iy == N(2)+2)
            iy = 1;
            ixt = ixt +1;
            iz = iz+1;
        end
        if(ixt == N(3)+1)
            ix = ix+1;
            ixt = 0;
        end
        if(iz == N(3)+2)
            iz = 1;
        end
        
    end
    %[Cubes,CubeNumbers] = CreateCubes(E,N);
    %Plot(E,N,Cubes)
return

function [Cubes,CubeNumbers] = CreateCubes(E,N)
Cubes = zeros((N(1))*(N(2))*(N(3)),6);
CubeNumbers = zeros(N(1)*N(2)*N(3),8);
it = 0;
ic = 1;
ik = 1;
for i = 1:N(1)*N(2)*N(3)
    if(it == N(2))
        it = 0;
        ic = ic + (N(2)+1);
    end
    if ic <= ik*(N(2)+1)*(N(3)+1) && ic >=  ik*(N(2)+1)*(N(3)+1)-(N(2)+1)
       ic = ic + (N(2)+1);
       ik = ik +1;
    end
    it = it + 1;
    CubeNumbers(i,1) = (it-1) + ic;
    CubeNumbers(i,2) = CubeNumbers(i,1)+1;
    CubeNumbers(i,3) = CubeNumbers(i,1)+(N(2)+1);
    CubeNumbers(i,4) = CubeNumbers(i,3)+1;
    CubeNumbers(i,5) = CubeNumbers(i,1)+(N(2)+1)*(N(3)+1);
    CubeNumbers(i,6) = CubeNumbers(i,5)+1;
    CubeNumbers(i,7) = CubeNumbers(i,5)+(N(2)+1);
    CubeNumbers(i,8) = CubeNumbers(i,7)+1;
end
for i = 1:N(1)*N(2)*N(3)
   Cubes(i,1) = E(CubeNumbers(i,1),1);
   Cubes(i,2) = E(CubeNumbers(i,5),1);
   Cubes(i,3) = E(CubeNumbers(i,1),2);
   Cubes(i,4) = E(CubeNumbers(i,6),2);
   Cubes(i,5) = E(CubeNumbers(i,1),3);
   Cubes(i,6) = E(CubeNumbers(i,3),3);
end
return

function Plot(E,N,Cubes)
    scatter3(E(:,1),E(:,2),E(:,3))
    A = [1:(N(1)+1)*(N(2)+1)*(N(3)+1)]';
    Astr = num2str(A);
    Acell = cellstr(Astr);
    dx = 0;
    dy = 0;
    text(E(:,1)+dx,E(:,2)+dy,E(:,3),Acell)

   for i = 1:N(1)*N(2)*N(3)
       plotCube(Cubes(i,:),i,N);
   end
   %axis([0 1.1 0 0.3 0 0.3])
return

function Next = Adjacent(N,D)  
    Next = zeros((N(1)+1)*(N(2)+1)*(N(3)+1),27);
    %1 - Itself
    %2 - Forward
    %3 - Backward
    %4 - Forward + Left
    %5 - Forward + Right
    %6 - Left
    %7 - Right
    %8 - Backward + Left
    %9 - Backward + Right
    for i = 1:N(1)+1
       for j = 1:N(2)+1
          for k = 1:N(3)+1
             Next(D(i,j,k),1) = D(i,j,k);
             if(i<N(1)+1)
                Next(D(i,j,k),2) = D(i+1,j,k);
             else
                 Next(D(i,j,k),2) = nan;
             end
             if(i>1)
                Next(D(i,j,k),3) = D(i-1,j,k);
             else
                 Next(D(i,j,k),3) = nan;
             end
             if(i<N(1)+1 && j < N(2)+1)
                Next(D(i,j,k),4) = D(i+1,j+1,k);
             else
                 Next(D(i,j,k),4) = nan;
             end
             if(i<N(1)+1 && j > 1)
                Next(D(i,j,k),5) = D(i+1,j-1,k);
             else
                 Next(D(i,j,k),5) = nan;
             end
             if(j < N(2)+1)
                Next(D(i,j,k),6) = D(i,j+1,k);
             else
                 Next(D(i,j,k),6) = nan;
             end
             if(j>1)
                Next(D(i,j,k),7) = D(i,j-1,k);
             else
                 Next(D(i,j,k),7) = nan;
             end
             if(i>1 && j < N(2)+1)
                Next(D(i,j,k),8) = D(i-1,j+1,k);
             else
                 Next(D(i,j,k),8) = nan;
             end
             if(i>1 && j>1)
                Next(D(i,j,k),9) = D(i-1,j-1,k);
             else
                 Next(D(i,j,k),9) = nan;
             end
             
             if(k < N(3)+1)
                 Next(D(i,j,k),10) = D(i,j,k+1);
                 if(i<N(1)+1)
                    Next(D(i,j,k),11) = D(i+1,j,k+1);
                 else
                     Next(D(i,j,k),11) = nan;
                 end
                 if(i>1)
                    Next(D(i,j,k),12) = D(i-1,j,k+1);
                 else
                     Next(D(i,j,k),12) = nan;
                 end
                 if(i<N(1)+1 && j < N(2)+1)
                    Next(D(i,j,k),13) = D(i+1,j+1,k+1);
                 else
                     Next(D(i,j,k),13) = nan;
                 end
                 if(i<N(1)+1 && j > 1)
                    Next(D(i,j,k),14) = D(i+1,j-1,k+1);
                 else
                     Next(D(i,j,k),14) = nan;
                 end
                 if(j < N(2)+1)
                    Next(D(i,j,k),15) = D(i,j+1,k+1);
                 else
                     Next(D(i,j,k),15) = nan;
                 end
                 if(j>1)
                    Next(D(i,j,k),16) = D(i,j-1,k+1);
                 else
                     Next(D(i,j,k),16) = nan;
                 end
                 if(i>1 && j < N(2)+1)
                    Next(D(i,j,k),17) = D(i-1,j+1,k+1);
                 else
                     Next(D(i,j,k),17) = nan;
                 end
                 if(i>1 && j>1)
                    Next(D(i,j,k),18) = D(i-1,j-1,k+1);
                 else
                     Next(D(i,j,k),18) = nan;
                 end
             else
                 Next(D(i,j,k),10) = nan;
                 Next(D(i,j,k),11) = nan;
                 Next(D(i,j,k),12) = nan; 
                 Next(D(i,j,k),13) = nan;
                 Next(D(i,j,k),14) = nan;
                 Next(D(i,j,k),15) = nan;
                 Next(D(i,j,k),16) = nan;
                 Next(D(i,j,k),17) = nan;
                 Next(D(i,j,k),18) = nan;
             end
             
             if(k>1)
                 Next(D(i,j,k),19) = D(i,j,k-1);
                 if(i<N(1)+1)
                    Next(D(i,j,k),20) = D(i+1,j,k-1);
                 else
                     Next(D(i,j,k),20) = nan;
                 end
                 if(i>1)
                    Next(D(i,j,k),21) = D(i-1,j,k-1);
                 else
                     Next(D(i,j,k),21) = nan;
                 end
                 if(i<N(1)+1 && j < N(2)+1)
                    Next(D(i,j,k),22) = D(i+1,j+1,k-1);
                 else
                     Next(D(i,j,k),22) = nan;
                 end
                 if(i<N(1)+1 && j > 1)
                    Next(D(i,j,k),23) = D(i+1,j-1,k-1);
                 else
                     Next(D(i,j,k),23) = nan;
                 end
                 if(j < N(2)+1)
                    Next(D(i,j,k),24) = D(i,j+1,k-1);
                 else
                     Next(D(i,j,k),24) = nan;
                 end
                 if(j>1)
                    Next(D(i,j,k),25) = D(i,j-1,k-1);
                 else
                     Next(D(i,j,k),25) = nan;
                 end
                 if(i>1 && j < N(2)+1)
                    Next(D(i,j,k),26) = D(i-1,j+1,k-1);
                 else
                     Next(D(i,j,k),26) = nan;
                 end
                 if(i>1 && j>1)
                    Next(D(i,j,k),27) = D(i-1,j-1,k-1);
                 else
                     Next(D(i,j,k),27) = nan;
                 end
             else
                 Next(D(i,j,k),19) = nan;
                 Next(D(i,j,k),20) = nan;
                 Next(D(i,j,k),21) = nan; 
                 Next(D(i,j,k),22) = nan;
                 Next(D(i,j,k),23) = nan;
                 Next(D(i,j,k),24) = nan;
                 Next(D(i,j,k),25) = nan;
                 Next(D(i,j,k),26) = nan;
                 Next(D(i,j,k),27) = nan;
             end
          end
       end
    end
return

function B = AdjacentType()%CHECKED
    B = zeros(27,27,16);
    %CHECKED
    B(1,1,:) = [2 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(1,2,:) = [2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(1,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(1,4,:) = [2 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(1,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(1,6,:) = [2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(1,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(1,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(1,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(1,10,:) = [2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(1,11,:) = [2 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(1,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(1,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(1,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(1,15,:) = [2 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(1,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(1,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(1,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(1,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(1,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(1,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(1,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(1,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(1,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(1,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(1,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(1,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(2,1,:) = [1 1 2 2 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(2,2,:) = [2 3 1 4 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(2,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(2,4,:) = [2 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(2,5,:) = [1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(2,6,:) = [2 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(2,7,:) = [1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(2,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(2,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(2,10,:) = [2 6 1 5 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(2,11,:) = [2 7 1 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(2,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(2,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(2,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(2,15,:) = [2 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(2,16,:) = [1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(2,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(2,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(2,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(2,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(2,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(2,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(2,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(2,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(2,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(2,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(2,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(3,1,:) = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(3,2,:) = [1 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(3,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(3,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(3,5,:) = [1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(3,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(3,7,:) = [1 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(3,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(3,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(3,10,:) = [1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(3,11,:) = [1 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(3,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(3,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(3,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(3,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(3,16,:) = [1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(3,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(3,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(3,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(3,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(3,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(3,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(3,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(3,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(3,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(3,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(3,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(4,1,:) = [2 2 6 6 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(4,2,:) = [2 3 6 7 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(4,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(4,4,:) = [2 4 6 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(4,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(4,6,:) = [2 1 6 5 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(4,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(4,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(4,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(4,10,:) = [2 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(4,11,:) = [2 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(4,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(4,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(4,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(4,15,:) = [2 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(4,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(4,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(4,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(4,19,:) = [6 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(4,20,:) = [6 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(4,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(4,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(4,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(4,24,:) = [6 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(4,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(4,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(4,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(5,1,:) = [1 1 2 2 5 5 6 6 0 0 0 0 0 0 0 0]; %Itself
    B(5,2,:) = [2 3 1 4 5 8 6 7 0 0 0 0 0 0 0 0];%Forward
    B(5,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(5,4,:) = [2 4 6 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(5,5,:) = [1 3 5 7 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(5,6,:) = [2 1 6 5 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(5,7,:) = [1 2 5 6 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(5,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(5,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(5,10,:) = [2 6 1 5 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(5,11,:) = [2 7 1 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(5,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(5,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(5,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(5,15,:) = [2 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(5,16,:) = [1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(5,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(5,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(5,19,:) = [5 1 6 2 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(5,20,:) = [5 4 6 3 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(5,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(5,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(5,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(5,24,:) = [6 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(5,25,:) = [5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(5,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(5,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(6,1,:) = [1 1 5 5 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(6,2,:) = [1 4 5 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(6,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(6,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(6,5,:) = [1 3 5 7 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(6,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(6,7,:) = [1 2 5 6 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(6,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(6,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(6,10,:) = [1 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(6,11,:) = [1 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(6,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(6,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(6,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(6,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(6,16,:) = [1 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(6,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(6,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(6,19,:) = [5 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(6,20,:) = [5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(6,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(6,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(6,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(6,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(6,25,:) = [5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(6,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(6,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(7,1,:) = [6 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(7,2,:) = [6 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(7,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(7,4,:) = [6 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(7,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(7,6,:) = [6 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(7,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(7,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(7,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(7,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(7,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(7,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(7,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(7,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(7,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(7,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(7,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(7,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(7,19,:) = [6 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(7,20,:) = [6 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(7,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(7,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(7,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(7,24,:) = [6 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(7,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(7,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(7,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(8,1,:) = [5 5 6 6 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(8,2,:) = [6 7 5 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(8,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(8,4,:) = [6 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(8,5,:) = [5 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(8,6,:) = [6 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(8,7,:) = [5 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(8,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(8,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(8,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(8,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(8,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(8,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(8,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(8,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(8,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(8,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(8,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(8,19,:) = [5 1 6 2 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(8,20,:) = [5 4 6 3 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(8,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(8,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(8,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(8,24,:) = [6 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(8,25,:) = [5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(8,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(8,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(9,1,:) = [5 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(9,2,:) = [5 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(9,3,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(9,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(9,5,:) = [5 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(9,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(9,7,:) = [5 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(9,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(9,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(9,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(9,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(9,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(9,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(9,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(9,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(9,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(9,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(9,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(9,19,:) = [5 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(9,20,:) = [5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(9,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(9,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(9,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(9,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(9,25,:) = [5 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(9,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(9,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    %CHECKED
    B(10,1,:) = [2 2 3 3 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(10,2,:) = [2 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(10,3,:) = [3 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(10,4,:) = [2 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(10,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(10,6,:) = [2 1 3 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(10,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(10,8,:) = [3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(10,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(10,10,:) = [2 6 3 7 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(10,11,:) = [2 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(10,12,:) = [3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(10,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(10,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(10,15,:) = [2 5 3 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(10,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(10,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(10,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(10,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(10,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(10,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(10,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(10,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(10,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(10,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(10,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(10,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(11,1,:) = [1 1 2 2 3 3 4 4 0 0 0 0 0 0 0 0]; %Itself
    B(11,2,:) = [2 3 1 4 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(11,3,:) = [3 2 4 1 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(11,4,:) = [2 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(11,5,:) = [1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(11,6,:) = [2 1 3 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(11,7,:) = [1 2 4 3 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(11,8,:) = [3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(11,9,:) = [4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(11,10,:) = [2 6 1 5 3 7 4 8 0 0 0 0 0 0 0 0];%Itself
    B(11,11,:) = [2 7 1 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(11,12,:) = [3 6 4 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(11,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(11,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(11,15,:) = [2 5 3 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(11,16,:) = [1 6 4 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(11,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(11,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(11,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(11,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(11,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(11,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(11,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(11,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(11,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(11,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(11,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(12,1,:) = [1 1 4 4 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(12,2,:) = [1 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(12,3,:) = [4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(12,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(12,5,:) = [1 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(12,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(12,7,:) = [1 2 4 3 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(12,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(12,9,:) = [4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(12,10,:) = [1 5 4 8 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(12,11,:) = [1 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(12,12,:) = [4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(12,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(12,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(12,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(12,16,:) = [1 6 4 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(12,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(12,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(12,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(12,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(12,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(12,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(12,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(12,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(12,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(12,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(12,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(13,1,:) = [2 2 3 3 6 6 7 7 0 0 0 0 0 0 0 0]; %Itself
    B(13,2,:) = [2 3 6 7 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(13,3,:) = [3 2 7 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(13,4,:) = [2 4 6 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(13,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(13,6,:) = [2 1 6 5 3 4 7 8 0 0 0 0 0 0 0 0];%Left
    B(13,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(13,8,:) = [3 1 7 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(13,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(13,10,:) = [2 6 3 7 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(13,11,:) = [2 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(13,12,:) = [3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(13,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(13,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(13,15,:) = [2 5 3 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(13,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(13,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(13,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(13,19,:) = [6 2 7 3 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(13,20,:) = [6 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(13,21,:) = [7 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(13,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(13,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(13,24,:) = [6 1 7 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(13,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(13,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(13,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    %CHECKED
    B(14,1,:) = [1 1 2 2 3 3 4 4 5 5 6 6 7 7 8 8]; %Itself
    B(14,2,:) = [1 4 2 3 5 8 6 7 0 0 0 0 0 0 0 0];%Forward
    B(14,3,:) = [4 1 3 2 8 5 7 6 0 0 0 0 0 0 0 0];%Backward
    B(14,4,:) = [2 4 6 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(14,5,:) = [1 3 5 7 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(14,6,:) = [2 1 6 5 3 4 7 8 0 0 0 0 0 0 0 0];%Left
    B(14,7,:) = [1 2 5 6 4 3 8 7 0 0 0 0 0 0 0 0];%Right
    B(14,8,:) = [3 1 7 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(14,9,:) = [4 2 8 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up%CHECKED
    B(14,10,:) = [1 5 2 6 3 7 4 8 0 0 0 0 0 0 0 0];%Itself
    B(14,11,:) = [2 7 1 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(14,12,:) = [3 6 4 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(14,13,:) = [2 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(14,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(14,15,:) = [2 5 3 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(14,16,:) = [1 6 4 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(14,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(14,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down%CHECKED
    B(14,19,:) = [5 1 6 2 7 3 8 4 0 0 0 0 0 0 0 0];%Itself
    B(14,20,:) = [5 4 6 3 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(14,21,:) = [8 1 7 2 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(14,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(14,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(14,24,:) = [6 1 7 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(14,25,:) = [5 2 8 3 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(14,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(14,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(15,1,:) = [1 1 4 4 5 5 8 8 0 0 0 0 0 0 0 0]; %Itself
    B(15,2,:) = [1 4 5 8 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(15,3,:) = [4 1 8 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(15,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(15,5,:) = [1 3 5 7 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(15,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(15,7,:) = [1 2 5 6 4 3 8 7 0 0 0 0 0 0 0 0];%Right
    B(15,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(15,9,:) = [4 2 8 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(15,10,:) = [1 5 4 8 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(15,11,:) = [1 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(15,12,:) = [4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(15,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(15,14,:) = [1 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(15,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(15,16,:) = [1 6 4 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(15,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(15,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(15,19,:) = [5 1 8 4 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(15,20,:) = [5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(15,21,:) = [8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(15,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(15,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(15,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(15,25,:) = [5 2 8 3 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(15,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(15,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(16,1,:) = [6 6 7 7 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(16,2,:) = [6 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(16,3,:) = [7 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(16,4,:) = [6 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(16,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(16,6,:) = [6 5 7 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(16,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(16,8,:) = [7 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(16,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(16,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(16,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(16,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(16,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(16,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(16,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(16,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(16,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(16,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(16,19,:) = [6 2 7 3 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(16,20,:) = [6 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(16,21,:) = [7 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(16,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(16,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(16,24,:) = [6 1 7 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(16,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(16,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(16,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(17,1,:) = [5 5 6 6 7 7 8 8 0 0 0 0 0 0 0 0]; %Itself
    B(17,2,:) = [5 8 6 7 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(17,3,:) = [8 5 7 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(17,4,:) = [6 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(17,5,:) = [5 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(17,6,:) = [6 5 7 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(17,7,:) = [5 6 8 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(17,8,:) = [7 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(17,9,:) = [8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(17,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(17,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(17,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(17,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(17,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(17,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(17,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(17,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(17,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(17,19,:) = [5 1 6 2 7 3 8 4 0 0 0 0 0 0 0 0];%Itself
    B(17,20,:) = [5 4 6 3 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(17,21,:) = [8 1 7 2 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(17,22,:) = [6 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(17,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(17,24,:) = [6 1 7 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(17,25,:) = [5 2 8 3 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(17,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(17,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(18,1,:) = [5 5 8 8 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(18,2,:) = [5 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(18,3,:) = [8 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(18,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(18,5,:) = [5 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(18,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(18,7,:) = [5 6 8 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(18,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(18,9,:) = [8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(18,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(18,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(18,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(18,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(18,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(18,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(18,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(18,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(18,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(18,19,:) = [5 1 8 4 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(18,20,:) = [5 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(18,21,:) = [8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(18,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(18,23,:) = [5 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(18,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(18,25,:) = [5 2 8 3 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(18,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(18,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   
    B(19,1,:) = [3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(19,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(19,3,:) = [3 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(19,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(19,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(19,6,:) = [3 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(19,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(19,8,:) = [3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(19,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(19,10,:) = [3 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(19,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(19,12,:) = [3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(19,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(19,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(19,15,:) = [3 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(19,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(19,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(19,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(19,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(19,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(19,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(19,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(19,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(19,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(19,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(19,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(19,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    
    B(20,1,:) = [3 3 4 4 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(20,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(20,3,:) = [3 2 4 1 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(20,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(20,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(20,6,:) = [3 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(20,7,:) = [4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(20,8,:) = [3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(20,9,:) = [4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(20,10,:) = [4 8 3 7 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(20,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(20,12,:) = [3 6 4 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(20,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(20,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(20,15,:) = [3 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(20,16,:) = [4 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(20,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(20,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(20,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(20,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(20,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(20,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(20,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(20,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(20,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(20,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(20,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    
    B(21,1,:) = [4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(21,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(21,3,:) = [4 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(21,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(21,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(21,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(21,7,:) = [4 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(21,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(21,9,:) = [4 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(21,10,:) = [4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(21,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(21,12,:) = [4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(21,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(21,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(21,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(21,16,:) = [4 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(21,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(21,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(21,19,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(21,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(21,21,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(21,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(21,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(21,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(21,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(21,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(21,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    
    B(22,1,:) = [3 3 7 7 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(22,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(22,3,:) = [3 2 7 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(22,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(22,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(22,6,:) = [7 8 3 4 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(22,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(22,8,:) = [3 1 7 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(22,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(22,10,:) = [3 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(22,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(22,12,:) = [3 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(22,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(22,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(22,15,:) = [3 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(22,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(22,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(22,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(22,19,:) = [7 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(22,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(22,21,:) = [7 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(22,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(22,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(22,24,:) = [7 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(22,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(22,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(22,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................

    B(23,1,:) = [3 3 4 4 7 7 8 8 0 0 0 0 0 0 0 0]; %Itself
    B(23,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(23,3,:) = [4 1 3 2 8 5 7 6 0 0 0 0 0 0 0 0];%Backward
    B(23,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(23,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(23,6,:) = [3 4 7 8 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(23,7,:) = [4 3 8 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(23,8,:) = [3 1 7 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(23,9,:) = [4 2 8 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(23,10,:) = [4 8 3 7 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(23,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(23,12,:) = [3 6 4 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(23,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(23,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(23,15,:) = [3 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(23,16,:) = [4 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(23,17,:) = [3 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(23,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(23,19,:) = [8 4 7 3 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(23,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(23,21,:) = [8 1 7 2 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(23,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(23,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(23,24,:) = [7 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(23,25,:) = [8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(23,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(23,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(24,1,:) = [4 4 8 8 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(24,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(24,3,:) = [4 1 8 5 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(24,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(24,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(24,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(24,7,:) = [4 3 8 7 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(24,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(24,9,:) = [4 2 8 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(24,10,:) = [4 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(24,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(24,12,:) = [4 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(24,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(24,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(24,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(24,16,:) = [4 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(24,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(24,18,:) = [4 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(24,19,:) = [8 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(24,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(24,21,:) = [8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(24,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(24,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(24,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(24,25,:) = [8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(24,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(24,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(25,1,:) = [7 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(25,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Forward
    B(25,3,:) = [7 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(25,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(25,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(25,6,:) = [7 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(25,7,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(25,8,:) = [7 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(25,9,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(25,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(25,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(25,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(25,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(25,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(25,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(25,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(25,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(25,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(25,19,:) = [7 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(25,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(25,21,:) = [7 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(25,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(25,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(25,24,:) = [7 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(25,25,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(25,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(25,27,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(26,1,:) = [7 7 8 8 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(26,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(26,3,:) = [8 5 7 6 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(26,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(26,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(26,6,:) = [7 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(26,7,:) = [8 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(26,8,:) = [7 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(26,9,:) = [8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(26,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(26,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(26,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(26,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(26,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(26,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(26,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(26,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(26,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(26,19,:) = [7 3 8 4 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(26,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(26,21,:) = [8 1 7 2 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(26,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(26,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(26,24,:) = [7 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(26,25,:) = [8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(26,26,:) = [7 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(26,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %.........................................................................................................
    B(27,1,:) = [8 8 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Itself
    B(27,2,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(27,3,:) = [8 5 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(27,4,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(27,5,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(27,6,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(27,7,:) = [8 7 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(27,8,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(27,9,:) = [8 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Up
    B(27,10,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(27,11,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(27,12,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(27,13,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(27,14,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(27,15,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(27,16,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(27,17,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(27,18,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
    %Down
    B(27,19,:) = [8 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Itself
    B(27,20,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward
    B(27,21,:) = [8 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward
    B(27,22,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Left
    B(27,23,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Forward + Right
    B(27,24,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Left
    B(27,25,:) = [8 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Right
    B(27,26,:) = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Left
    B(27,27,:) = [8 2 0 0 0 0 0 0 0 0 0 0 0 0 0 0];%Backward + Right
return

function [K11,K12,K13,K22,K23,K33,M] = SmallMatrix(Delta)
syms x
syms y
syms z
Q = [1 x y z x*y x*z y*z x*y*z];

K11 = zeros(size(Q,2),'single')*x*y*z;
K12 = zeros(size(Q,2),'single')*x*y*z;
K13 = zeros(size(Q,2),'single')*x*y*z;
K22 = zeros(size(Q,2),'single')*x*y*z;
K23 = zeros(size(Q,2),'single')*x*y*z;
K33 = zeros(size(Q,2),'single')*x*y*z;

M   = zeros(size(Q,2),'single')*x*y*z;

for i = 1:size(Q,2)
   for j = 1:size(Q,2)
      K11(j,i) = diff(Q(j),x)*diff(Q(i),x);
      K12(j,i) = diff(Q(j),y)*diff(Q(i),x);
      K13(j,i) = diff(Q(j),z)*diff(Q(i),x);
      K22(j,i) = diff(Q(j),y)*diff(Q(i),y);
      K23(j,i) = diff(Q(j),z)*diff(Q(i),y);
      K33(j,i) = diff(Q(j),z)*diff(Q(i),z);
      M(j,i)   = Q(j)*Q(i);
      
   end
end
K11 = int(int(int(K11,x,[0,1]),y,[0,1]),z,[0,1]);
K12 = int(int(int(K12,x,[0,1]),y,[0,1]),z,[0,1]);
K13 = int(int(int(K13,x,[0,1]),y,[0,1]),z,[0,1]);
K22 = int(int(int(K22,x,[0,1]),y,[0,1]),z,[0,1]);
K23 = int(int(int(K23,x,[0,1]),y,[0,1]),z,[0,1]);
K33 = int(int(int(K33,x,[0,1]),y,[0,1]),z,[0,1]);
M   = int(int(int(M,x,[0,1]),y,[0,1]),z,[0,1]);

T = MATRIX_T(Q);
Tinv = inv(T);

K11 = (Tinv)'*K11*Tinv;
K12 = (Tinv)'*K12*Tinv;
K13 = (Tinv)'*K13*Tinv;
K22 = (Tinv)'*K22*Tinv;
K23 = (Tinv)'*K23*Tinv;
K33 = (Tinv)'*K33*Tinv;
M   = (Tinv)'*M*Tinv;

K11 = double(K11*Delta(2)*Delta(3)/Delta(1));
K12 = double(K12*Delta(3));
K13 = double(K13*Delta(2));
K22 = double(K22*Delta(1)*Delta(3)/Delta(2));
K23 = double(K23*Delta(1));
K33 = double(K33*Delta(1)*Delta(2)/Delta(3));
M   = double(M*Delta(1)*Delta(2)*Delta(3));
return

function T = MATRIX_T(Q)
syms x;
syms y;
syms z;

n = size(Q,2);

T = zeros(n);
for j = 1:n
   T(j,1) = subs(Q(j),[x,y,z],[0,1,0]);
   T(j,2) = subs(Q(j),[x,y,z],[0,0,0]);
   T(j,3) = subs(Q(j),[x,y,z],[1,0,0]);
   T(j,4) = subs(Q(j),[x,y,z],[1,1,0]);
   T(j,5) = subs(Q(j),[x,y,z],[0,1,1]);
   T(j,6) = subs(Q(j),[x,y,z],[0,0,1]);
   T(j,7) = subs(Q(j),[x,y,z],[1,0,1]);
   T(j,8) = subs(Q(j),[x,y,z],[1,1,1]);
end
T = T';
return;

function [B11,B12,B13,B22,B23,B33,BM] = AddMatrix(K11,K12,K13,K22,K23,K33,M)
B = AdjacentType();
B11 = zeros(27,27);
B12 = zeros(27,27);
B13 = zeros(27,27);
B22 = zeros(27,27);
B23 = zeros(27,27);
B33 = zeros(27,27);
BM = zeros(27,27);
for i = 1:27
   for j = 1:27
      for k = 1:2:15
          if B(i,j,k) ~= 0
            B11(i,j) = B11(i,j) + K11(B(i,j,k),B(i,j,k+1));
            B12(i,j) = B12(i,j) + K12(B(i,j,k),B(i,j,k+1));
            B13(i,j) = B13(i,j) + K13(B(i,j,k),B(i,j,k+1));
            B22(i,j) = B22(i,j) + K22(B(i,j,k),B(i,j,k+1));
            B23(i,j) = B23(i,j) + K23(B(i,j,k),B(i,j,k+1));
            B33(i,j) = B33(i,j) + K33(B(i,j,k),B(i,j,k+1));
            BM(i,j) = BM(i,j) + M(B(i,j,k),B(i,j,k+1));
          else
              break
          end
      end
   end
end
return;

function [K11,K12,K13,K22,K23,K33,M,D,E] = Matrices(Delta,N,mwb,method)
[K11q,K12q,K13q,K22q,K23q,K33q,Mq] = SmallMatrix(Delta);
[B11,B12,B13,B22,B23,B33,BM] = AddMatrix(K11q,K12q,K13q,K22q,K23q,K33q,Mq);
%{
K11 = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
K12 = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
K13 = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
K22 = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
K23 = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
K33 = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
M = spalloc((N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1)*27);
%}
[D,E] = Domain(N,Delta);
A = Adjacent(N,D);
T = Type(N,A); 
%

mwb.Update(2, 1, 0, ['Matrix Creation ' num2str(0) '%']);
%{
for i = 1:(N(1)+1)*(N(2)+1)*(N(3)+1)
    mwb.Update(2, 1, i/((N(1)+1)*(N(2)+1)*(N(3)+1)+1), ['Matrix Creation ' num2str(i/((N(1)+1)*(N(2)+1)*(N(3)+1)+1)*100) '%']);
    for j = 1:27
       k = 1;
       while (B(T(i),j,k) ~= 0 && ~isnan(A(i,j)))
            K11(i,A(i,j)) = K11(i,A(i,j)) + K11q(B(T(i),j,k),B(T(i),j,k+1));
            K12(i,A(i,j)) = K12(i,A(i,j)) + K12q(B(T(i),j,k),B(T(i),j,k+1));
            K13(i,A(i,j)) = K13(i,A(i,j)) + K13q(B(T(i),j,k),B(T(i),j,k+1));
            K22(i,A(i,j)) = K22(i,A(i,j)) + K22q(B(T(i),j,k),B(T(i),j,k+1));
            K23(i,A(i,j)) = K23(i,A(i,j)) + K23q(B(T(i),j,k),B(T(i),j,k+1));
            K33(i,A(i,j)) = K33(i,A(i,j)) + K33q(B(T(i),j,k),B(T(i),j,k+1));
          
            M(i,A(i,j)) = M(i,A(i,j)) + Mq(B(T(i),j,k),B(T(i),j,k+1));
            k = k + 2;
            if(k >= 16)
               break; 
            end
       end
    end
end
%}
n = 0;
for i = 1:(N(1)+1)*(N(2)+1)*(N(3)+1)
    for j = 1:27
      if ~isnan(A(i,j))
         n = n +1; 
      end
    end
end

%{

%}
%{

NAN_A = ~isnan(A);
A2 = repmat(A(:,1),1,size(A,2));
A3 = repmat([1:27]',1,size(A,1))';
iy = A(NAN_A);
ix = A2(NAN_A);
iz = A3(NAN_A);
[ix iy iz];

for i = 1:size(iy,1)
    mwb.Update(2, 1, i/(size(iy,1)), ['Matrix Creation ' num2str(i/size(iy,1)*100) '%']);
    k = 1;
    while(k < 16 && B(T(ix(i)),iz(i),k) ~= 0)
        K11s(i) = K11s(i) + K11q(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1));
        K12s(i) = K12s(i) + K12q(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1)); 
        K13s(i) = K13s(i) + K13q(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1)); 
        K22s(i) = K22s(i) + K22q(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1)); 
        K23s(i) = K23s(i) + K23q(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1)); 
        K33s(i) = K33s(i) + K33q(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1));
        Ms(i) = Ms(i) + Mq(B(T(ix(i)),iz(i),k),B(T(ix(i)),iz(i),k+1)); 
        k = k + 2;
    end
end
%}
%%{

if (method == 2)
    mwb.Update(2, 1, 0.1, ['Matrix Creation ' num2str(10) '%']);
    InvA = A';
    NAN_A = ~isnan(InvA);
    A2 = repmat([1:27]',1,size(InvA,2))';
    A3 = repmat(1:size(A,1),size(A,2),1)';
    InvA3 = A3';
    InvA2 = A2';
    iy = InvA(NAN_A);
    ix = InvA3(NAN_A);
    iz = InvA2(NAN_A);

    Typex = T(ix);

    %BAdd = B(Typex,Typey,:);
mwb.Update(2, 1, 0.2, ['Matrix Creation ' num2str(20) '%']);
    K11s = B11(sub2ind(size(B11),Typex,iz));
    mwb.Update(2, 1, 0.3, ['Matrix Creation ' num2str(30) '%']);
    K12s = B12(sub2ind(size(B12),Typex,iz));
    mwb.Update(2, 1, 0.4, ['Matrix Creation ' num2str(40) '%']);
    K13s = B13(sub2ind(size(B13),Typex,iz));
    mwb.Update(2, 1, 0.5, ['Matrix Creation ' num2str(50) '%']);
    K22s = B22(sub2ind(size(B22),Typex,iz));
    mwb.Update(2, 1, 0.6, ['Matrix Creation ' num2str(60) '%']);
    K23s = B23(sub2ind(size(B23),Typex,iz));
    mwb.Update(2, 1, 0.7, ['Matrix Creation ' num2str(70) '%']);
    K33s = B33(sub2ind(size(B33),Typex,iz));
    mwb.Update(2, 1, 0.8, ['Matrix Creation ' num2str(80) '%']);
    Ms = BM(sub2ind(size(BM),Typex,iz));
    mwb.Update(2, 1, 0.9, ['Matrix Creation ' num2str(90) '%']);
elseif (method == 1)
    B = AdjacentType();
    K11s = zeros(n,1);
    K12s = zeros(n,1);
    K13s = zeros(n,1);
    K22s = zeros(n,1);
    K23s = zeros(n,1);
    K33s = zeros(n,1);
    Ms = zeros(n,1);
    ix = zeros(n,1);
    iy = zeros(n,1);
    ii = 1;
    for i = 1:(N(1)+1)*(N(2)+1)*(N(3)+1)
        mwb.Update(2, 1, i/((N(1)+1)*(N(2)+1)*(N(3)+1)+1), ['Matrix Creation ' num2str(i/((N(1)+1)*(N(2)+1)*(N(3)+1)+1)*100) '%']);
        for j = 1:27
          if ~isnan(A(i,j))
             ix(ii) = A(i,1);
             iy(ii) = A(i,j);
             k = 1;
             while (B(T(i),j,k) ~= 0 && ~isnan(A(i,j)))

                K11s(ii) = K11s(ii) + K11q(B(T(i),j,k),B(T(i),j,k+1));
                K12s(ii) = K12s(ii) + K12q(B(T(i),j,k),B(T(i),j,k+1)); 
                K13s(ii) = K13s(ii) + K13q(B(T(i),j,k),B(T(i),j,k+1)); 
                K22s(ii) = K22s(ii) + K22q(B(T(i),j,k),B(T(i),j,k+1)); 
                K23s(ii) = K23s(ii) + K23q(B(T(i),j,k),B(T(i),j,k+1)); 
                K33s(ii) = K33s(ii) + K33q(B(T(i),j,k),B(T(i),j,k+1)); 
                Ms(ii) = Ms(ii) + Mq(B(T(i),j,k),B(T(i),j,k+1)); 
                k = k + 2;
                if(k >= 16)
                   break; 
                end
             end
             Tempi(ii) = T(i);
                 Tempj(ii) = j;
             ii = ii +1;
          end
        end
    end
end

%}
K11 = sparse(ix,iy,K11s,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
K12 = sparse(ix,iy,K12s,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
K13 = sparse(ix,iy,K13s,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
K22 = sparse(ix,iy,K22s,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
K23 = sparse(ix,iy,K23s,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
K33 = sparse(ix,iy,K33s,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
M = sparse(ix,iy,Ms,(N(1)+1)*(N(2)+1)*(N(3)+1),(N(1)+1)*(N(2)+1)*(N(3)+1));
%{
B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,2:2:16)
K11q(B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,1:2:15),B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,2:2:16))

sum(K11q(B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,1:2:15)>0,B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,2:2:16))>0)
K11(1:(N(1)+1)*(N(2)+1)*(N(3)+1),1:(N(1)+1)*(N(2)+1)*(N(3)+1)) = sum(B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,1:2:15))
K11(1,:)
size(K11)
K11(1:(N(1)+1)*(N(2)+1)*(N(3)+1),A(~isnan(A(1:(N(1)+1)*(N(2)+1)*(N(3)+1),1:27)))) =  sum(K11q(B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,1:2:15)>0,B(T(1:(N(1)+1)*(N(2)+1)*(N(3)+1)),1:27,2:2:1))>0);
%}

mwb.Update(2, 1, 1, ['Matrix Creation ' num2str(100) '%']);
return;

function T = Type(N,A)
T = zeros((N(1)+1)*(N(2)+1)*(N(3)+1),1);
TEST = [1    10   nan    11   nan     2   nan   nan   nan     4    13   nan    14   nan     5   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        2    11   nan    12    10     3     1   nan   nan     5    14   nan    15    13     6     4   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        3    12   nan   nan    11   nan     2   nan   nan     6    15   nan   nan    14   nan     5   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        4    13   nan    14   nan     5   nan   nan   nan     7    16   nan    17   nan     8   nan   nan   nan     1    10   nan    11   nan     2   nan   nan   nan;
        5    14   nan    15    13     6     4   nan   nan     8    17   nan    18    16     9     7   nan   nan     2    11   nan    12    10     3     1   nan   nan;
        6    15   nan   nan    14   nan     5   nan   nan     9    18   nan   nan    17   nan     8   nan   nan     3    12   nan   nan    11   nan     2   nan   nan;
        7    16   nan    17   nan     8   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan     4    13   nan    14   nan     5   nan   nan   nan;
        8    17   nan    18    16     9     7   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan     5    14   nan    15    13     6     4   nan   nan;
        9    18   nan   nan    17   nan     8   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan     6    15   nan   nan    14   nan     5   nan   nan;
        10    19     1    20   nan    11   nan     2   nan    13    22     4    23   nan    14   nan     5   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        11    20     2    21    19    12    10     3     1    14    23     5    24    22    15    13     6     4   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        12    21     3   nan    20   nan    11   nan     2    15    24     6   nan    23   nan    14   nan     5   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        13    22     4    23   nan    14   nan     5   nan    16    25     7    26   nan    17   nan     8   nan    10    19     1    20   nan    11   nan     2   nan;
        14    23     5    24    22    15    13     6     4    17    26     8    27    25    18    16     9     7    11    20     2    21    19    12    10     3     1;
        15    24     6   nan    23   nan    14   nan     5    18    27     9   nan    26   nan    17   nan     8    12    21     3   nan    20   nan    11   nan     2;
        16    25     7    26   nan    17   nan     8   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan    13    22     4    23   nan    14   nan     5   nan;
        17    26     8    27    25    18    16     9     7   nan   nan   nan   nan   nan   nan   nan   nan   nan    14    23     5    24    22    15    13     6     4;
        18    27     9   nan    26   nan    17   nan     8   nan   nan   nan   nan   nan   nan   nan   nan   nan    15    24     6   nan    23   nan    14   nan     5;
        19   nan    10   nan   nan    20   nan    11   nan    22   nan    13   nan   nan    23   nan    14   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        20   nan    11   nan   nan    21    19    12    10    23   nan    14   nan   nan    24    22    15    13   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        21   nan    12   nan   nan   nan    20   nan    11    24   nan    15   nan   nan   nan    23   nan    14   nan   nan   nan   nan   nan   nan   nan   nan   nan;
        22   nan    13   nan   nan    23   nan    14   nan    25   nan    16   nan   nan    26   nan    17   nan    19   nan    10   nan   nan    20   nan    11   nan;
        23   nan    14   nan   nan    24    22    15    13    26   nan    17   nan   nan    27    25    18    16    20   nan    11   nan   nan    21    19    12    10;
        24   nan    15   nan   nan   nan    23   nan    14    27   nan    18   nan   nan   nan    26   nan    17    21   nan    12   nan   nan   nan    20   nan    11;
        25   nan    16   nan   nan    26   nan    17   nan   nan   nan   nan   nan   nan   nan   nan   nan   nan    22   nan    13   nan   nan    23   nan    14   nan;
        26   nan    17   nan   nan    27    25    18    16   nan   nan   nan   nan   nan   nan   nan   nan   nan    23   nan    14   nan   nan    24    22    15    13;
        27   nan    18   nan   nan   nan    26   nan    17   nan   nan   nan   nan   nan   nan   nan   nan   nan    24   nan    15   nan   nan   nan    23   nan    14];
for i = 1:(N(1)+1)*(N(2)+1)*(N(3)+1)
   for j = 1:27
      bflag = true;
      for k = 1:27
         if(isnan(A(i,k))~= isnan(TEST(j,k)))
            bflag = false; 
         end
      end
      if(bflag == true)
         T(i) = j;
         break;
      end
   end
end
return