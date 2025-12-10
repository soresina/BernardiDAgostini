%%%%%%%%%%%%%%%%%%%%
% PDE model simulation for the interaction between cells and medium

% time
t0 = 0;
tend = 2500;
dt = 0.05;
time = t0:dt:tend;
tdim = length(time);

% space
x0 = 0;
xend = 4;
dx = 0.05;
xall = x0:dx:xend;
xdim = length(xall);
% the space domain is a square [0,xend]x[0,xend]

% diffusion parameters
Dn = 1.8e-6;
Dmu = 2.36;
%DG = 0.236;    % value estimated with Eistein-Stokes
DG=0;           % assumption

% arbitrary parameters
Vol=0.7;
% mu0,G0 in [0,1]
mu0 = 0.25;    % RGM: 0.25 ca.
G0 = 0.0001;   % RGM: 0.1

% convenient choice for the last parameters
S=1;
a_G=0;         % assumption
a_mu = 0.8;

% experimentally fitted parameters
a = 2.2;
h1 = 3.2;
h2 = 11.2;
mn=0.013;

% random seed
rng(2);

% reaction functions
f1 =@(n,mu, G) (a .* G.*mu) ./ ((h1.*mu +1).*( h2.*G + 1)) .*n - mn.*n;
f2 =@(n,mu, G) (-a_mu .*mu .*n .*G) ./(h1.*mu+1) + Vol.*(mu0-mu);
f3 =@(n,G) (-a_G .* G.*n) ./(h2.*G+1) + S.*n + Vol.*(G0-G);


%% matrices for the implicit scheme

[X,Y] = meshgrid(0:dx:xend, 0:dx:xend);
N= xdim -1;         % number of intervals in the space doamin
M = N;              % For a squared space domain
dy=dx;

% auxiliary matrices
R= zeros (M+1,M+1);R(1 ,1) = -1; R(1 ,2) =1; R(M+1,M)=1; R(M+1,M+1) = -1; 
i=2; j=1;
while i <=M && j <=M
    R(i,j)=1; R(i,j+1) = -2; R(i,j+2) =1;
    i=i+1; j=j+1;
end
R =1/( dx)^2*R;

S =1/( dy)^2* eye(M+1);
T=-S;
V=2*T;

% final matrix
A= zeros ((M+1) *(N+1) , (M+1) *(N+1));
A(1:M+1 ,1:M+1)=R+T; A(1:M+1, M+2:M+2+M)=S;
i=M+2; j=1;
while i <=(M+1) *(N) && j <=(M+1) *(N+1)
    A(i:i+M, j:j+M)=S;
    A(i:i+M, j+M+1:j+M+1+M)=R+V;
    A(i:i+M,j+M+1+M+1:j+M+1+M+1+M)=S;
    i=i+M+1; j=j+M+1;
end
A((M+1) *(N+1) -M:(M+1) *(N+1) , (M+1) *(N+1) -M:(M+1) *(N+1))=R+T;
A((M+1) *(N+1) -M:(M+1) *(N+1) , (M+1) *(N+1) -M-M -1:( M+1) *(N+1) -M -1)=S;
A= sparse (A);

%% numerical scheme

% initial conditions
n0 = 2;
sd = 1;
noise = randi([-sd,sd], (M+1) *(N+1),1);
n = n0 + noise;
mu = mu0 * ones((M+1) *(N+1),1);
G = G0 * ones((M+1) *(N+1),1);

n1 = n0;
mu1 = mu0;
G1 = G0;
k=1;        % used to distinguish for the first iteration
z=0;        % time of the simulation
nt = zeros(tdim,1);     % contains the space integral at each time step


% for the animation
figure(10);
parN = n;
Nsup = reshape(parN,[M+1,N+1])';

h = surf(X,Y,Nsup);
shading interp
view(2)
axis equal

clim([0 max(Nsup(:))]);
clim([0 3.5]);      % per tenere la barra fissa, non so se è molto bello
hcb = colorbar;
hcb.Label.String = 'cell density n(x,y;t)';

xlabel('X'); ylabel('Y');
title(sprintf('n(x,t) at t = %.2f, V = %.2f', z, Vol));
drawnow;
step = 0;

while z <=tend
    F1 = f1(n,mu,G);
    F2 = f2(n,mu,G);
    F3 = f3(n,G);
    F1_1 = f1(n1,mu1,G1);
    F2_1 = f2(n1,mu1,G1);
    F3_1 = f3(n1,G1);

    n2 = n;
    mu2 = mu;
    G2 = G;

    if k==1     % the first iteration
        n = (speye ((M+1) *(N+1))-Dn*dt*A)\(n+dt* F1);
        mu = (speye ((M+1) *(N+1))-Dmu*dt*A)\(mu+dt* F2);
        G = (speye ((M+1) *(N+1))-DG*dt*A)\(G+dt* F3);
        k=2;
    else
        n = (3*speye ((M+1) *(N+1))-2*Dn*dt*A)\(4*n-n1+2*dt* (2*F1-F1_1));
        mu = (3*speye ((M+1) *(N+1))-2*Dmu*dt*A)\(4*mu-mu1+2*dt* (2*F2-F2_1));
        G = (3*speye ((M+1) *(N+1))-2*DG*dt*A)\(4*G-G1+2*dt* (2*F3-F3_1));
        k=k+1;
    end
    n1=n2;
    mu1=mu2;
    G1=G2;

    % integral over space
    nt_tot = sum(n(:,:));
    nt(k) = nt_tot;
    
    % update the plot 
    step = step + 1;
    if mod(step,150) == 0
        parN = n;
        Nsup = reshape(n,[M+1,N+1])';

        % update surface
        set(h,'ZData',Nsup);

        % update title
        title(sprintf('n(x,y;t) at t = %.2f, V = %.2f', z, Vol));

        drawnow;
    end


    z=z+dt;
end

%% integral over space of n

figure()
hold on
box on
title('n(t) total over space','FontSize',18)
xlabel('Time[h]','FontSize',16)
ylabel('n(t)','FontSize',16)
plot(time,nt(1:length(time)), 'LineWidth',2.5)
ylim([0,max(nt)])
ax = gca;
ax.FontSize = 16; 
hold off



