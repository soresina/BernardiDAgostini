%%%%%%%%%%%%%%%%%%%%
% PDE model simulation for the interaction between cells and medium

% method:
explicit = 1
implicit = 0


% time
t0 = 0;
tend = 2500;
dt = 0.005;
time = t0:dt:tend;
tdim = length(time);

% space
x0 = 0;
xend = 4;
dx = 0.2;
xall = x0:dx:xend;
xdim = length(xall);

% diffusion parameters
Dn = 1.8e-6;
Dmu = 2.36;
%DG = 0.236;    % value estimated with Eistein-Stokes
DG=1e-8;           % assumption

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

% variables of the model
nxt = zeros(length(xall),length(time));
Gxt = zeros(length(xall),length(time));
muxt = zeros(length(xall),length(time));

% initial conditions
n0 = 1;
sd = 1;
noise = randi([-sd,sd], 1, length(xall));
for i = 1:length(xall)
    nxt(i,1) = n0 + noise(i);
    Gxt(i,1) = G0;
    muxt(i,1) = mu0;
end


%% explicit scheme
if explicit
    fprintf('stability dt<0.5*dx^2/Dn: dt=%.2f, 1/2*dx^2/Dn=%.2f \n', dt, 1/2*dx^2/Dn)
    fprintf('stability dt<0.5*dx^2/Dmu: dt=%.3f, 1/2*dx^2/Dmu=%.4f \n', dt, 0.5*dx^2/Dmu)
    if dt>0.5*dx^2/Dmu
        error('Error: CFL condition not verified');

    else
        for i= 2:tdim
            %nxt(i,1) = n0;%+noise(i);
            nxt(1,i) = nxt(1,i-1) + dt*f1(nxt(1,i-1), muxt(1,i-1), Gxt(1, i-1))+ 2*dt *Dn/(dx^2) * (nxt(2,i-1)-nxt(1,i-1)); % Neumann BC
            muxt(1,i) = muxt(1,i-1) + dt*f2(nxt(1,i-1), muxt(1,i-1), Gxt(1,i-1))+ 2*dt *Dmu/(dx^2) * (muxt(2,i-1)-muxt(1,i-1));
            Gxt(1,i) = Gxt(1,i-1) + dt*f3(nxt(1,i-1), Gxt(1,i-1))+ 2*dt *DG/(dx^2) * (Gxt(2,i-1)-Gxt(1,i-1));
            for j =2:xdim-1
                nxt(j,i) = nxt(j,i-1) + dt*f1(nxt(j,i-1), muxt(j,i-1), Gxt(j,i-1)) + dt* Dn/(dx^2) * (nxt(j-1,i-1)-2*nxt(j,i-1)+nxt(j+1,i-1));
                muxt(j,i) = muxt(j,i-1) + dt*f2(nxt(j,i-1), muxt(j,i-1), Gxt(j,i-1)) + dt* Dmu/(dx^2) * (muxt(j-1,i-1)-2*muxt(j,i-1)+muxt(j+1,i-1));
                Gxt(j,i) = Gxt(j,i-1) + dt*f3(nxt(j,i-1), Gxt(j,i-1)) + dt* DG/(dx^2) * (Gxt(j-1,i-1)-2*Gxt(j,i-1)+Gxt(j+1,i-1));
                if isnan(nxt(j,i))
                    x=j
                    t=i
                    error("some value for n is NaN");
                end
            end
            nxt(xdim,i) = nxt(xdim,i-1) + dt*f1(nxt(xdim,i-1), muxt(xdim,i-1), Gxt(xdim, i-1)) + 2*dt* Dn/(dx^2) * (nxt(xdim-1,i-1) - nxt(xdim, i-1));
            muxt(xdim,i) = muxt(xdim,i-1) + dt*f2(nxt(xdim,i-1), muxt(xdim,i-1),Gxt(xdim, i-1)) + 2*dt* Dmu/(dx^2) * (muxt(xdim-1,i-1) - muxt(xdim, i-1));
            Gxt(xdim,i) = Gxt(xdim,i-1) + dt*f3(nxt(xdim,i-1), Gxt(xdim, i-1)) + 2*dt* DG/(dx^2) * (Gxt(xdim-1,i-1) - Gxt(xdim, i-1));
        end
    end
end


%% implicit scheme

if implicit
    % matrices
    diag_n = dx^2/(dt*Dn);
    diag_mu = dx^2/(dt*Dmu);
    diag_G = dx^2/(dt*DG);
    v1 = ones(xdim-1,1);
    v0 = (diag_n+2)*ones(xdim,1);
    M_n = -diag(v1, -1) + diag(v0,0) -diag(v1,1);
    M_n(1,2) = -2;
    M_n(xdim, xdim-1) = -2;
    
    v0 = (diag_mu+2)*ones(xdim,1);
    M_mu = -diag(v1, -1) + diag(v0,0) -diag(v1,1);
    M_mu(1,2) = -2;
    M_mu(xdim, xdim-1) = -2;
    
    v0 = (diag_G+2)*ones(xdim,1);
    M_G = -diag(v1, -1) + diag(v0,0) -diag(v1,1);
    M_G(1,2) = -2;
    M_G(xdim, xdim-1) = -2;
end

if implicit
    % method algorithm
    v1 = zeros(xdim, 1);
    v2 = zeros(xdim, 1);
    v3 = zeros(xdim, 1);
    for i=2:tdim
        for j=1:xdim
            v1(j) = dx^2 / Dn * (nxt(j,i-1)/dt + f1(nxt(j,i-1), muxt(j,i-1), Gxt(j,i-1)));
            v2(j) = dx^2 / Dmu * (muxt(j,i-1)/dt + f2(nxt(j,i-1), muxt(j,i-1), Gxt(j, i-1)));
            v3(j) = dx^2 / DG * (Gxt(j,i-1)/dt + f3(nxt(j,i-1), Gxt(j,i-1)));
        end
        nxt(:,i) = linsolve(M_n, v1);
        muxt(:,i) = linsolve(M_mu, v2);
        Gxt(:,i) = linsolve(M_G, v3);
    end

end


%% figure
figure
pos = get(gcf, 'Position');  % [x0, y0, width, height]

pos(2) = 40;
pos(3) = 550;   
pos(4) = 750;   % così la figura è più lunga

set(gcf, 'Position', pos);

% n(x,t)
ax1 = subplot(3,1,1);   
h_1 = plot(ax1, xall, nxt(:,1), 'LineWidth', 1.5);  
xlabel(ax1, 'space')
ylabel(ax1, 'n(x,t)')
ylim(ax1, [0, max(nxt(:))])
title(ax1, 'n(x,t)')
grid(ax1, 'on')

% mu(x,t)
ax2 = subplot(3,1,2);
h_2 = plot(ax2, xall, muxt(:,1), 'LineWidth', 1.5);
xlabel(ax2, 'space')
ylabel(ax2, '\mu(x,t)')
ylim(ax2, [0, 1])     
title(ax2, '\mu(x,t)')
grid(ax2, 'on')

% G(x,t)
ax3 = subplot(3,1,3);
h_3 = plot(ax3, xall, Gxt(:,1), 'LineWidth', 1.5);
xlabel(ax3, 'space')
ylabel(ax3, 'G(x,t)')
ylim(ax3, [0, 1.2])    
title(ax3, 'G(x,t)')
grid(ax3, 'on')

if esplicit
    step=500; % to satisfy CFL dt must be very small
else
    step=100;
end

for i = 1:tdim
    if mod(i,step) == 0
        
        t_current = time(i);
        if isvalid(h_1)
            set(h_1, 'YData', nxt(:,i));
            title(ax1, sprintf('n(x,t) at t = %.2f', t_current))
        end
        if isvalid(h_2)
            set(h_2, 'YData', muxt(:,i));
            title(ax2, sprintf('\\mu(x,t) at t = %.2f', t_current))
        end
        if isvalid(h_3)
            set(h_3, 'YData', Gxt(:,i));
            title(ax3, sprintf('G(x,t) at t = %.2f', t_current))
        end

        drawnow limitrate
    end
end
hold off


%% integral over space of n
n_vec=zeros(1,tdim);
t_vec = [];
if implicit
    jump = 1;
else
    jump=100;
    % to speed up, for the explicit schemes with small dt
end

for i = 1:tdim
    n_vec(i) = sum(nxt(:,i));
    t_vec = [t_vec, time(i)];
end
n_vec = n_vec(n_vec ~= 0);

figure()
hold on
title('n(t) total over space', 'FontSize',18)
ax = gca;
ax.FontSize = 16; 
xlabel('Time[h]', 'FontSize',16)
ylabel('n(t)', 'FontSize',16)
plot(t_vec,n_vec, 'LineWidth',2.5)
hold off


