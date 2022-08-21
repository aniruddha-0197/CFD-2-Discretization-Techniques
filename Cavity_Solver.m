%% Navier-Stokes CFD Solver for a Lid-Driven Cavity Flow
% (c) Aaron Sequeira, Aniruddha Panajape, 2022

clc 
clearvars
close all
set(0,'defaultTextInterpreter','latex');
warning('off','MATLAB:nearlySingularMatrix')

%% Settings

% General
N   = 31;	    % Number of Divisions
Re  = 1000;	    % Reynolds Number
tol = 1e-05;    % Covergence Criterion

mdisp = 'on';      % Show mesh?

% Boundary Conditions
u_L = 0;		        % u - Left Wall
u_R = 0;		        % u - Right Wall
u_B = 0;		        % u - Bottom Wall
u_T = -1;	            % u - Top Wall
v_L = 0;		        % v - Left Wall
v_R = 0;		        % v - Right Wall
v_B = 0;		        % v - Bottom Wall
v_T = 0;		        % v - Top Wall

%% Generate Mesh

% Primal Points and Spacing
x_p = 0.5*(1 - cos(pi*linspace(0,1,N+1)));
h_p = x_p(2:N+1) - x_p(1:N);

% Dual Points and Spacing
x_d = (x_p(2:N+1) + x_p(1:N))/2;
x_d = [0 x_d 1];
h_d = x_d(2:N+2) - x_d(1:N+1);

[Xp,Yp] = meshgrid(x_p);            % 2D Primal mesh points
[Xd,Yd] = meshgrid(x_d(2:end-1));   % 3D Dual mesh points

% Plot Mesh
if strcmpi(mdisp,'on')
    figure('Name','Mesh')
    plot(Xp,Yp,'k',Xp',Yp','k')
    %plot(Xp,Yp,'k',Xd,Yd,'r',Xp',Yp','k',Xd',Yd','r')
    xlabel('$x$')
    ylabel('$y$')
    title('Primal Mesh')
    set(gca,'TickLength',[0 0])
    drawnow()
end

%% Discretization

% Compute Incidence Matrices
[tE21, E10, E_norm] = incidence_1(N);
[E21, tE10, E_tan]  = incidence_2(N);

% Compute Hodge Matrices
[Ht1_1, H1_t1, Ht0_2] = hodge(N,h_p,h_d);

% Apply Boundary Conditions
[V_norm, V_tan] = boundary(u_L,u_R,u_B,u_T,v_L,v_R,v_B,v_T,N,h_d);
F_div           = E_norm*V_norm.';
F_vort          = E_tan*V_tan.';

%% Solver

tic     % Start the clock

% Construct Pressure Matrix
A           = tE21*Ht1_1*E10;         % Pressure Laplacian with Neumann BC
[eV,lambda] = eig(full(A));           % Get Eigenspectrum of the Pressure Matrix

% Initialize Velocity and Pressure
u = zeros(2*N*(N+1),1);
p = zeros(N);

% Initialize Convective Terms
ux_xi = zeros((N+1)*(N+1),1);
uy_xi = zeros((N+1)*(N+1),1);
convective = zeros(2*N*(N+1),1);

% Time Stepping
h_min = min(h_d);    						% Minimum Dual Edge Length
dt 	  = 5*min([h_min,0.5*Re*h_min^2]);	    % Time Step Size
res   = 1;									% Initial Residual
itr   = 1;									% Iteration Counter
iuns  = 1;                                  % Unsteady Counter

while res >= tol

	% Calculate (Dual) Vorticity 
	xi = E21*u;

	% Diffusive Term
	diffusive = H1_t1*(tE10*(Ht0_2*(xi + F_vort)))/Re;

	% Convective Term
	xi_c = Ht0_2*(xi + F_vort);
    
    for i=1:N+1
        for j=1:N+1
            k = j + (i-1)*(N+1); 
            if j==1
                ux_xi(k) = u_B*xi_c(i+(j-1)*(N+1));    
                uy_xi(k) = v_L*xi_c(j+(i-1)*(N+1));
            elseif j==N+1
                ux_xi(k) = u_T*xi_c(i+(j-1)*(N+1));
                uy_xi(k) = v_R*xi_c(j+(i-1)*(N+1));
            else
                ux_xi(k) = (u(i+(j-1)*(N+1))+u(i+(j-2)*(N+1)))*xi_c(i+(j-1)*(N+1))/(2*h_d(i));                      
                uy_xi(k) = (u(N*(N+1)+j+(i-1)*N) + u(N*(N+1)+j-1+(i-1)*N))*xi_c(k)/(2*h_d(i));
            end
        end
    end

    for  i=1:N
        for j=1:N+1
            convective(j+(i-1)*(N+1))     = -(uy_xi(j+(i-1)*(N+1))+uy_xi(j+i*(N+1)))*h_d(j)/2;
            convective(N*(N+1)+i+(j-1)*N) = (ux_xi(j+(i-1)*(N+1))+ux_xi(j+i*(N+1)))*h_d(j)/2;
        end
    end

	% Assemble RHS function
	f = tE21*( Ht1_1*( u/dt - convective - diffusive ) ) + F_div/dt ;

	% Solve Poisson Eqn for (Total)Pressure
	P = A\f;

	% Update Velocity
	u0 = u;
	u  = u - dt*( convective + E10*P + diffusive );

	% Monitor Solution
	res 	= max(abs(u - u0)/dt);              % Time Derivative
	maxdiv	= max(abs(tE21*Ht1_1*u + F_div));   % Check for Mass Conservation

    if mod(itr,100)==0
        fprintf("Time step "+itr+"    Time: "+round(itr*dt,3)+"s \n")
        fprintf("Max Res: "+res+"\n")
        fprintf("Max Div: "+maxdiv+"\n\n")

        u_uns(:,iuns) = u;              %#ok
        t_uns(iuns)   = round(itr*dt);  %#ok
        iuns          = iuns+1;
    end
		
	itr = itr + 1;      % Rinse and Repeat

end

fprintf("Converged! ("+round(toc,2)+" s) \n")

%% Post Processing

% Velocity Field
[Vel,U,V] = vel_int(u,h_d,N);
velocity  = reshape(Vel,[N,N]).';
U         = reshape(U,[N,N]).';    U_xmid = [0,U(:,(N+1)/2).',-1]; U_ymid = [0,U((N+1)/2,:),0];
V         = reshape(V,[N,N]).';    V_xmid = [0,V(:,(N+1)/2).',0]; V_ymid  = [0,V((N+1)/2,:),0];

% Static/Total Pressure Field
P         = P - P((N^2+1)/2);
Tpressure = reshape(P(1:N^2),[N,N]).';
pressure  = Tpressure - 0.5*velocity.^2;

P_xmid = [P(N^2+2*N+(N+1)/2),Tpressure(:,(N+1)/2).',P(N^2+3*N+(N+1)/2)]; P_ymid = [P(N^2+N),Tpressure((N+1)/2,:),P(N^2+N+1)];
p_xmid = P_xmid - 0.5*(U_xmid.^2 + V_xmid.^2); p_ymid = P_ymid - 0.5*(U_ymid.^2 + V_ymid.^2);

% Vorticity Field
vorticity = reshape(xi_c,[N+1,N+1]).';
xi_xmid   = vorticity(:,(N+1)/2).'; xi_ymid   = vorticity((N+1)/2,:);
xi_int    = sum(xi + F_vort);

% Stream Function
A_psi = Ht0_2*E21*H1_t1*tE10;
f_psi = Ht0_2*E21*u;

psi    = A_psi(2:end,2:end)\f_psi(2:end);
stream = reshape([0;psi],[N+1,N+1]).';

% Contour Levels from Botella & Peyret
p_levels   = [0.3 0.17 0.12 0.11 0.09 0.07 0.05 0.02 0 -0.002];
xi_levels  = [5,4,3,2,1,0.5,0,-0.5,-1,-2,-3];
psi_levels = [0.1175 0.115 0.11 0.1 9e-2 7e-2 ...
              5e-2 3e-2 1e-2 1e-4 1e-5 1e-10 0 ...
              -1e-6 -1e-5 -5e-5 -1e-4 -2.5e-4 ...
              -5e-4 -1e-3 -1.5e-3];

% % Saving...
% save("Data\results_N"+N+".mat",'velocity','U','V','pressure','vorticity','stream',...
%     'U_xmid','V_xmid','p_xmid','xi_xmid','U_ymid','V_ymid','p_ymid','xi_ymid','xi_int','x_p','x_d','N');

%% Plotting

% Import Botella & Peyret Data
BP = readtable("BP.csv");

% Contours
figure('Name','Fields')
colormap parula

subplot(2,2,1)
contourf(Xd,Yd,velocity,20,'Linestyle','none')
xlabel('$x$')
ylabel('$y$')
title('Velocity Magnitude')
colorbar

subplot(2,2,2)
contourf(Xd,Yd,pressure,p_levels)
xlabel('$x$')
ylabel('$y$')
title('Pressure')
colorbar

subplot(2,2,3)
contourf(Xp,Yp,stream,psi_levels)
xlabel('$x$')
ylabel('$y$')
title('Stream Function')
colorbar

subplot(2,2,4)
contourf(Xp,Yp,vorticity,xi_levels);
xlabel('$x$')
ylabel('$y$')
title('Vorticity Magnitude')
colorbar

% Centerline Plots
figure('Name','Lines')

subplot(3,2,1)
plot(U_xmid,x_d,'r',BP.u,BP.y,'ko')
ylabel('$y$')
xlabel('$u$')
title("$u$-component of Velocity at $x$ = 0.5")
grid on

subplot(3,2,2)
plot(x_d,V_ymid,'r',BP.x,BP.v,'ko')
xlabel('$x$')
ylabel('$v$')
title("$v$-component of Velocity at $y$ = 0.5")
grid on

subplot(3,2,3)
plot(p_xmid,x_d,'r',BP.pv,BP.y,'ko')
ylabel('$y$')
xlabel('$p$')
title("Static Pressure at $x$ = 0.5")
grid on

subplot(3,2,4)
plot(x_d,p_ymid,'r',BP.x,BP.ph,'ko')
xlabel('$x$')
ylabel('$p$')
title("Static Pressure at $y$ = 0.5")
grid on

subplot(3,2,5)
plot(xi_xmid,x_p,'r',BP.xiv,BP.y,'ko')
ylabel('$y$')
xlabel('$\xi$')
title("Vorticity at $x$ = 0.5")
grid on

subplot(3,2,6)
plot(x_p,xi_ymid,'r',BP.x,BP.xih,'ko')
xlabel('$x$')
ylabel('$\xi$')
title("Vorticity at $y$ = 0.5")
grid on

%% Unsteady Animation

n_uns = size(u_uns);

figure('Name','Unsteady')
colormap jet

for ii = 1:n_uns(2)

    clf()

    [Vel_uns,U_uns,V_uns] = vel_int(u_uns(:,ii),h_d,N);
    velocity_uns  = reshape(Vel_uns,[N,N]).'; U_uns = reshape(U_uns,[N,N]).'; V_uns = reshape(V_uns,[N,N]).';

    contourf(Xd,Yd,velocity_uns,20,'Linestyle','none')
    hold on
    quiver(Xd,Yd,U_uns,V_uns,'k')
    xlabel('$x$')
    ylabel('$y$')
    title("Velocity Field at $t = $ "+t_uns(ii)+" s")
    colorbar
    xlim([0 1])
    ylim([0 1])
    hold off

    pause(0.1)

end
