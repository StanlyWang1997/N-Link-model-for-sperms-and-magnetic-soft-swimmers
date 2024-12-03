% N+2 parameters solutions in a system of ODEs

clear all; close all;

NT = 40;  % time discretization

Nk = 5;   % number of beat cycles
L_tail = 100e-6; % total length of the tail [m]

N = 31; % the number of discretized links 

l_head = 10e-6; % the length of the head; the first segment is the head. [m]
w_head = 10e-6; % the minor length of the head [m]
t_head = 1e-7; % the thickness of the head [m]

l_tail = L_tail/(N-1); % the length of discrete segment [m] 
w_tail = 1e-6; % the width of discrete segment [m]
t_tail = 1e-7; % the thickness of discrete segment (along z direction) [m]


E = 1e8; % Young's modulus N.m^(-2) 1e8
J = w_tail*t_tail^3/12; % the moment of intertial [m^4]

%kappa = E*J/l_tail; % the spring constant between two adjacent segments

kappa = [0, ones(N-1,1)'*E*J/l_tail];

zeta_head = 0.05; % drag coefficient of the head [N.s.m^(-2)]
zetan_tail = 12.4e-3; % normal drag coefficient along the discrete segment [N.s.m^(-2)]
zetat_tail = 6.2e-3; % tangent drag coefficient along the discrete segment [N.s.m^(-2)]
gamma = zetan_tail/zetat_tail;

% Initial value (N+2) parameters, (x1, y1, theta1, alpha2, ... , alphaN)
% (x1, y1)-the coordinate of the head front; theta1-orientation of the head
% (1) Straight line
z0=[-l_head/2;0;0;zeros(N-1,1)];

% (2) Half-circle
% z0=[0;0;-pi/2;pi/(N-1)*ones(N-1,1)];

f = 1; % actuation frequency [Hz]
omega = 2*pi*f;   % angular frequency
T = 1/f;  % time for one beat cycle


%% case 1: 
% each segment, excluding only the first one describing the head, is constantly magnetized, and we make the simplifying assumption that the magnetization on each segment stays permanently aligned with the segment axis.
% Paper: Can Magnetic Multilayers Propel Artificial Microswimmers Mimicking Sperm Cells?

M_s = 8e5; % saturation magnetization [A.m^(-1)]
Mh_s = 8e7;
M_unitlen = M_s*t_tail*w_tail; % the magnetization per unit length of each segment [A.m]
Mi = M_unitlen*l_tail; % the (total) magnetization of the ith segment [A.m^2]

Mh = Mh_s*t_head*w_head*l_head; % the (total) magnetization of the head [A.m^2]

Mag = [0 ones(1,N-1)*Mi];  % tail magnetized, excluding head segment
% Mag = [Mh zeros(1, N-1)];  % head magnetized
% Mag = [Mh ones(N-1,1)'*Mi];  % head and tail magnetized

l = [l_head, l_tail*ones(1,N-1)];
kexi = [zeta_head, zetan_tail*ones(1,N-1)];
eta = [zeta_head, zetat_tail*ones(1,N-1)];

Sp = L_tail*(omega*zetan_tail/(E*J))^(1/4);  % sperm number characterizes the ratio between viscous force and elastic force.

%% Oscillating magnetic field
Bx = @(t) 0.01; % the magnetic field strength along x [T];
By = @(t) 0.01*sin(omega*t); % the magnetic field strength along y [T];

%% Circle motion
% Tmax = 20; % [s]
% theta = @(t) 2*pi*t/Tmax;
% Bx = @(t) 0.01*sin(theta(t))+0.01*sin(omega*t)*(cos(theta(t)));
% By = @(t) 0.01*cos(theta(t))+0.01*sin(omega*t)*(-sin(theta(t)));

%% Turning abruptly; thete(t) experiences a sudden jump from 0 to pi/2 aroud t = Tmax/2
% Tmax = 8; % [s]
% theta = @(t) pi/4*(1+tanh(30*(t/Tmax-1/2)));
% By = @(t) 0.01*cos(theta(t))+0.01*sin(omega*t)*(-sin(theta(t)));
% Bx = @(t) 0.01*sin(theta(t))+0.01*sin(omega*t)*(cos(theta(t)));

%% draw theta versus t
% fplot(theta, [0, Tmax]);
% xlabel('Time t');
% ylabel('\theta(t)');
% title('Plot of \theta(t) vs. Time t');

%% case: magnetic sperm microswimmer
% dZ = @(t, z) magnetism(t,z,N,Bx,By,l,kexi,eta,kappa,Mag);

%% case: sperm oscillation is for the standard case with a pinned end and angular actuation
% % ---- Choose an angular amplitude in rad
amp=0.4362;
% kappa_s=0;
kappa_s = 0.06*(E*J/l_tail)/l_tail; % an effective resistance to sliding between the sliding filaments 4.5e-11
dZ=@(t,z) oscillation(t,f,z,N,l,kexi,eta,amp,kappa,kappa_s);

tps = linspace(0,Nk*T,Nk*NT);

[tps,traj]=ode15s(dZ,tps,z0);  % calculate NT*NK times

% ---- Graphic visualisation 
% figure;
for i = 1:length(tps)
     [X,Y,TH(:,i)]=coordinates_filament(traj(i,:),N,l);
     centerX(i) = (X(1)+X(2))/2;
     centerY(i) = (Y(1)+Y(2))/2;

     % Calculate the semi-major and semi-minor axes
     a = norm([X(2) - X(1), Y(2) - Y(1)]) / 2;
     b = w_head/2; % Change this value as needed
     % Create a rotation angle if needed (set to 0 for no rotation)
     angle = atan((Y(2) - Y(1))/(X(2) - X(1)));
end

figure (5);
contourf(linspace(0,Nk*2,Nk*NT), linspace(0,(L_tail+l_head)*1e6,N), TH);
colormap jet;
% colorbar;
h=colorbar(gca,'Limits',[-0.8 0.8]);
set(get(h,'Title'),'string','\theta (rad)');
ylabel("\it L (um)",'FontSize',10);
xlabel("\it t (\pi)",'FontSize',10);

figure(2);
% TH: N*(Nk * NT) matrix
[X,Y,TH(:,1)]=coordinates_filament(traj(1,:),N,l);  % calculate X, Y, and theta for the initial solution at the first beat cycle
[X_end,Y_end,TH(:,end)]=coordinates_filament(traj(end,:),N,l);  % calculate X, Y, and theta for the last time at the last beat cycle
% Calculate the semi-major and semi-minor axes
a = l_head / 2;
b = w_head / 2; % Change this value as needed
% Create a rotation angle if needed (set to 0 for no rotation)
angle = atan((Y(2) - Y(1))/(X(2) - X(1)));
angle_end = atan((Y_end(2) - Y_end(1))/(X_end(2) - X_end(1)));
% Create an ellipse using the ellipse function
ellipse(a*1e6, b*1e6, angle, centerX(1)*1e6, centerY(1)*1e6); hold on;
plot(X(2:end)*1e6,Y(2:end)*1e6,'k','LineWidth',0.5)
plot(centerX*1e6, centerY*1e6);
plot(centerX(1)*1e6, centerY(1)*1e6,"*");

ellipse(a*1e6, b*1e6, angle_end, centerX(end)*1e6, centerY(end)*1e6);
plot(X_end(2:end)*1e6,Y_end(2:end)*1e6,'k','LineWidth',0.5)
plot(centerX(end)*1e6, centerY(end)*1e6,"o");

% xlim([-10 250]);
axis equal;

figure(3);
Np=5;  % number of tail deformation in one beat cycle
k=1;
for i = 1:round(NT/Np):NT
     k=k+1;
     [X_T,Y_T,TH_T(:,i)]=coordinates_filament(traj(i,:),N,l);
        
     angle_T = atan((Y_T(2) - Y_T(1))/(X_T(2) - X_T(1))); CenterX_T = (X_T(2) + X_T(1))/2; CenterY_T = (Y_T(2) + Y_T(1))/2;
     ellipse(a*1e6, b*1e6, angle_T, CenterX_T*1e6, CenterY_T*1e6, [1-k/(Np+3),1-k/(Np+3),1-k/(Np+3)],[],5);hold on;
     plot(X_T(2:end)*1e6,Y_T(2:end)*1e6,'k','LineWidth',5,'color',[1-k/(Np+3),1-k/(Np+3),1-k/(Np+3)]); 
%      ylim([-0.01 0.01]);
     xlabel("\it s",'FontSize',10);
     title(['T = ',num2str(tps(i)) ' s']);
     axis equal;
     drawnow
end

figure(4);
for i = 1:Nk
     [X,Y,TH(:,i)]=coordinates_filament(traj((i)*NT,:),N,l);
     center_X = (X(1)+X(2))/2;
     center_Y = (Y(1)+Y(2))/2;

     % Calculate the semi-major and semi-minor axes
     a = l_head / 2;
     b = w_head/2; % Change this value as needed
     % Create a rotation angle if needed (set to 0 for no rotation)
     angle = atan((Y(2) - Y(1))/(X(2) - X(1)));

     % Create an ellipse using the ellipse function
     ellipse(a*1e6, b*1e6, angle, center_X*1e6, center_Y*1e6); hold on;
     plot(X(2:end)*1e6,Y(2:end)*1e6,'k','LineWidth',0.5)
%      axis([-1 1 -1 1])
%      xlim([-10 250]);
     axis equal;
     title(['T = ',num2str(tps(i))]);
     drawnow
end

%% define the system of ODEs that determines the motion of the magnetoelastic swimmer, t and z are input parameters
function B=magnetism(t,z,N,Bx,By,l,kexi,eta,kappa,Mag)

th=zeros(N,1);
th(1) = z(3);
for i =2:N
    th(i) = th(i-1)+z(i+2);
end
F0 = zeros(N+2,1);
F1 = zeros(N+2,1);
F2 = zeros(N+2,1);
F1(N+2) = Mag(N)*sin(th(N));
F2(N+2) = -Mag(N)*cos(th(N));

% adding the magnetic torque
for i = N-2:-1:0
    F0(4+i) = kappa(i+2)*(th(i+2)-th(i+1)); %correct
    F1(3+i) = F1(3+i+1)+Mag(i+1)*sin(th(i+1));
    F2(3+i) = F2(3+i+1)-Mag(i+1)*cos(th(i+1));
end
FF = F0+F1*Bx(t)+F2*By(t);
M=matrixNparam(t,z,N,l,kexi,eta);

% % solve pseudoinverse using svd decomposition
% [U,S,V] = svd(M); % SVD decomposition of M 
% T=S;
% T(find(S~=0)) = 1./S(find(S~=0));  %reciprocal of nonzero elements in M
% svdInvM = V * T' * U';
% B = svdInvM*FF;

% % solve pseudoinverse using QR decomposition
% [Q,R] = qr(M);
% InvR =  inv(R'*R)*R';
% qrInvM =InvR*Q';
% B = qrInvM*FF;

B = M\FF;
% B = pinv(M)*FF; %Moore-Penrose pseudoinverse

end

% **************************************
function B=oscillation(t,f,z,N,l,kexi,eta,amp,kappa,kappa_s)
% This function solves the linear system of sperm flagellum
% from Appendix VII-C in the case of a pinned end with a forced
% angular actuation
% z - (N+2)X1 vector; x1_dot, y1_dot, theta1_dot, alpha2_dot, ...
% alphaN_dot
th=zeros(N,1);
th(1)=z(3);
for i=2:N
    th(i)=th(i-1)+z(i+2);
end

F0=zeros(N+2,1);
F1=zeros(N+2,1);
F1(N+2)=kappa_s*l(N)*(th(N)-th(1)); %% revisit

% for forced angular actuation
a0p=(amp*2*pi*f)*cos(2*pi*f*t);

for i=N-2:-1:0
    F0(4+i) = kappa(i+2)*(th(i+2)-th(i+1)); % elastic spring
    F1(3+i)= F1(4+i)+kappa_s*l(i+1)*(th(1+i)-th(1)); %% revisit
end

FF = F0+F1;
FF(3)=a0p;

M=matrixNparam_oscillation(t,z,N,l,kexi,eta);

% % solve pseudoinverse using svd decomposition
% [U,S,V] = svd(M); % SVD decomposition of M 
% T=S;
% T(find(S~=0)) = 1./S(find(S~=0));  %reciprocal of nonzero elements in M
% svdInvM = V * T' * U';
% B = svdInvM*FF;

% % solve pseudoinverse using QR decomposition
% [Q,R] = qr(M);
% InvR =  inv(R'*R)*R';
% qrInvM =InvR*Q';
% B = qrInvM*FF;

% B=M\FF;
B = mldivide(M,FF);
% B = lsqminnorm(M,FF);

end

function  M=matrixNparam_oscillation(t,z,N,l,kexi,eta)
% This function is similar to matrixNparam but in the 
% case of a pinned proximal end with a forced angular actuation for sperm cells.
% z3 x1_dot,...,xN_dot; y1_dot, ..., yN_dot; 

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);   % th2...thN
    z3(i)=z3(i-1)+cos(z3(2*N+i-1))*l(i-1);  % x2....xN
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1))*l(i-1);  % y2 ... yN
end

M3=matrix3Nparameters(z3,N,l,kexi,eta);

C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;  % >= has higher preference than '=' .
    end
end

for i=2:N
    C1(i,i-1)=-sin(z3(2*N+i-1))*l(i-1);
    C2(i,i-1)=cos(z3(2*N+i-1))*l(i-1);
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-sin(z3(2*N+j))*l(j);
        C2(i,j)=C2(i,j+1)+cos(z3(2*N+j))*l(j);
    end
end
C=[ones(N,1),zeros(N,1),C1;zeros(N,1),ones(N,1),C2;zeros(N,2),C3];   % calculate Q
M=M3*C;

% --- this is where the first two equations are being changed to 
% implement the pinned proximal end (boundary conditions)

M(1,:)=[1,zeros(1,N+1)];  % xdot = 0
M(2,:)=[0,1,zeros(1,N)];  % ydot = 0
M(3,:)=[0,0,1,zeros(1,N-1)]; % thetadot = amp*cos(t)
end


function M=matrixNparam(t,z,N,l,kexi,eta)
% This function fills the matrix Q as defined in the text, covert 3N
% parameters to (N+2) parameters for magentic sperm microswimmers
% and returns the product A*Q

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);   % th2...thN
    z3(i)=z3(i-1)+cos(z3(2*N+i-1))*l(i-1);  % x2....xN
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1))*l(i-1);  % y2 ... yN
end

M3=matrix3Nparameters(z3,N,l,kexi,eta);

C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;  % >= has higher preference than '=' .
    end
end
for i=2:N
    C1(i,i-1)=-sin(z3(2*N+i-1))*l(i-1);
    C2(i,i-1)=cos(z3(2*N+i-1))*l(i-1);
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-sin(z3(2*N+j))*l(j);
        C2(i,j)=C2(i,j+1)+cos(z3(2*N+j))*l(j);
    end
end
C=[ones(N,1),zeros(N,1),C1;zeros(N,1),ones(N,1),C2;zeros(N,2),C3];   % calculate Q
M=M3*C;
end   % C1, C2, and C3 are correct

function M=matrix3Nparameters(z,N,l,kexi,eta)
% This function fills the matrix defined as A in the text 
% and returns A （N+2）X 3N

x=z(1:N);
y=z(N+1:2*N);
th=z(2*N+1:3*N);
F=zeros(2,3*N);
T=zeros(N,3*N);
Tx=zeros(N,N);
Ty=zeros(N,N);
Tth=zeros(N,N);
Fxx = zeros(1,N);
Fxy = zeros(1,N);
Fyx = zeros(1,N);
Fyy = zeros(1,N);
Fxth = zeros(1,N);
Fyth = zeros(1,N);

for i = 1:N
    cos_th = cos(th(i));
    sin_th = sin(th(i));
    cos2_th = cos_th^2;
    sin2_th = sin_th^2;

    Fxx(i) = -l(i) * (kexi(i) * cos2_th + eta(i) * sin2_th);
    Fyx(i) = -l(i) * (kexi(i) - eta(i)) * sin_th * cos_th;
    Fxy(i) = Fyx(i);
    Fyy(i) = -l(i) * (kexi(i) * sin2_th + eta(i) * cos2_th);
    Fxth(i) = l(i)^2 / 2 * eta(i) * sin_th;
    Fyth(i) = -l(i)^2 / 2 * eta(i) * cos_th;

end

for i = 1:N
%     Fxx(i) = -l(i)*(kexi(i)*cos(th(i))^2+eta(i)*sin(th(i))^2);
%     Fyx(i) = -l(i)*(kexi(i)-eta(i))*sin(th(i))*cos(th(i));
%     Fxy(i) = Fyx(i);
%     Fyy(i) = -l(i)*(kexi(i)*sin(th(i))^2+eta(i)*cos(th(i))^2);
%     Fxth(i) = l(i)^2/2*eta(i)*sin(th(i));
%     Fyth(i) = -l(i)^2/2*eta(i)*cos(th(i)); % F correct
    for j = i:N
        dx = x(j) - x(i);
        dy = y(j) - y(i);
        cos_th_j = cos(th(j));
        sin_th_j = sin(th(j));
        cos2_th_j = cos_th_j^2;
        sin2_th_j = sin_th_j^2;

%         Tx(i,j) = eta(j)*l(j)^2/2*sin(th(j))-(x(j)-x(i))*l(j)*(kexi(j)-eta(j))*sin(th(j))*cos(th(j))+(y(j)-y(i))*l(j)*(kexi(j)*cos(th(j))^2+eta(j)*sin(th(j))^2);
%         Ty(i,j) = -eta(j)*l(j)^2/2*cos(th(j))-(x(j)-x(i))*l(j)*(eta(j)*cos(th(j))^2+kexi(j)*sin(th(j))^2)+(y(j)-y(i))*l(j)*(kexi(j)-eta(j))*sin(th(j))*cos(th(j));
%         Tth(i,j) = -eta(j)*l(j)^3/3-(x(j)-x(i))*l(j)^2/2*eta(j)*cos(th(j))-(y(j)-y(i))*l(j)^2/2*eta(j)*sin(th(j));
        Tx(i, j) = eta(j) * l(j)^2 / 2 * sin_th_j - dx * l(j) * (kexi(j) - eta(j)) * sin_th_j * cos_th_j + dy * l(j) * (kexi(j) * cos2_th_j + eta(j) * sin2_th_j);
        Ty(i, j) = -eta(j) * l(j)^2 / 2 * cos_th_j - dx * l(j) * (eta(j) * cos2_th_j + kexi(j) * sin2_th_j) + dy * l(j) * (kexi(j) - eta(j)) * sin_th_j * cos_th_j;
        Tth(i, j) = -eta(j) * l(j)^3 / 3 - dx * l(j)^2 / 2 * eta(j) * cos_th_j - dy * l(j)^2 / 2 * eta(j) * sin_th_j;  % correct
    end
end

F = [Fxx Fxy Fxth; Fyx Fyy Fyth];
T = [Tx Ty Tth];

M=[F;T];
end

function [X,Y,TH]=coordinates_filament(z,N,l)
% This function computes the '3N coordinates' -- X_3N in the text
% from the 'N+2 coordinates' -- X in the text 

% --- input : N+2 coordinates, number of links N
% --- output : X, Y coordinates of the N links
% TH orientation of each link

X=zeros(N+1,1);
Y=zeros(N+1,1);
TH=zeros(N,1);

X(1)=z(1);
Y(1)=z(2);
TH(1)=z(3);                              

for i=2:N
    X(i)=X(i-1)+cos(TH(i-1))*l(i-1);
    Y(i)=Y(i-1)+sin(TH(i-1))*l(i-1);
    TH(i)=TH(i-1)+z(i+2);
end

X(N+1)=X(N)+cos(TH(N))*l(N);
Y(N+1)=Y(N)+sin(TH(N))*l(N);

end

