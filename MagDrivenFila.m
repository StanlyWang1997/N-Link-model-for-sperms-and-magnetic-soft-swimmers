% *******************************************
% Main function: integrating the dynamics equation of an elastic microfilament
% *******************************************
%  The results in the reference titled "Elastohydrodynamic propulsion of a filament magnetically driven at both ends" are reproduced throughout the code
%
%  Returns N+2 parameters for each time (contained in the 'traj' array) : 
%   * x and y are the coordinates of the end of the first link
%   * theta is the orientation of the first link
%   * alpha2 to alphaN are the 'shape angles' : angle between i+1-th and i-th link
% *******************************************
%
clear all; close all;
tic; % start timing

% ---- Choose number of links
N=50;

% Declaring global parameters
global gamma Sp

% ---- Choose 'sperm number' characterizing the 'floppiness' of the filament - Sp=L*(zeta*omega/kappa)^(1/4)
% typically between 2 and 16, see the section 2.3
Sp=3; % >=4.5, the calculation is very slow
gamma=2; % ratio between the hydrodynamics drag coefficients - gamma = xi/eta, normal drag coefficient / tangent drag coefficient

% ---- Choose initial condition

% z0=zeros(N+2,1);
z0=[-1/2; zeros(N+1,1)]; % the centroid is (0,0)

% ---- Choose final time and time step
NT = 10;
Nt = 50; % 50-100
T=2*pi; %time for one period
tps=linspace(0,NT*T,NT*Nt); %time step
t_end = tps(end);

% magnetism is for the magnetised case
% % ---- Choose a magnetisation for the filament (N values)
M0 = 0; % magnetic moment strength at the left end
ML = 1-M0;
Mag = Sp^4*[-M0, zeros(1,N-2), ML]; % 
% % ---- Choose time-varying magnetic fields Hx and Hy (see section 2.2)
% % (I recommend not to use fields with absolute value greater than 1)
Hx=@(t) 1;
Hy=@(t) 1*sin(t);
dZ=@(t,z) magnetism(t,z,N,Mag,Hx,Hy,t_end);
% ************

% ---- SOLVING
[tps,traj]=ode15s(dZ,tps,z0);
centroid_x=zeros(length(tps),1);
centroid_y=zeros(length(tps),1);

% ---- Graphic visualisation 
figure(1);
for i = 1:length(tps)
     [X,Y,TH(:,i)]=coordinates_filament(traj(i,:),N);
     centroid_x(i) = X(N/2+1);
     centroid_y(i) = Y(N/2+1);

     plot(X,Y,'k','LineWidth',0.5); hold on;

      % plot centroid
     plot(centroid_x(1:i), centroid_y(1:i), 'r-', 'LineWidth', 0.5);
%      axis([-1 1 -0.5 0.5]);
     xlim([-1 1]);
     ylim([-0.5 0.5]);
%      axis equal;
     title(['T = ',num2str(tps(i))]);
     hold off;
     drawnow
end

D = zeros(1,NT-1);
for i=1:NT-1
    D(i) = centroid_x(i*Nt+1)-centroid_x((i-1)*Nt+1);
end
figure(4);
plot(1:(NT-1), D,'ro');
ylim([0 0.12]);
xlabel('Number of cycle');
ylabel('\itD');
title('Dimensionless displacement of the filament in a single cycle');

figure(2);
Np=6;
k=1;
for i = round(length(tps)/NT/Np+(NT-1)*length(tps)/NT):round(length(tps)/NT/Np):length(tps)
     k=k+1;
     [X,Y,TH(:,i)]=coordinates_filament(traj(i,:),N);
     plot(X,Y,'k','LineWidth',2,'color',[1-k/(Np+3),1-k/(Np+3),1-k/(Np+3)]); hold on;
     ylim([-0.2 0.2]);
     axis equal;
     xlabel("\it s",'FontSize',10);
     title(['T = ',num2str(tps(i))]);
end
plot(centroid_x(1:i), centroid_y(1:i), 'r-', 'LineWidth', 1); hold off;

for i = 1:length(tps)
    [X,Y,TH(:,i)]=coordinates_filament(traj(i,:),N);
end

figure (3);
contourf(linspace(0,NT,Nt*NT), linspace(0,1,N), TH);
colormap jet;
colorbar;
colorbar(gca,'Limits',[-0.44 0.44]);
ylabel("\its",'FontSize',10);
xlabel("\itt (2*pi)",'FontSize',10);

elapsedTime = toc; % stop timing
hours = floor(elapsedTime / 3600);
minutes = floor(mod(elapsedTime, 3600) / 60);
seconds = round(mod(elapsedTime, 60));

disp(['Total execution time: ' num2str(hours,'%02d'), 'h ', num2str(minutes,'%02d'), 'm ', num2str(seconds,'%02d'), 's.']);

% *********** END OF MAIN SCRIPT ***********
% ***** SUMMARY OF FOLLOWING FUNCTIONS *****
% magnetism
% matrixNparam
% matrix3Nparameters
% coordinates_filament

% **************************************
function B=magnetism(t,z,N,Mag,Hx,Hy, t_end)
% This function is similar to second_member, but with the magnetic effects added
% (see Appendix VII-E, Eq. 24)
% Hx and Hy are the external magnetic fields, defined in the MAIN file.
th=zeros(N+1,1);
th(2)=z(3);
for i=3:N+1
    th(i)=th(i-1)+z(i+1);
end
B1=zeros(N+2,1);
B2=zeros(N+2,1);
B3=zeros(N+2,1);
B1(N+2)=N*(th(N+1)-th(N));

% adding the magnetic effects 
B2(N+2)=Mag(N)*sin(th(N));
B3(N+2)=-Mag(N)*cos(th(N));
for i=N-2:-1:0
    B1(3+i)=N*(th(i+2)-th(i+1));
    B2(3+i)=B2(3+i+1)+Mag(i+1)*sin(th(i+2));
    B3(3+i)=B3(3+i+1)-Mag(i+1)*cos(th(i+2));
end
B1(3)=0;

BB=(B1+Hx(t)*B2+Hy(t)*B3);
M=matrixNparam(t,z,N);
%B = mldivide(M,BB);
% B = lsqminnorm(M,BB);
% B=M\BB;
M = sparse(M);

if issymmetric(M) && isposdef(M)
%     R = chol(M); 
%     B = R' \ (R \ BB);
    [L,U,P] = lu(M);
    B = U \ (L \ (P * BB));
else
    B = M \ BB;
end

% [L,U,P] = lu(M);
% B = U \ (L \ (P * BB));

% fprintf('\rThe calculation is finished in %.2f%%.', t/t_end*100);
disp(['The calculation is finished in ', num2str(t/t_end*100,'%.3f'), '%.']);

clear th B1 B2 B3 M BB

end

% **************************************
function M=matrixNparam(t,z,N)
% This function fills the matrix Q as defined in the text
% and returns the product Sp^4*A*Q

z3=zeros(3*N,1);
z3(1)=z(1);
z3(N+1)=z(2);
z3(2*N+1)=z(3);
for i=2:N
    z3(2*N+i)=z3(2*N+i-1)+z(i+2);
    z3(i)=z3(i-1)+cos(z3(2*N+i-1))/N; %correct
    z3(N+i)=z3(N+i-1)+sin(z3(2*N+i-1))/N; % correct
end

M3=matrix3Nparameters(t,z3,N);

C1=zeros(N,N);C2=zeros(N,N);C3=zeros(N,N);
for i=1:N
    for j=1:N
        C3(i,j)=i>=j;  % >= has higher preference than '=' .
    end
end
for i=2:N
    C1(i,i-1)=-sin(z3(2*N+i-1));
    C2(i,i-1)=cos(z3(2*N+i-1));
end
for i=N:-1:3
    for j=i-2:-1:1
        C1(i,j)=C1(i,j+1)-sin(z3(2*N+j));
        C2(i,j)=C2(i,j+1)+cos(z3(2*N+j));
    end
end
C=[ones(N,1),zeros(N,1),C1/N;zeros(N,1),ones(N,1),C2/N;zeros(N,2),C3];
M=M3*C;

clear z3 M3 C1 C2 C3 C

end

% **************************************
function M=matrix3Nparameters(t,z,N)
% This function fills the matrix defined as A in the text 
% and returns Sp^4 * A

global gamma Sp

x=z(1:N);
y=z(N+1:2*N);
th=z(2*N+1:3*N);
F=zeros(2,3*N);
T=zeros(N,3*N);

for i=1:N
    u=cos(th(i));
    v=sin(th(i));
    F(1,i)=-(u^2+gamma*v^2);
    F(1,N+i)=(gamma-1)*u*v;
    F(2,i)=(gamma-1)*u*v;
    F(2,N+i)=-(gamma*u^2+v^2);
    F(1,2*N+i)=gamma/2*v/N;
    F(2,2*N+i)=-gamma/2*u/N;
end

F=Sp^4*F/N/gamma; % here divides by N, dividing by gamma：Sp^4/gamma = L^4*(kesi*omega)/E

for i=1:N
    for j=i:N
        u=cos(th(j));
        v=sin(th(j));
        A=(x(j)-x(i));
        B=(y(j)-y(i));
        T(i,j)=gamma/2*v/N+...
            A*(gamma-1)*v*u+...
            B*(gamma*v*v+u*u);
        T(i,N+j)=-gamma/2*u/N+...
            B*(1-gamma)*v*u-...
            A*(gamma*u*u+v*v);
        T(i,2*N+j)=-gamma/3/N/N-...
            gamma/2*A*u/N...
            -gamma/2*B*v/N;
    end
end

T=Sp^4*T/N/gamma; % here divides by N

M=[F;T];

clear x y th F T u v A B
end

% **************************************
function [X,Y,TH]=coordinates_filament(z,N)
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
    X(i)=X(i-1)+cos(TH(i-1))/N;
    Y(i)=Y(i-1)+sin(TH(i-1))/N;
    TH(i)=TH(i-1)+z(i+2);
end

X(N+1)=X(N)+cos(TH(N))/N;
Y(N+1)=Y(N)+sin(TH(N))/N;

end
