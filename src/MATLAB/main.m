%% main.m
%
% Description:
%   Program entry point. Calls various functions.
%
% Dependencies:
%   Modern Robotics
%
% Inputs:
%   none
%
% Outputs:
%   none

function main
%% Initialization
clear;
close all;
clc;

init_env;
body_name = 'rect_prism';
dimensions.rect_prism.l = 1e-1;
dimensions.rect_prism.w = 1e-1;
dimensions.rect_prism.h = 1e-1;
params = init_params(body_name,dimensions);

%% Set initial conditions and numerical integration parameters
Rwb_init = eye(3);
Twb_init = [Rwb_init,   zeros(3,1);
            zeros(1,3), 1]; % CoM frame coincident with world frame
Vb_init = [1;2;3;1;2;3];
Fb = zeros(6,1);
        
dt = 1e-3; % timestep size
t0 = 0;
tf = 10; % simulate for 1 second
t = t0:dt:tf;

%% Simulate via numerical integration
Twb{1} = Twb_init;
Vb{1} = Vb_init;
for k = 1:numel(t)-1
    [~,Twb{k+1},Vb{k+1}] = cg4(@rbdyn,dt,t(k),Twb{k},Vb{k},Fb,params);
end

%% Postprocess
% kinematics:
Twb_hist = NaN(4,4,numel(t));
Vb_hist = NaN(6,numel(t));
p_hist = NaN(3,numel(t));
v_hist = NaN(3,numel(t));

% wrench in {b} and {c} frames:
Fc_hist = NaN(6,numel(t));
Fb_hist = NaN(6,numel(t));

% determinant of Rwb:
detRwb_hist = NaN(1,numel(t));
distToSE3_hist = NaN(1,numel(t));

% energy and momenta:
E_hist = NaN(1,numel(t)); % mechanical energy
Lb_hist = NaN(3,numel(t)); % angular momentum in body frame
Pb_hist = NaN(3,numel(t)); % linear momentum in body frame
Lw_hist = NaN(3,numel(t)); % angular momentum in world frame
Pw_hist = NaN(3,numel(t)); % linear momentum in world frame

for k = 1:numel(t)
    Rwb = Twb{k}(1:3,1:3); % orientation of {b} in {w}
    pwb = Twb{k}(1:3,4);   % position of origin of {b} in {s}
    
    % energy
    Gb = params.dyn.spatial_inertia_matrix_CoM;
    KE = 0.5*transpose(Vb{k})*Gb*Vb{k};
    PE = params.dyn.mass*params.dyn.grav*pwb(3);
    E = KE + PE;
    
    % momentum:
    Pb = Gb*Vb{k};
    Lb_hist(:,k) = [eye(3),zeros(3)]*Pb; % angular momentum in {b}
    Pb_hist(:,k) = [zeros(3),eye(3)]*Pb; % linear momentum in {b}
    Tbw = TransInv(Twb{k});
    Pw = transpose(Adjoint(Tbw))*Pb;
    Lw_hist(:,k) = [eye(3),zeros(3)]*Pw; % angular momentum in {w}
    Pw_hist(:,k) = [zeros(3),eye(3)]*Pw; % linear momentum in {w}
    
    % velocity of CoM in world frame:
    v = Rwb*Vb{k}(4:6);
    
    Twb_hist(:,:,k) = Twb{k};
    Vb_hist(:,k) = Vb{k};
    p_hist(:,k) = pwb;
    v_hist(:,k) = v;
    detRwb_hist(k) = det(Rwb);
    distToSE3_hist(k) = DistanceToSE3(Twb{k});
    E_hist(k) = E;
end

%% Visualize simulation results
% handy diagnostic functions:
%   total energy
%   angular momenta
%   linear momenta in x and y (but not z, b/c of gravity)
%   determinant of Rwb (how close is it to 1?)
%   distance to SE(3) manifold (based on Frobenius norm of Twb)

figure;
plot(t,detRwb_hist,'k:',t,distToSE3_hist,'k-','LineWidth',2)
legend('$\mathrm{det}\left(R_\mathrm{wb}\right)$',...
    'distance to $SE(3)$','Location','Best')
xlabel('time $t$ [s]')
axis([t(1) t(end) -0.5 1.5])

figure;
plot(t,E_hist,'k-','LineWidth',2)
ylabel('Energy [J]')
xlabel('time $t$ [s]')
axis([t(1) t(end) 0 1.5*max(E_hist)])

figure;
subplot(2,2,1)
plot(t,Lb_hist(1,:),'r-',...
     t,Lb_hist(2,:),'g-',...
     t,Lb_hist(3,:),'b-',...
    'LineWidth',2)
legend('$L_x$','$L_y$','$L_z$','Location','Best')
ylabel('Angular momentum $L_\mathrm{b}$ [kg m$^2$/s]')
title('Momentum in body frame')

subplot(2,2,3)
plot(t,Pb_hist(1,:),'r-',...
     t,Pb_hist(2,:),'g-',...
     t,Pb_hist(3,:),'b-',...
    'LineWidth',2)
legend('$P_x$','$P_y$','$P_z$','Location','Best')
ylabel('Linear momentum $P_\mathrm{b}$ [kg m/s]')
xlabel('time $t$ [s]')

subplot(2,2,2)
plot(t,Lw_hist(1,:),'r-',...
     t,Lw_hist(2,:),'g-',...
     t,Lw_hist(3,:),'b-',...
    'LineWidth',2)
legend('$L_x$','$L_y$','$L_z$','Location','Best')
ylabel('Angular momentum $L_\mathrm{w}$ [kg m$^2$/s]')
title('Momentum in world frame')

subplot(2,2,4)
plot(t,Pw_hist(1,:),'r-',...
     t,Pw_hist(2,:),'g-',...
     t,Pw_hist(3,:),'b-',...
    'LineWidth',2)
legend('$P_x$','$P_y$','$P_z$','Location','Best')
ylabel('Linear momentum $P_\mathrm{w}$ [kg m/s]')
xlabel('time $t$ [s]')

figure;
subplot(2,2,1)
plot(t,squeeze(Twb_hist(1,1,:)),'r-',...
     t,squeeze(Twb_hist(1,2,:)),'g-',...
     t,squeeze(Twb_hist(1,3,:)),'b-',...
     t,squeeze(Twb_hist(2,1,:)),'r--',...
     t,squeeze(Twb_hist(2,2,:)),'g--',...
     t,squeeze(Twb_hist(2,3,:)),'b--',...
     t,squeeze(Twb_hist(3,1,:)),'r:',...
     t,squeeze(Twb_hist(3,2,:)),'g:',...
     t,squeeze(Twb_hist(3,3,:)),'b:',...
     'LineWidth',2);
legend('$R_{11}$','$R_{12}$','$R_{13}$',...
       '$R_{21}$','$R_{22}$','$R_{23}$',...
       '$R_{31}$','$R_{32}$','$R_{33}$',...
       'Location','Best')
ylabel('orientation $R_\mathrm{wb}$');

subplot(2,2,3)
plot(t,Vb_hist(1,:),'r-',...
     t,Vb_hist(2,:),'g-',...
     t,Vb_hist(3,:),'b-',...
    'LineWidth',2);
legend('$\omega_x$','$\omega_y$','$\omega_z$','Location','Best')
ylabel('angular velocity [rad/s] in $\left\{\mathrm{b}\right\}$');
xlabel('time $t$ [s]')

subplot(2,2,2)
plot(t,p_hist(1,:),'r-',...
     t,p_hist(2,:),'g-',...
     t,p_hist(3,:),'b-',...
     'LineWidth',2);
legend('$p_x$','$p_y$','$p_z$','Location','Best')
ylabel('COM position in $\left\{\mathrm{w}\right\}$ [m]');

subplot(2,2,4)
plot(t,v_hist(1,:),'r-',...
     t,v_hist(2,:),'g-',...
     t,v_hist(3,:),'b-',...
    'LineWidth',2);
legend('$v_x$','$v_y$','$v_z$','Location','Best')
ylabel('COM velocity [m/s] in $\left\{\mathrm{w}\right\}$');
xlabel('time $t$ [s]')

end