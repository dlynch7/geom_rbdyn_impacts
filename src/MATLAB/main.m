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
body_name = 'ellipsoid';
dimensions.ellipsoid.a = 0.05;
dimensions.ellipsoid.b = 0.05;
dimensions.ellipsoid.c = 0.1;
params = init_params(body_name,dimensions);

%% Set initial conditions and numerical integration parameters
Rwb_init = eye(3);
% Rwb_init = [0 0 1;
%             0 1 0;
%            -1 0 0];
Twb_init = [Rwb_init,   zeros(3,1);
            zeros(1,3), 1]; % CoM frame coincident with world frame
Twb_init(3,4) = 0.5; % initial height above ground [m]
% Vb_init = [10;10;10;0;0;0];
Vb_init = [1;1;1;-1;0;0];
Fb = zeros(6,1);
        
dt = 1e-2; % timestep size [s]
t0 = 0;
tf = 5; % simulation duration [s]
t = t0:dt:tf;

% event location parameters:
bis_err_tol = 1e-6;
bis_its_max = 100;

ode_fun = @(t,T,V) rbdyn(t,T,V,Fb,params);
event_fun = @(T) dist_to_impact(T,params);

%% Simulate via numerical integration
Twb{1} = Twb_init;
Vb{1} = Vb_init;
for k = 1:numel(t)-1
    [~,Twb{k+1},Vb{k+1}] = cg4(ode_fun,dt,t(k),Twb{k},Vb{k});
    
    [t_e,Twb_e,Vb_e] = event_locator(event_fun,...
                                     ode_fun,...
                                     t(k),t(k+1),...
                                     Twb{k},Twb{k+1},...
                                     Vb{k},Vb{k+1},...
                                     bis_err_tol,bis_its_max);
                        
    if ~isempty(t_e)
        fprintf('Impact at t = %f s\n',t_e);
        
        % express contact point in {w}:
        pwc = closest_point_on_ellipsoid(Twb_e,params);
        
        % construct semi-arbitrary Rwc s.t. z-axis of {c} is unit normal of
        % plane (angle about unit normal doesn't matter:
        Rwc = MatrixExp3(VecToso3(params.geom.plane.unit_normal));
        
        % combine Rwc and pwc into Twc:
        Twc = RpToTrans(Rwc,pwc);
        
        [Twb_post,Vb_post,impulse_mag] = impact_map(Twb_e,...
                                                    Vb_e,...
                                                    Twc,...
                                                    params);
                                              
        % simulate to end of timestep using (Twb_post,Vb_post) as initial
        % conditions:
        [~,Twb{k+1},Vb{k+1}] = cg4(ode_fun,t(k+1)-t_e,t_e,Twb_post,Vb_post);
    end
end

%% Animate simulation
options.trace_CoM = true;
options.record_video = true;
options.video_name = 'bounce_with_infinite_friction.avi';
animate(Twb,t,params,options);

%% Postprocess
% kinematics:
Twb_hist = NaN(4,4,numel(t));
Vb_hist = NaN(6,numel(t));
p_hist = NaN(3,numel(t));
v_hist = NaN(3,numel(t));

% event function (distance to impact):
dist_hist = NaN(1,numel(t));

Gb = params.dyn.spatial_inertia_matrix_CoM;

% energy and momenta:
E_hist = NaN(1,numel(t)); % mechanical energy
Lb_hist = NaN(3,numel(t)); % angular momentum in body frame
Pb_hist = NaN(3,numel(t)); % linear momentum in body frame
Lw_hist = NaN(3,numel(t)); % angular momentum in world frame
Pw_hist = NaN(3,numel(t)); % linear momentum in world frame

for k = 1:numel(t)
    Rwb = Twb{k}(1:3,1:3); % orientation of {b} in {w}
    pwb = Twb{k}(1:3,4);   % position of origin of {b} in {s}
    
    % velocity of CoM in world frame:
    v = Rwb*Vb{k}(4:6);
    
    % momentum:
    Pb = Gb*Vb{k};
    Lb_hist(:,k) = [eye(3),zeros(3)]*Pb; % angular momentum in {b}
    Pb_hist(:,k) = [zeros(3),eye(3)]*Pb; % linear momentum in {b}
    Tbw = TransInv(Twb{k});
    Pw = transpose(Adjoint(Tbw))*Pb;
    Lw_hist(:,k) = [eye(3),zeros(3)]*Pw; % angular momentum in {w}
    Pw_hist(:,k) = [zeros(3),eye(3)]*Pw; % linear momentum in {w}
    
    % energy
    KE = 0.5*transpose(Vb{k})*Gb*Vb{k};
    PE = params.dyn.mass*params.dyn.grav*pwb(3);
    E = KE + PE;
    
    % distance to impact:
    dist_hist(k) = dist_to_impact(Twb{k},params);
    
    Twb_hist(:,:,k) = Twb{k};
    Vb_hist(:,k) = Vb{k};
    p_hist(:,k) = pwb;
    v_hist(:,k) = v;
    E_hist(k) = E;
end

%% Visualize simulation results
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
subplot(2,1,1)
plot(t,p_hist(3,:),'k:',t,dist_hist,'k-','LineWidth',2)
legend('$p_z$','bottom of ellipsoid','Location','Best')
ylabel('position [m]')
subplot(2,1,2)
plot(t,p_hist(3,:)- dist_hist,'k-','LineWidth',2)
legend('$p_z$ - bottom of ellipsoid','Location','Best')
ylabel('distance [m]')
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