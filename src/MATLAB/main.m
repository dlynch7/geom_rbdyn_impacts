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
dimensions.ellipsoid.a = 0.1;
dimensions.ellipsoid.b = 0.1;
dimensions.ellipsoid.c = 1;
params = init_params(body_name,dimensions);

%% Set initial conditions and numerical integration parameters
Rwb_init = eye(3);
Twb_init = [Rwb_init,   zeros(3,1);
            zeros(1,3), 1]; % CoM frame coincident with world frame
Twb_init(3,4) = 1; % initial height above ground [m]
Vb_init = [10;10;10;0;0;0];
Fb = zeros(6,1);
        
dt = 1e-3; % timestep size [s]
t0 = 0;
tf = 1; % simulation duration [s]
t = t0:dt:tf;

% event location parameters:
bis_err_tol = 1e-6;
bis_its_max = 100;

%% Simulate via numerical integration
Twb{1} = Twb_init;
Vb{1} = Vb_init;
for k = 1:numel(t)-1
    t_start = t(k);
    t_end = t(k+1);
    [~,Twb{k+1},Vb{k+1}] = cg4(@rbdyn,dt,t(k),Twb{k},Vb{k},Fb,params);
    
    event_time = event_locator(@(T) dist_to_impact(T,params),...
                                [t(k),t(k+1)],...
                                Twb{k},...
                                Vb{k},...
                                bis_err_tol,...
                                bis_its_max);
    
    % check for impact at t and t + dt:
    h_start = dist_to_impact(Twb{k},params);
    h_end = dist_to_impact(Twb{k+1},params);
    
    Twb_start = Twb{k};
    Vb_start = Vb{k};
    
    if sign(h_end) ~= sign(h_start)
        % bisection search:
        bis_err = 10*bis_err_tol;
        bis_its = 0;
        while ((bis_err > bis_err_tol) && (bis_its < bis_its_max))
            t_mid = t_start + 0.5*(t_end - t_start);
            [~,Twb_mid,Vb_mid] = cg4(@rbdyn,t_mid - t_start,t_start,Twb_start,Vb_start,Fb,params);
            h_mid = dist_to_impact(Twb_mid,params);
            
            % update inputs to termination criteria:
            bis_err = abs(h_mid);
            bis_its = bis_its + 1;

            if sign(h_mid) ~= sign(h_start)
%                 fprintf('Iteration %d: impact occurred between t = %f and %f s.\n',bis_its,t_start,t_mid);
                t_end = t_mid;
                h_end = h_mid;
            elseif sign(h_mid) ~= sign(h_end)
%                 fprintf('Iteration %d: impact occurred between t = %f and %f s.\n',bis_its,t_mid,t_end);
                t_start = t_mid;
                h_start = h_mid;
                Twb_start = Twb_mid;
                Vb_start = Vb_mid;
            elseif h_mid == 0
                break;
            else % no event detected
               error('Failed to locate event.');
            end
        end
        fprintf('Impact at t = %f s, %d iterations, error = %f.\n',...
            t_mid,bis_its,bis_err);
    end
end

%% Postprocess
% kinematics:
Twb_hist = NaN(4,4,numel(t));
Vb_hist = NaN(6,numel(t));
p_hist = NaN(3,numel(t));
v_hist = NaN(3,numel(t));

% event function (distance to impact):
dist_hist = NaN(1,numel(t));

for k = 1:numel(t)
    Rwb = Twb{k}(1:3,1:3); % orientation of {b} in {w}
    pwb = Twb{k}(1:3,4);   % position of origin of {b} in {s}
    
    % velocity of CoM in world frame:
    v = Rwb*Vb{k}(4:6);
    
    % distance to impact:
    dist_hist(k) = dist_to_impact(Twb{k},params);
    
    Twb_hist(:,:,k) = Twb{k};
    Vb_hist(:,k) = Vb{k};
    p_hist(:,k) = pwb;
    v_hist(:,k) = v;
end

%% Visualize simulation results
figure;
subplot(2,1,1)
plot(t,p_hist(3,:),'k:',t,dist_hist,'k-','LineWidth',2)
legend('$p_z$','distance to impact','Location','Best')
ylabel('position [m]')
subplot(2,1,2)
plot(t,p_hist(3,:)- dist_hist,'k-','LineWidth',2)
legend('$p_z$ - distance to impact','Location','Best')
ylabel('position [m]')
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