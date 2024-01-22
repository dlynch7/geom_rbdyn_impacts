%% init_params.m
%
% Description:
%   This function initializes a number of parameters required to simulate
%   rigid-body dynamics.
%
% Inputs:
%   body_name: a string specifying the rigid body type. Options are
%       'rect_prism'
%       'cylinder'
%       'ellipsoid'
%   dimensions: a struct with entries specifying the dimensions of the
%       rigid body. Entry names are specific to each rigid body.
%
% Outputs:
%   params: a struct with many fields

function params = init_params(body_name,dimensions)

params.dyn.grav = 9.81; % acceleration [m/s^2] due to gravity
params.dyn.mass = 1; % rigid body mass [kg]

switch body_name
    case 'rect_prism'
        % extract dimensions from input struct:
        l = dimensions.rect_prism.l; % length [m]
        w = dimensions.rect_prism.w; % width  [m]
        h = dimensions.rect_prism.h; % height [m]
        
        % inertia [kg*m^2] about each principal axis:
        Ixx = params.dyn.mass*(w^2 + h^2)/12;
        Iyy = params.dyn.mass*(l^2 + h^2)/12;
        Izz = params.dyn.mass*(l^2 + w^2)/12;
    case 'cylinder'
        % extract dimensions from input struct:
        r = dimensions.cylinder.r; % radius [m]
        h = dimensions.cylinder.h; % height [m]
        
        % inertia [kg*m^2] about each principal axis:
        Ixx = params.dyn.mass*(3*(r^2) + h^2)/12;
        Iyy = Ixx;
        Izz = params.dyn.mass*(r^2)/2;
    case 'ellipsoid'
        % geometry and inertia components for an ellipsoid:
        a = dimensions.ellipsoid.a; % ellispoid x-axis length [m]
        b = dimensions.ellipsoid.b; % ellispoid y-axis length [m]
        c = dimensions.ellipsoid.c; % ellispoid z-axis length [m]
        Ixx = params.dyn.mass*(b^2 + c^2)/5;
        Iyy = params.dyn.mass*(a^2 + c^2)/5;
        Izz = params.dyn.mass*(a^2 + b^2)/5;
    otherwise
        error("Unrecognized input for 'body_name' argument.");
end

% inertia matrix (3x3) for center-of-mass (CoM) frame aligned with
% principal axes of inertia:
params.dyn.inertia_matrix = diag([Ixx,Iyy,Izz]);

% 6x6 spatial inertia matrix for CoM frame aligned with principal axes of
% inertia:
params.dyn.spatial_inertia_matrix_CoM = [params.dyn.inertia_matrix, zeros(3);
                                         zeros(3), params.dyn.mass*eye(3)];

end