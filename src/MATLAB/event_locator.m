%% event_locator.m
%
% Description:
%   This function uses bisection to search for the time at which an event
%   occurs (a zero-crossing of a user-defined event function) some time in
%   tspan.
%   Although bisection may require more iterations than a Newton-Raphson
%   root-finder, bisection respects the constraint that the event time be
%   within tspan, whereas Newton-Raphson may converge to an event time
%   outside tspan if the underlying event function is sufficiently
%   pathological.
%
% Inputs:
%   event_fun: handle to the function whose zero-crossing we are finding
%   ode_fun:
%   t_init: lower bound on event time
%   t_term: upper bound on event time
%   Twb_init: configuration of body frame {b} w.r.t. world frame {w} at
%       t_init
%   Twb_term: configuration of body frame {b} w.r.t. world frame {w} at
%       t_term
%   Vb_init: body twist at t_init
%   Vb_term: body twist at t_term
%   err_tol: convergence critertion; absolute value of event function must
%       be less than this number for event_locator to return exitflag = 1
%   max_its: maximum iterations the bisection method may take; if bisection
%       method exceeds this limit, event_locator returns exitflag = 0
%
% Outputs:
%   t_e: the time at which the event function first crosses zero
%   Twb_e: configuration of body frame {b} w.r.t. world frame {w} at t_e
%   Vb_e: body twist at t_e
%   exitflag (optional): 1 (success), 0 (exceeded max_its), or -1 (fail)
%   err (optional): value of event function at event function; compare to 0
%   iterations (optional): number of iterations taken by bisection method

function [t_e,varargout] = event_locator(event_fun,...
                                                odefun,...
                                                t_init,...
                                                t_term,...
                                                Twb_init,...
                                                Twb_term,...
                                                Vb_init,...
                                                Vb_term,...
                                                err_tol,...
                                                max_its)
                                            
    % initialize outputs:
    t_e = [];
    Twb_e = Twb_init;
    Vb_e = Vb_init;
    exitflag = [];
    err = [];
    its = [];

    % check for impact at t and t + dt:
    h_init = feval(event_fun,Twb_init);
    h_term = feval(event_fun,Twb_term);
    
    if sign(h_term) ~= sign(h_init)
        % bisection search:
        err = 10*err_tol;
        its = 0;
        while ((err > err_tol) && (its < max_its))
            % calculate state (Twb, Vb) and event function (h) midway
            % between t_init and t_term:
            t_mid = t_init + 0.5*(t_term - t_init);
            dt = t_mid - t_init;
            [~,Twb_mid,Vb_mid] = cg4(odefun,dt,t_init,Twb_init,Vb_init);
            h_mid = feval(event_fun,Twb_mid);
            
            % update inputs to termination criteria:
            err = abs(h_mid);
            its = its + 1;

            % check for zero crossing in between t_init and t_mid and
            % between t_mid and t_term:
            if sign(h_mid) ~= sign(h_init)
                t_term = t_mid;
                h_term = h_mid;
            elseif sign(h_mid) ~= sign(h_term)
                t_init = t_mid;
                h_init = h_mid;
                Twb_init = Twb_mid;
                Vb_init = Vb_mid;
            elseif h_mid == 0
                break;
            else
                continue;
            end
        end
        t_e = t_mid;
        Twb_e = Twb_mid;
        Vb_e = Vb_mid;
        
         % set exitflag:
        if ((err < err_tol) && (its < max_its))
            exitflag = 1; % success
        elseif (its >= max_its)
            exitflag = 0; % reached iteration limit
        else
            exitflag = -1;
        end
    end
    
    varargout{1} = Twb_e;
    varargout{2} = Vb_e;
    varargout{3} = exitflag;
    varargout{4} = err;
    varargout{5} = its;
                                            
end