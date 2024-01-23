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
%   tspan: time window [t_init, t_end]; outer bounds on event time
%   Twb_init: configuration of body frame {b} w.r.t. world frame {w} at
%       t_init
%   Vb_init: body twist at t_init
%   err_tol: convergence critertion; absolute value of event function must
%       be less than this number for event_locator to return exitflag = 1
%   max_its: maximum number of iterations the bisection method may take; if
%       bisection method exceeds this limit, event_locator returns exitflag
%       = 0
%
% Outputs:
%   event_time: the time at which the event function first crosses zero
%   exitflag (optional): 1 (success), 0 (exceeded max_its), or -1 (fail)
%   err (optional): value of event function at event function; compare to 0
%   iterations (optional): number of iterations taken by bisection method

function [event_time,varargout] = event_locator(event_fun,...
                                                tspan,...
                                                Twb_init,...
                                                Vb_init,...
                                                err_tol,...
                                                max_its)

end