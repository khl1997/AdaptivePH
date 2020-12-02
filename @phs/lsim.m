function y = lsim(sys, u, t)
% LSIM - simulates the port-Hamiltonian system sys for given input u and
%        time t
%
% Syntax:
%   y = LSIM(sys, u, t)
%
% Description:
%       y = lsim(sys, u, t) simulates the system for input u
%       and returns the system response y corresponding to the time 
%       vector t
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs-object
%       - u:        input vector
%       - t:        time vector
%
% Output Arguments:
%       - y: 		system response
%
% Examples:
%       The following code creates a simple phs object and runs the
%       lsim-function for a sine input signal:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       t = linspace(0, 100, 1000);
%       u = sin(2*pi*0.1*t);
%       lsim(sys, u, t);
%
% See Also:
%       phs, lsim
%
%------------------------------------------------------------------
% This file is part of <a href="matlab:docsearch phsmor">phsmor</a>, a port-Hamiltonian system and reduction
% Toolbox developed at the Chair of Automatic Control, Technische Universitaet Muenchen.
% For updates and further information please visit <a href="https://www.rt.mw.tum.de/?phsmor">www.rt.mw.tum.de/?phsmor</a>
% For any suggestions, submission and/or bug reports, mail us at
%                   -> <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a> <-
%
% More Toolbox Info by searching <a href="matlab:docsearch phsmor">phsmor</a> in the Matlab Documentation
%
%------------------------------------------------------------------
% Authors:      Julius Durmann
% Email:        <a href="mailto:morlab@rt.mw.tum.de">morlab@rt.mw.tum.de</a>
% Website:      <a href="https://www.rt.mw.tum.de/?phsmor">www.rt.mw.tum.de/?phsmor</a>
% Work Adress:  Technische Universitaet Muenchen
% Last Change:  13 May 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    dynSys = phs2ss(sys);
    y = lsim(dynSys, u, t);
    plot(t,y(:,1,1));
end