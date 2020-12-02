function [y, t] = step(varargin)
% STEP - plots and returns the step-response of the port-Hamiltonian system
%        sys
%
% Syntax:
%   STEP(sys)
%   STEP(sys,tfinal)
%   [y, t] = STEP(sys)
%   [y, t] = STEP(sys, tfinal)
%
% Description:
%       step(sys) plots the step-response of the phs-object sys.
%       In case of MIMO-systems only the first output stimulated by a step
%       in the first input is plotted!
%
%       [y, t] = step(sys, tfinal) plots and returns the step-response of
%       the system sys calculated to tfinal.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs-object
%       *Optional Input Arguments:*
%       - tfinal:   end time of simulation
%
% Output Arguments:
%       - y:        vector/matrix of step-response
%       - t:        vector of time corresponding to y
%
% Examples:
%       The following code creates a simple phs object and runs the
%       step-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       step(sys);
%
% See Also:
%       phs, ss/step
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
% Last Change:  12 May 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    if nargin == 1
        sys = varargin{1};
        dynSys = phs2ss(sys);
        [y, t] = step(dynSys);
        plot(t, y(:,1,1));
    elseif nargin == 2
        sys = varargin{1};
        tfinal = varargin{2};
        dynSys = phs2ss(sys);
        [y, t] = step(dynSys, tfinal);
        plot(t,y(:,1,1))
    else
        disp('Wrong number of arguments! (Expected 1 or 2)')
    end
end