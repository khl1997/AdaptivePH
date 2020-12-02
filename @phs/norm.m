function [n, FPEAK] = norm(sys, type, varargin)
% NORM - calculates the H2 or L-infinity norm of the system sys
%
% Syntax:
%   n = NORM(sys)
%   n = NORM(sys, 2)
%   n = NORM(sys, inf)
%   n = NORM(sys, inf, TOL)
%   [GPEAK, FPEAK] = norm(sys, inf)
%
% Description:
%       phs/norm is a wrapper for DynamicSystem/norm!
%
%       n = norm(sys) returns the H2-norm of sys
%
%       n = norm(sys, 2) / n = norm(sys, inf)  returns the H2- or the
%       H-infinity norm respectively
%
%       n = norm(sys, inf, TOL) specifies a relative accuracy TOL for the
%       computed infinity norm (default: TOl = 1e-2)
%
%       [GPEAK,FPEAK] = norm(sys, inf) returns the peak value GPEAK at the
%       corresponding frequency FPEAK
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:      phs-object
%       *Optional Input Arguments:*
%       - type: 	specify which norm should be computed
%       - TOL:		tolerance (default: 1e-2)
%
% Output Arguments:
%       - n:                norm
%       - [GPEAK,FPEAK]: 	see 'Description'
%
% Examples:
%       The following code creates a simple phs object and runs the
%       norm-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       norm(sys);
%
% See Also:
%       DynamicSystem/norm, phs
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
% Last Change:  20 May 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

narginchk(1,3)
if nargin < 2
    type = 2;
end

sys = ss(sys);      %transfer to ss system for later call of DynamicSystem/norm
if nargout == 0 || nargout == 1
    n = norm(sys, type, varargin);
elseif nargout == 2
    [n, FPEAK] = norm(sys, inf);
else
    error('wrong output')
end
end