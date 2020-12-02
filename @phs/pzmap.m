function [p, z] = pzmap(sys)
% PZMAP - plots and returns the poles and zeros of the port-Hamiltonian
%         system sys
%
% Syntax:
%   PZMAP(sys)
%   [p, z] = PZMAP(sys)
%
% Description:
%       pzmap(sys) creates a plot of the poles and zeros of the PH system
%       sys
%
%       [p, z] = pzmap(sys) returns and plots the poles and zeros of the PH
%       system sys
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs-object
% 
% Output Arguments:
%       - p:        poles of sys
%       - z:        zeros of sys
%
% Examples:
%       The following code creates a simple phs object and runs the
%       pzmap-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       pzmap(sys);
%
% See Also:
%       phs, ss/pzmap
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

    dynSys = phs2ss(sys);
    [p, z] = pzmap(dynSys);
    pzmap(dynSys);
end