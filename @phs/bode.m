function bode(sys)
% BODE - This function plots the bode-plot of the port-Hamiltonian system sys
%
% Syntax:
%   BODE(sys)
%
% Description:
%       bode(sys) plots the bode diagram of the port-Hamiltonian system
%       sys. It therefore transforms the phs-object into a ss-object and
%       then runs the function ss/bode.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:       phs-object, system of interest
%
% Output Arguments:
%       /
%
% Examples:
%       The following code creates a simple phs object and runs the
%       bode-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       bode(sys);
%
% See Also:
%   phs, ss/bode
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
% Copyright (c) 2015 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
    dynSys = phs2ss(sys);
    bode(dynSys);
end