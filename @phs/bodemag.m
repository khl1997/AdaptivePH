function bodemag(sys)
% BODEMAG - plots only the magnitude of the bode plot of the 
%           port-Hamiltonian system sys
%
% Syntax:
%   BODEMAG(sys)
%
% Description:
%       bodemags(sys) plots the magnitude part of the bode diagram of the
%       PH-system sys
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:       phs-object
% 
% Output Arguments:
%       /
%
% Examples:
%       The following code creates a simple phs object and runs the
%       bodemag-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       bodemag(sys);
%
% See Also:
%       phs, ss/bodemag
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
    bodemag(dynSys);
end