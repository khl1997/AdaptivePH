function sys = makeExplicit(sys)
% MAKEEXPLICIT - Makes Descriptor system sys explicit by applying
%       state space transformation z = Ex
%
% Syntax:
%   sysExp = MAKEEXPLICIT(sys)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:  phs object
%
% Output Arguments:
%       - sysExp:   phs object of explicit PH system
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
% Last Change:  30 August
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

    if sys.isDescriptor
        if sys.isDAE
            error("Cannot make PHDAE system explicit: E is not invertible!");
        end
        sys.Q = sys.Q/sys.E;   
        sys.E = eye(size(sys.E));
        % This does not change the original system since it is 'call by value' instead of 'call by reference'
        % (phs class does NOT inherit from Matlab 'handle' class!)
    else
        disp("System is already in explicit representation!");
    end
end