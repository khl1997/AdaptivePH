function [y,t] = impulse(varargin)
% IMPULSE - plots and returns the impulse response of the 
%           port-Hamiltonian system sys
%
% Syntax:
%   IMPULSE(sys)
%   IMPULSE(sys, tfinal)
%   [y, t] = IMPULSE(sys)
%   [y, t] = IMPULSE(sys, tfinal)
%
% Description:
%       impulse(sys) plots the impulse response of the PH system sys.
%       In case of MIMO-systems, only the impulse response of the first
%       output variable to an impulse on the first input variable is
%       plotted!
%
%       [y, t] = impulse(sys, tfinal) plots the impulse response of the PH
%       system sys and returns the vectors t and y of the system response.
%       It therefore transforms the phs-object into a ss-object and
%       then runs the function ss/impulse.
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys:		phs-object
%       *Optional Input Arguments:*
%       - tfinal:	end time of simulation
% 
% Output Arguments:
%       - y: 		vector of system response (output-variable y)
%       - t: 		vector of time corresponding to vector y
%
% Examples:
%       The following code creates a simple phs object and runs the
%       impulse-function:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); G = [1; 0];
%       sys = phs(J, R, Q, G);
%       impulse(sys);
%
% See Also:
%       phs, ss/impulse
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
        [y,t] = impulse(dynSys);
        plot(t, y(:,1,1));
    elseif nargin == 2
        sys = varargin{1};
        tfinal = varargin{2};
        dynSys = phs2ss(sys);
        [y,t] = impulse(dynSys, tfinal);
        plot(t,y(:,1,1))
    else
        disp('Wrong number of arguments! (Expected 1 or 2)')
    end
end