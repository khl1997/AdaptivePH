function errSys = errorSystem(sys1, sys2)
% ERRORSYSTEM - Returns the system model (sss-object!) of the system which
%               subtracts output of sys2 from output of sys1.
%
% Syntax:
%   errSys = ERRORSYTEM(sys1, sys2)
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       - sys1: System 1, "positive", phs-object
%       - sys2: System 2, "negative", phs-object
% 
% Output Arguments:
%       - errSys:   error system ("sys1 - sys2"), sss-object
%
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
% Last Change:  16 September 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
    % calculate dimensions
    dim1 = sys1.dim; dim2 = sys2.dim;
    dim_err = dim1 + dim2;
    s = size(sys1.B);
    
    % initialize error system matrices with zeros
    A_err = zeros(dim_err);
    E_err = zeros(dim_err);
    B_err = zeros(dim_err, s(2));
    C_err = zeros(s(2), dim_err);
    D_err = zeros(s(2), s(2));
    
    % create error system matrices
    A1 = (sys1.J-sys1.R)*sys1.Q;
    A2 = (sys2.J-sys2.R)*sys2.Q;
    A_err(1:dim1, 1:dim1) = A1;
    A_err(dim1+1:dim1+dim2, dim1+1:dim1+dim2) = A2;
    
    E_err(1:dim1, 1:dim1) = sys1.E;
    E_err(dim1+1:dim1+dim2, dim1+1:dim1+dim2) = sys2.E;
    
    B_err(1:dim1, 1:s(2)) = sys1.B-sys1.P;
    B_err(dim1+1:dim1+dim2, 1:s(2)) = sys2.B-sys2.P;
    
    C_err(1:s(2), 1:dim1) = (sys1.B+sys1.P)'*sys1.Q;
    C_err(1:s(2), dim1+1:dim1+dim2) = -(sys2.B+sys2.P)'*sys2.Q;
    
    D_err(1:s(2), 1:s(2)) = sys1.S+sys1.N-sys2.S-sys2.N;
    
    % create error system
    if isequal(E_err, eye(dim_err))
        errSys = ss(A_err,B_err,C_err,D_err);
    else
        errSys = dss(A_err,B_err,C_err,D_err,E_err);
    end
    
end