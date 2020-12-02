function result = isPositiveDefiniteHermitian(A,tol,varargin)
% ISPOSITIVEDEFINITEHERMITIAN - Determines whether the matrix A is positive
%       (semi-) definite. Uses eig or eigs (depending on dimension of
%       matrix A) to determine the smallest realpart of the eigenvalues of
%       the eigenvalues of A. If this value is smaller than zero - tol, the
%       function will return false, otherwise it will return true.
%
% Syntax:
%   result = ISPOSITIVEDEFINITEHERMITIAN(A, tol)
%   result = ISPOSITIVEDEFINITEHERMITIAN(A, tol, 'semi')
%
% Description:
%       result = isPositiveDefiniteHermitian(A) returns logical value
%       'result' which is 1 for positive definite A
%
%       result = isPositiveDefiniteHermitian(A, 'semi') returns logical
%       value 'result' which is 1 for positive semi-definite A
%
% Input Arguments:
%       *Required Input Arguments:*
%       - A:		hermitian matrix
%       *Optional Input Arguments:*
%       - 'semi':   qualifier if semi-definiteness is requested
%
% Output Arguments:
%       - result:   logical value
%
%
% See also:
%       eig, eigs
%
% References:
%       [1] documentation on eigs
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
% Last Change:  3 July 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------

% alternative to  Sylvester's criterion: check if there exists a Cholesky
% decomposition

    if length(A) > 5e3
        use = 'eigs'; %'eig' or 'eigs'
    else
        use = 'eig';
    end

    narginchk(2,3);
    result = true;
    
    switch use
        case 'eig'
            x = min(real(eig(A)));
            if nargin == 2
            % positive definite?
               if isnan(x)
                    warning('Could not make sure that matrix is positive definite (eig threw NaN)!')
                else
                    if x <= 0 - tol
                        result = false;
                    elseif x<= 0
                        warning(['Matrix might not be completely positive definite, but within tolerance (-> numerical accuracy); smallest calculated eigenvalue: ', num2str(x)]);
                    end
                end
            else
            % positive semidefinite?
                if isnan(x)
                    warning('Could not make sure that matrix is positive semidefinite (eig threw NaN)!')
                else
                    if x < 0 - tol
                        result = false;
                    elseif x < 0
                        warning('phs:isPositiveDefiniteHermitian:semiDefinitenessUnsure',['Matrix might not be completely positive semidefinite (computed eigenvalue is smaller than zero), but within tolerance (-> numerical accuracy); smallest calculated eigenvalue: ', num2str(x)]);
                    end
                end

            end
            
        case 'eigs'
            if nargin == 2
            % positive definite?
                x = eigs(A,1,'smallestreal');
                if isnan(x)
                    for i = 1:length(A)
                        d = det(A(1:i,1:i));
                        if isnan(d)
                            result = true;
                            warning('Could not make sure that matrix is positive definite due to failed convergence of eigs and NaN as determinant!');
                            return
                        end
                        if ~(d > 0)
                            result = false;
                            return
                        end
                    end
                else
                    if x <= 0 - tol
                        result = false;
                    elseif x<= 0
                        warning('Matrix might not be completely positive definite, but within tolerance (-> numerical accuracy)');
                    end
                end
            else
            % positive semidefinite?
                x = eigs(A,1,'smallestreal');
                if isnan(x)
                    warning('Could not make sure that matrix is positive semidefinite due to failed convergence of eigs!')
                else
                    if x < 0 - tol
                        result = false;
                    elseif x < 0
                        warning('phs:isPositiveDefiniteHermitian:semiDefinitenessUnsure','Matrix might not be completely positive semidefinite (computed eigenvalue is smaller than zero), but within tolerance (-> numerical accuracy)');
                    end
                end
            end
    end
    
    if result == false
        warning(['Calculated eigenvalue is ', num2str(x)]);
    end

end