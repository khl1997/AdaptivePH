classdef phsRed < phs
% PHSRED (class) - offers same functionalities as phs-class and adds some
%                   features for reduced models
%
% Syntax:
%   redSys = PHSRED(J,R,Q,B)
%   redSys = PHSRED(J,R,Q,B,method)
%   redSys = PHSRED(J,R,Q,B,method,params)
%   redSys = PHSRED(J,R,Q,B,E)
%   redSys = PHSRED(J,R,Q,B,E,method)
%   redSys = PHSRED(J,R,Q,B,E,method,params)
%   redSys = PHSRED(J,R,Q,B,E,P,S,N)
%   redSys = PHSRED(J,R,Q,B,E,P,S,N,method,params)
%
% Description:
%       This class extends the phs-class. It is useful for reduction
%       algorithms: It stores all the attributes and offers all functions
%       provided by the phs-class. In addition, it stores the method and
%       parameters that were used for deriving the reduced model (if
%       information is provided by reduction algorithm)
%
% Input Arguments:
%       *Required Input Arguments:*
%       - J,R,Q,B:  system matrices
%       *Optional Input Arguments:*
%       - method:	method used for reduction
%       - params:	structure with reduction parameters (depends on method)
%       - E,P,S,N:  further system matrices
%
% Output Arguments:
%       - redSys:   phsRed-object
%
% See Also:
%       phs
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
% Last Change:  19 June 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    
    properties (Access = public)
        parameters
        method
    end
    
    methods
        function redSys = phsRed(J, R, Q, B, varargin)
            narginchk(4,10);
            [E, P, S, N, method, params,opts] = phsRed.parseInputs(Q, B, varargin{:});
            redSys = redSys@phs(J,R,Q,B,E,P,S,N,opts);
            redSys.method = method;
            redSys.parameters = params;
        end %of function phsRed (constructor)
    end %methods public
    
    methods (Static)
        function [E, P, S, N, method, params, opts] = parseInputs(Q, B, varargin)
            n_varargin = nargin-2;
            if isempty(varargin)
                E = eye(size(Q,1));
                P = zeros(size(B));
                S = zeros(size(B,2), size(B,2));
                N = zeros(size(B,2), size(B,2));
                method = 'unknown';
                params = struct();
                opts = struct();
                return
            end
            methodProvided = false;
            paramsProvided = false;
            if n_varargin == 1
                if isstring(varargin{end}) || ischar(varargin{end})
                % method as last input argument
                    method = varargin{end};
                    methodProvided = true;
                    params = struct();
                elseif isstruct(varargin{end})
                    params = varargin{end};
                    paramsProvided = true;
                else 
                % matrix as last input argument
                    method = 'unknown';
                    params = struct();
                end
            else            
                if isstring(varargin{end}) || ischar(varargin{end})
                % method as last input argument
                    method = varargin{end};
                    methodProvided = true;
                    params = struct();
                elseif (isstring(varargin{end-1}) || ischar(varargin{end-1})) && isstruct(varargin{end})
                % method and params as last input arguments
                    method = varargin{end-1};
                    methodProvided = true;
                    params = varargin{end};
                    paramsProvided = true;
                else
                % matrix as last input argument
                    method = 'unknown';
                    params = struct();
                end                    
            end
            
            opts = struct();
            if paramsProvided
                if isfield(params, 'inputValidation')
                    opts.inputValidation = params.inputValidation;
                    params = rmfield(params, 'inputValidation');
                end
                if isfield(params, 'inputTolerance')
                    opts.inputTolerance = params.inputTolerance;
                    params = rmfield(params, 'inputTolerance');
                end                
            end
            
            n_matrices = 0;
            if ~methodProvided && ~paramsProvided
                n_matrices = n_varargin;
            elseif methodProvided && ~paramsProvided
                n_matrices = n_varargin-1;
            elseif methodProvided && paramsProvided
                n_matrices = n_varargin-2;
            else
                error('phsmor:phsRed:wrongInput', 'input pattern not supported!');
            end
            
            switch n_matrices
                case 0
                    E = eye(size(Q));
                    P = zeros(size(B));
                    S = zeros(size(B,2),size(B,2));
                    N = zeros(size(B,2),size(B,2));
                case 1
                    E = varargin{1};
                    P = zeros(size(B));
                    S = zeros(size(B,2),size(B,2));
                    N = zeros(size(B,2),size(B,2));
                case 2
                    E = varargin{1};
                    P = varargin{2};
                    S = zeros(size(B,2),size(B,2));
                    N = zeros(size(B,2),size(B,2));
                case 3
                    E = varargin{1};
                    P = varargin{2};
                    S = varargin{3};
                    N = zeros(size(B,2),size(B,2));
                case 4
                    E = varargin{1};
                    P = varargin{2};
                    S = varargin{3};
                    N = varargin{4};
                otherwise
                    error('phsmor:phsRed:programmingError', 'something is programmed wrong... Please contact support ;)');
            end
            
        end %of function parseInputs
    end %methods private
    
end