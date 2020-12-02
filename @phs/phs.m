classdef phs
% PHS (class) - provides a class and methods for handling port-Hamiltonian systems
%
% Syntax:
%   sys = PHS(J, R, Q, B)
%   sys = PHS(J, R, Q, B, E)
%   sys = PHS(J, R, Q, B, opts)
%   sys = PHS(J, R, Q, B, E, opts)
%   sys = phs(J, R, Q, B, E, P, S, N)
%   sys = phs(J, R, Q, B, E, P, S, N, opts)
%
% Description:
%       sys = phs(J, R, Q, B) returns a phs-object representing the
%       port-Hamiltonian system sys defined by
%           dx/dt = (J - R)*Q*x(t) + B*u(t)
%               y = B'*Q*x
%            H(x) = x'*Q*x      (Hamiltonian)
%               f = - dx/dt     (flow)
%               e = grad_x(H)   (effort)        [1,2]
%
%       sys = phs(J, R, Q, B, E) returns a phs-object representing the
%       port-Hamiltonian system sys defined by
%           E*dx/dt = (J - R)*Q*x(t) + B*u(t)
%                 y = B'*Q*x
%              H(x) = x'*E'*Q*x   (Hamiltonian)
%                 f = - dx/dt     (flow)
%                 e = grad_x(H)   (effort)      [3,4]
%
%       sys = phs(J, R, Q, B, E, P, S, N) returns a phs-object representing 
%       the port-Hamiltonian system sys defined by
%           E*dx/dt = (J - R)*Q*x(t) + (B - P)*u(t)
%                 y = (B + P)'*Q*x + (S + N)*u
%              H(x) = x'*E'*Q*x   (Hamiltonian)
%                 f = - dx/dt     (flow)
%                 e = grad_x(H)   (effort)      [3,4]
%
%       For detailed information on matrix properties, see [3, Definition 5].
%
%
% Input Arguments:
%       *Required Input Arguments:*
%       - J:        SKEW-SYMMETRIC matrix describing the interconnection
%                   structure of the system
%       - R:        POSITIVE SEMI-DEFINITE SYMMETRIC matrix describing the
%                   dissipative behaviour of the system
%       - Q:        POSITIVE SEMI-DEFINITE, SYMMETRIC matrix describing the
%                   energy storage behaviour of the system. (For
%                   descriptor-systems, restrictions differ slightly. -> see
%                   'E')
%       - B:        input/output matrix describing the
%                   input/output-behaviour of the system.
%       *Optional Input Arguments:*
%       - E:        (non-)singular descriptor matrix (Q'*E must be
%                   symmetric and positive semi-definite)
%       - P:        input/output matrix describing the
%                   input/output-behaviour of the system
%       - S:        SYMMETRIC feedthrough matrix
%       - N:        SKEW-SYMMETRIC feedthrough matrix
%       - opts:     structure with options:
%           - .inputValidation:     switch input validation on (true) or off
%                                   (false); This may speed up instantiation
%                                   but faulty input is not detected any
%                                   more. Default: true
%           - .inputTolerance:      Tolerance for input validation. If,
%                                   for example, Q has eigenvalues close to
%                                   zero in the left half plane, it will
%                                   still be considered positive definite
%                                   if their (negative) real part is within
%                                   opts.inputTolerance. Default: 1e-15
% 
% Output Arguments:
%       - sys:      phs-object
%
% Examples:
%       The following code creates a simple phs object:
%
%       J = [0 -1; 1 0]; R = [0 0; 0 1]; Q = eye(2); B = [1; 0];
%       sys = phs(J, R, Q, B);
%
% See Also:
%       ss, sss
%
% References:
%       [1] Duindam et al. (2009), 
%               "Modeling and Control of Complex Physical Systems - The
%               Port-Hamiltonian approach"
%       [2] van der Schaft (2017), 
%               "L2-Gain and Passivity Techniques in Nonlinear Control"
%       [3] Beattie et al. (2018), 
%               "Linear port-Hamiltonian descriptor systems"
%       [4] Mehl et al. (2018), 
%               "Linear algebra properties of dissipative Hamiltonian
%               descriptor systems"
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
% Last Change:  12 September 2020
% Copyright (c) 2020 Chair of Automatic Control, TU Muenchen
%------------------------------------------------------------------
    %% properties
    properties (Access = public)
        J
        R
        Q
        B
        E
        P
        S
        N
        opts
    end
    properties (Access = private)
        flag_descriptor
        flag_MIMO
        flag_DAE
    end
    properties (Dependent)
        dim (1,1)
    end
    %end of properties
    
    %% methods
    methods
        
        function sys = phs(J, R, Q, B, varargin)
         %creates phs object
         
            % empty constructor (required by MATLAB)
            if nargin == 0
                sys.J = 0; sys.R = 0; sys.Q = 0; sys.B = 0; sys.E = 0; 
                sys.P = 0; sys.S = 0; sys.N = 0; 
                sys.opts = struct();
                sys.opts.inputValidation = true;
                sys.opts.inputTolerance = 1e-15;
                return
            end
            
            narginchk(4,9);
            
            [J, R, Q, B, E, P, S, N, opts] = phs.parseInputs(J, R, Q, B, varargin{:});
                   
            if opts.inputValidation
                phs.inputValidation(J,R,Q,B,E,P,S,N,opts);
                % inputValidation function is static (see methods (Static))
                % and specified in external file in @phs class folder
            else
                warning('No input validation for correctness of port-Hamiltonian system!')
            end
                
            % save system matrices
            sys.opts = opts;
            sys.J = J;
            sys.R = R;
            sys.Q = Q;
            sys.B = B;
            sys.E = E;
            sys.P = P;
            sys.S = S;
            sys.N = N;
            % set flags
            if isequal(E, eye(length(J)))
                sys.flag_descriptor = false;
                sys.flag_DAE = false;
            else
                sys.flag_descriptor = true;
                if rank(E) < length(E)
                    sys.flag_DAE = true;
                else
                    sys.flag_DAE = false;
                end
            end
            if size(B,2) > 1
                sys.flag_MIMO = true;
            else
                sys.flag_MIMO = false;
            end
        end
        
        
        %% Getter and Setter
                
        function J = get.J(sys); J = sys.J; end
        function R = get.R(sys); R = sys.R; end
        function Q = get.Q(sys); Q = sys.Q; end
        function B = get.B(sys); B = sys.B; end
        function E = get.E(sys); E = sys.E; end
        function P = get.P(sys); P = sys.P; end
        function S = get.S(sys); S = sys.S; end
        function N = get.N(sys); N = sys.N; end
        function dim = get.dim(sys); dim = length(sys.J); end
        
        function sys = set.J(sys, J)
            if ~isempty(sys.J)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.J = J;
        end
        function sys = set.R(sys, R)
            if ~isempty(sys.R)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.R = R; 
        end
        function sys = set.Q(sys, Q)
            if ~isempty(sys.Q)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.Q = Q;
        end
        function sys = set.B(sys, B)
            if ~isempty(sys.B)
            warning('phs:phs:ChangedProperty',...
                'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.B = B;
            if size(B,2) > 1
                sys.flag_MIMO = true;
            else
                sys.flag_MIMO = false;
            end
        end
        function sys = set.E(sys, E)
            if ~isempty(sys.E)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.E = E;  
            if isequal(E, eye(sys.dim))
                sys.flag_descriptor = false;
                sys.flag_DAE = false;
            else
                sys.flag_descriptor = true;
                if rank(E) < length(E)
                    sys.flag_DAE = true;
                else
                    sys.flag_DAE = false;
                end
            end
        end
        function sys = set.P(sys, P)
            if ~isempty(sys.P)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.P = P;
        end
        function sys = set.S(sys, S)
             if ~isempty(sys.S)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
             end
             sys.S = S; 
        end
        function sys = set.N(sys, N)
            if ~isempty(sys.N)
                warning('phs:phs:ChangedProperty',...
                    'Property has been changed successfully.\nConsider running phs.inputValidation(sys,sys.opts) for validating new properties');
            end
            sys.N = N; 
        end
        function sys = set.dim(sys, dim); warning('You cannot set the dimension of the system explicitly!'); end
        
        function [J, R, Q, B, E, P, S, N] = getMatrices(sys)
        % GETMATRICES - returns the matrices of the system
        %   Syntax:
        %       [J, R, Q, B, E, P, S, N] = GETMATRICES(sys)
            J = sys.J;
            R = sys.R;
            Q = sys.Q;
            B = sys.B;
            E = sys.E;
            P = sys.P;
            S = sys.S;
            N = sys.N;
        end
        
        function isDescriptor = isDescriptor(sys)
            isDescriptor =  sys.flag_descriptor;
        end
        
        function isMIMO = isMIMO(sys)
            isMIMO = sys.flag_MIMO;
        end
        
        function isDAE = isDAE(sys)
            isDAE = sys.flag_DAE;
        end
        
        %% Converter
        
        function sys = phs2ss(sys)
        % PHS2SS - transforms a phs system to regular ss system
        %   Syntax: 
        %       sys_ss = phs2ss(sys)
        %   Input arguments:    sys - phs-object
        %   Output arguments:   sys_ss - ss-object
            if ~sys.flag_descriptor
                sys = ss((sys.J - sys.R)*sys.Q, sys.B-sys.P, (sys.B+sys.P)'*sys.Q, sys.S+sys.N);
                %sys = ss((sys.J - sys.R)*sys.Q, sys.B, sys.B'*sys.Q, 0);
            else
                disp('System is descriptor system')
                sys = dss((sys.J - sys.R)*sys.Q, sys.B-sys.P, (sys.B+sys.P)'*sys.Q, sys.S+sys.N, sys.E);
            end
        end
        
        function sys = ss(sys)
        % SS - transforms a phs system to regular ss system
        %   Syntax: 
        %       sys_ss = ss(sys)
        %   Input arguments:    sys - phs-object
        %   Output arguments:   sys_ss - ss-object
            sys = phs2ss(sys);
        end
        
        function sys = phs2sss(sys)
        % PHS2SSS - transforms a phs system to a sss system
        %   Syntax: 
        %       sys_sss = phs2sss(sys)
        %   Input arguments:    sys - phs-object
        %   Output arguments:   sys_sss - sss-object
            sys = sss((sys.J - sys.R)*sys.Q, sys.B-sys.P, (sys.B+sys.P)'*sys.Q, sys.S+sys.N, sys.E);
        end
        
        function sys = sss(sys)
        % SSS - transforms a phs system to a sss system
        %   Syntax: 
        %       sys_sss = sss(sys)
        %   Input arguments:    sys - phs-object
        %   Output arguments:   sys_sss - sss-object
            sys = phs2sss(sys);
        end
        
        %% Operator overloading
        
        function sys = minus(sys1,sys2)
            sys = errorSystem(sys1,sys2);
            %Alternative: sys = ss(sys1) - ss(sys2)
        end
        
        %% Methods in external files
        [y, t] = step(varargin);        %step response
        [y, t] = impulse(varargin);     %impulse response
        bode(sys);                      %bode plot
        bodemag(sys);                   %bode plot (magnitude only)
        disp_(sys);                     %display information on system
        [p, z] = pzmap(sys);            %poles and zeros of the system
        y = lsim(sys, u, t);            %response to arbitrary input
        varargout = norm(varargin);     %norm (H2, H-inf)
        sys = makeExplicit(sys);        %make Descriptor-system explicit by state space transformation
        
    end
    
    methods (Access = private)
        
    end
    
    methods (Static)
        
        function [J, R, Q, B, E, P, S, N, opts] = parseInputs(J, R, Q, B, varargin)
            narginchk(4, 9)
            s = size(B);
            
            if ~isempty(varargin) && isstruct(varargin{end})
                opts = varargin{end};
                if ~isfield(opts, 'inputValidation')
                    opts.inputValidation = true;
                    disp('opts.inputValidation has been set to true')
                end
                if ~isfield(opts, 'inputTolerance')
                    opts.inputTolerance = 1e-14;
                end
                switch length(varargin)
                    case 1
                        E = eye(length(J));
                        P = zeros(s);
                        S = zeros(s(2),s(2));
                        N = zeros(s(2),s(2));
                    case 2
                        E = varargin{1};
                        P = zeros(s);
                        S = zeros(s(2),s(2));
                        N = zeros(s(2),s(2));
                    case 5
                        E = varargin{1};
                        P = varargin{2};
                        S = varargin{3};
                        N = varargin{4};
                    otherwise
                        error('phsmor:phs:parseInput', 'Input combination is not specified!')
                end
            else
                opts.inputValidation = true;
                opts.inputTolerance = 1e-14;
                switch length(varargin)
                    case 0
                        E = eye(length(J));
                        P = zeros(s);
                        S = zeros(s(2),s(2));
                        N = zeros(s(2),s(2));
                    case 1
                        E = varargin{1};
                        P = zeros(s);
                        S = zeros(s(2),s(2));
                        N = zeros(s(2),s(2));
                    case 4
                        E = varargin{1};
                        P = varargin{2};
                        S = varargin{3};
                        N = varargin{4};
                    otherwise
                        error('phsmor:phs:parseInput', 'Input combination is not specified!')
                end
            end
        end
                
        %% Methods in external files
        result = isPositiveDefiniteHermitian(A,varargin);    %returns 1 if hermitian A is positive definite
        isPH = inputValidation(varargin);                   % check for PH structure of the system
        
    end
    %end of methods
    
end