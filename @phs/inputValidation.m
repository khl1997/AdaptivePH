function isPH = inputValidation(varargin)
% INPUTVALIDATION - checks if the system sys (alternatively the
%       passed matrices) is / define a PH system
%
% Syntax:
%   isPH = INPUTVALIDATION(sys)
%   isPH = INPUTVALIDATION(J, R, Q, B, E, P, S, N, opts)
% 
% Input arguments:
%   - sys:                      phs object
%   - opts:                     Struct that defines inputTolerance
%   - J, R, Q, B, E, P, S, N:   system matrices
% Output arguments:
%   - isPH:                     true (system is PH system) 
%                               false (system is not PH) (should
%                               normally throw error)


    if nargin > 2
    % -------------------------------------------------------------------------------------
    % check for matrices before system is created
    opts = varargin{end};
        % check if dimensions are correct
            sJ = size(varargin{1});
            % J not square
            if (sJ(1) ~= sJ(2))
                error('phsmor:phs:wrongInput', 'dimensions of J are not correct');
            end
            sR = size(varargin{2});
            % R not square or dimension different to J
            if ( (sR(1) ~= sR(2)) || ~isequal(sJ, sR) )     
                error('phsmor:phs:wrongInput', 'dimensions of R are not correct');
            end
            sQ = size(varargin{3});
            % Q not square or dimension different to J
            if ( (sQ(1) ~= sQ(2)) || ~isequal(sJ, sQ) )     
                error('phsmor:phs:wrongInput', 'dimensions of Q are not correct');
            end
            sE = size(varargin{5});
            % E not square or dimension different to J
            if ( (sE(1) ~= sE(2)) || ~isequal(sJ, sE) )     
                error('phsmor:phs:wrongInput', 'dimensions of E are not correct');
            end
            sB = size(varargin{4});
            % B has wrong state dimension
            if (sB(1) ~= sJ(1))             
                error('phsmor:phs:wrongInput', 'dimensions of B are not correct');
            end
            % different dimensions of P and B
            if (~isequal(sB, size(varargin{6})))      
                error('phsmor:phs:wrongInput', 'dimensions of P are not correct');
            end
            % dimensions of S and N not equal
            if (~isequal(size(varargin{8}), size(varargin{7}))) 
                error('phsmor:phs:wrongInput', 'dimensions of S and N must be equal');
            end
            % dimensions of S and N not correct (compared to B)
            if (size(varargin{8},1) ~= sB(2) || size(varargin{8},2) ~= sB(2))   
                error('phsmor:phs:wrongInput', 'dimensions of S and N are not correct');
            end

        % check if E, J, R, Q, S and N define correct PH system
        % see [1], definition 5
            errorMessage = "\nConsider changing opts.inputTolerance to avoid this error when caused by numerical inaccuracy";

            % check symmetry of S
            if (~issymmetric(varargin{7}))
                norm_S = norm(varargin{7}-varargin{7}')/norm(varargin{7});
                if (norm_S > opts.inputTolerance)
                    errorMessage = strcat("S is not symmetric:\n norm(S-S')/norm(S) = ", num2str(norm_S),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            % check skew-symmetry of N
            if (~issymmetric(varargin{8},'skew'))
                norm_N = norm(varargin{8}+varargin{8}')/norm(varargin{8});
                if  (norm_N > opts.inputTolerance)
                    errorMessage = strcat("N is not skew-symmetric:\n norm(N+N')/norm(N) = ", num2str(norm_N),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            % check symmetry of Q'*E
            if (~issymmetric(varargin{3}'*varargin{5}))
                % norm_QE = norm(Q'*E - E'*Q)/norm(Q'*E)
                norm_QE = norm(varargin{3}'*varargin{5} - varargin{5}'*varargin{3})/norm(varargin{3}'*varargin{5});
                if (norm_QE > opts.inputTolerance)
                    errorMessage = strcat("Operator defined by E, Q and J is not skew-adjoint (Q'*E must be symmetric):\n norm(Q'*E - E'*Q)/norm(Q'*E) = ", num2str(norm_QE),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            % check skew symmetry of Q'*J*Q
            if (~issymmetric(varargin{3}'*varargin{1}*varargin{3},'skew'))
                % norm_QJQ = norm(Q'*J*Q + Q'*J'*Q)/norm(Q'*J*Q)
                norm_QJQ = norm(varargin{3}'*varargin{1}*varargin{3}+varargin{3}'*varargin{1}'*varargin{3})/norm(varargin{3}'*varargin{1}*varargin{3});
                if (norm_QJQ > opts.inputTolerance)
                    errorMessage = strcat("Operator defined by E, Q and J is not skew-adjoint (Q'*J*Q must be skew-symmetric):\n norm(Q'*J*Q+Q'*J'*Q)/norm(Q'*J*Q) = ", num2str(norm_QJQ),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            % check positive semi definiteness of Q'*E
            if (~phs.isPositiveDefiniteHermitian(varargin{3}'*varargin{5},opts.inputTolerance,'semi'))
                errorMessage = strcat("Q'*E is not positive semi-definite (Hamiltonian function must be bounded from below!)",errorMessage);
                error('phsmor:phs:wrongInput', errorMessage)
            end

            % W = [Q'*R*Q, Q'*P; P'*Q, S]
            W = [varargin{3}'*varargin{2}*varargin{3}, varargin{3}'*varargin{6}; varargin{6}'*varargin{3}, varargin{7}];
            
            % check symmetry of W
            if (~issymmetric(W))
                % W_transp = W' (--> avoid MATLAB numerical errors with W')
                W_transp = [varargin{3}'*varargin{2}'*varargin{3}, varargin{3}'*varargin{6}; varargin{6}'*varargin{3}, varargin{7}'];
                norm_W = norm(W-W_transp)/norm(W);
                if (norm_W > opts.inputTolerance)
                    errorMessage = strcat("W = [Q'*R*Q, Q'*P; P'*Q, S] is not symmetric!:\n norm(W-W')/norm(W) = ", num2str(norm_W),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage)
                end
            end

            % check positive semi definiteness of W
            if (~phs.isPositiveDefiniteHermitian(W,opts.inputTolerance, 'semi'))
                errorMessage = strcat("W = [Q'*R*Q, Q'*P; P'*Q, S] is not positive semi-definite",errorMessage);
                error('phsmor:phs:wrongInput', errorMessage)
            end

            isPH = true;
    % -------------------------------------------------------------------------------------
    else
    % -------------------------------------------------------------------------------------
    % check for system after creation
        sys = varargin{1};
        opts = sys.opts;
        % check if dimensions are correct
            sJ = size(sys.J);
            % J not square
            if (sJ(1) ~= sJ(2))
                error('phsmor:phs:wrongInput', 'dimensions of J are not correct');
            end
            sR = size(sys.R);
            % R not square or dimension different to J
            if ( (sR(1) ~= sR(2)) || ~isequal(sJ, sR) )     
                error('phsmor:phs:wrongInput', 'dimensions of R are not correct');
            end
            sQ = size(sys.Q);
            % Q not square or dimension different to J
            if ( (sQ(1) ~= sQ(2)) || ~isequal(sJ, sQ) )     
                error('phsmor:phs:wrongInput', 'dimensions of Q are not correct');
            end
            sE = size(sys.E);
            % E not square or dimension different to J
            if ( (sE(1) ~= sE(2)) || ~isequal(sJ, sE) )     
                error('phsmor:phs:wrongInput', 'dimensions of E are not correct');
            end
            sB = size(sys.B);
            % B has wrong state dimension
            if (sB(1) ~= sJ(1))             
                error('phsmor:phs:wrongInput', 'dimensions of B are not correct');
            end
            % different dimensions of P and B
            if (~isequal(sB, size(sys.P)))      
                error('phsmor:phs:wrongInput', 'dimensions of P are not correct');
            end
            % dimensions of S and N not equal
            if (~isequal(size(sys.N), size(sys.S))) 
                error('phsmor:phs:wrongInput', 'dimensions of S and N must be equal');
            end
            % dimensions of S and N not correct (compared to B)
            if (size(sys.N,1) ~= sB(2) || size(sys.N,2) ~= sB(2))   
                error('phsmor:phs:wrongInput', 'dimensions of S and N are not correct');
            end

        % check if E, J, R, Q, S and N define correct PH system
        % see [1], definition 5
            errorMessage = "\nConsider changing opts.inputTolerance to avoid this error when caused by numerical inaccuracy";

            if (~issymmetric(sys.S))
                norm_S = norm(sys.S-sys.S')/norm(sys.S);
                if (norm_S > opts.inputTolerance)
                    errorMessage = strcat("S is not symmetric:\n norm(S-S')/norm(S) = ", num2str(norm_S),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            if (~issymmetric(sys.N,'skew'))
                norm_N = norm(sys.N+sys.N')/norm(sys.N);
                if  (norm_N > opts.inputTolerance)
                    errorMessage = strcat("N is not skew-symmetric:\n norm(N+N')/norm(N) = ", num2str(norm_N),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            if (~issymmetric(sys.Q'*sys.E))
                norm_QE = norm(sys.Q'*sys.E - sys.E'*sys.Q)/norm(sys.Q'*sys.E);
                if (norm_QE > opts.inputTolerance)
                    errorMessage = strcat("Operator defined by E, Q and J is not skew-adjoint (Q'*E must be symmetric):\n norm(Q'*E - E'*Q)/norm(Q'*E) = ", num2str(norm_QE),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            if (~issymmetric(sys.Q'*sys.J*sys.Q,'skew'))
                norm_QJQ = norm(sys.Q'*sys.J*sys.Q+sys.Q'*sys.J'*sys.Q)/norm(sys.Q'*sys.J*sys.Q);
                if (norm_QJQ > opts.inputTolerance)
                    errorMessage = strcat("Operator defined by E, Q and J is not skew-adjoint (Q'*J*Q must be skew-symmetric):\n norm(Q'*J*Q+Q'*J'*Q)/norm(Q'*J*Q) = ", num2str(norm_QJQ),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage);
                end
            end

            if (~sys.isPositiveDefiniteHermitian(sys.Q'*sys.E,opts.inputTolerance,'semi'))
                errorMessage = strcat("Q'*E is not positive semi-definite (Hamiltonian function must be bounded from below!)",errorMessage);
                error('phsmor:phs:wrongInput', errorMessage)
            end

            W = [sys.Q'*sys.R*sys.Q, sys.Q'*sys.P; sys.P'*sys.Q, sys.S];
            if (~issymmetric(W))
                % W_transp = W' (--> avoid MATLAB numerical errors with W')
                W_tranps = [sys.Q'*sys.R'*sys.Q, sys.Q'*sys.P; sys.P'*sys.Q, sys.S'];
                norm_W = norm(W-W')/norm(W);
                if (norm_W > opts.inputTolerance)
                    errorMessage = strcat("W = [Q'*R*Q, Q'*P; P'*Q, S] is not symmetric!:\n norm(W-W')/norm(W) = ", num2str(norm_W),errorMessage);
                    error('phsmor:phs:wrongInput', errorMessage)
                end
            end

            if (~phs.isPositiveDefiniteHermitian(W,opts.inputTolerance, 'semi'))
                errorMessage = strcat("W = [Q'*R*Q, Q'*P; P'*Q, S] is not positive semi-definite",errorMessage);
                error('phsmor:phs:wrongInput', errorMessage)
            end

            isPH = true;
    % -------------------------------------------------------------------------------------
    end
end % of inputValidation()