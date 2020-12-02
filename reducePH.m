function [redSys,V] = reducePH(sys,w,varargin)
% reduce the PH system with the interpolation frequency in w
%
% Syntax:
%  redSys = reducePH(sys,w)
%  redSys = reducePH(sys,w,opts)
%  redSys = reducePH(sys,w,b)
%  redSys = reducePH(sys,w,b,opts)
% 
% Arguments 
% Inputs:
%  - Required Arguments: 
%     sys    : phs object
%     w      : vector of interpolation frequencies (can also be a singular frequency)
%  - Optional Arguments:
%     b      : matrix of the tangential directions as columns
%     opts   : options
%        .orth    : Gram Schmidt method used while forming V (['none'/'mgs'/'2mgs'], Default:'2mgs')
%        .reorth  : use MATLAB function orth to reorthogonalize V (Default: true)
% Outputs:
%     redSys : phs object
%     V      : corresponding Krylov subspace matrix

    narginchk(2, 4);
    [w,b,opts] = parseInputs(w,varargin{:});

    % initialize system matrices
    J = sys.J;
    R = sys.R;
    Q = sys.Q;
    B = sys.B;
    P = sys.P;
    E = sys.E;

    % determine Krylov subspace matrix V
    if ~sys.isMIMO % SISO system
        V = getV((J-R)*Q,B-P,E,w,opts,'SISO');
    else
        V = getV((J-R)*Q,B-P,E,w,opts,'MIMO',b);
    end

    % reorthogonalize V
    if opts.reorth
        V = orth(V);
    end

    % create reduced model
    W = Q*V/(V'*Q*V);

    J_r = W'*J*W;
    R_r = W'*R*W;
    Q_r = V'*Q*V;

    if sys.isDescriptor
        E_r = W'*E*V;
        [J_r,R_r] = correctSymmetry(J_r,R_r);
    else
        E_r = eye(length(J_r));
        [J_r,R_r,Q_r] = correctSymmetry(J_r,R_r,Q_r);
    end
    B_r = W'*B;
    P_r = W'*P;

    % create phsRed object
    parameters = opts;
    parameters.shifts = w;
    if sys.isMIMO
        parameters.tangentDirections = b;
    end
    redSys = phsRed(J_r,R_r,Q_r,B_r,E_r,P_r,'1',parameters);
end

% Supporting function---------------------------------------------------------------

function [w,b,opts] = parseInputs(w,varargin)
    % opts
    if ~isempty(varargin) && isstruct(varargin{end})
        optsProvided = true;
        opts = varargin{end};
        if ~isfield(opts,'reorth')
            opts.reorth = true;
        end

        if ~isfield(opts,'orth')
            opts.orth = '2mgs';
        else
            switch opts.orth
            case 'none'
            case 'mgs'
            case '2mgs'
            otherwise
                opts.orth = '2mgs';
            end
        end
    else
        optsProvided = false;
        opts = struct();
        opts.orth = '2mgs';
        opts.reorth = true;
    end

    % w and b
    if min(size(w)) ~= 1
        error('w must be a vector.');
    end

    if optsProvided && nargin == 3
        b = varargin{1};
    elseif optsProvided && nargin == 2
        b = [];
    elseif ~optsProvided && nargin == 2
        b = varargin{1};
    else
        b = [];
    end
end