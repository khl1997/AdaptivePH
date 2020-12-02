function [redSys,f,info] = adaptPH(sys,w,varargin)
%
% AdaptPH - calculates the order system redSys by application of the adaptive algorithm
%
% Syntax:
%   redSys = adaptPH(sys,w)
%   redSys = adaptPH(sys,w,opts)
%   redSys = adaptPH(sys,w,b)
%   redSys = adaptPH(sys,w,b,opts)
%   [redSys,f,info] = adaptPH(sys,w,...)
%
% Inputs Arguments:
% sys : phs object
% w   : vector containing the initial interpolation frequencies (all real number)
% b   : matrix containing the corresponding initial tangential direction (as columns)
% opts : struct containing options
%   
%   .maxiter        : max number of iteration allowed.
%   .tolz           : relative tolerance of the change of the 
%                        optimal frequencies between 2 consecutive iterations.
%   .tolf           : relative tolerance of the change of the computed 
%                        L-infinity norms between 2 consecutive iterations.
%   .interactive    :enables the interactive mode which allows to stop and
%                    plot the result of every iteration. (Default: false)
% 
% Outputs Arguments:
% redSys :  struct containing reduced port-Hamiltonian system.
%        redSys.J, redSys.R, redSys.B, redSys.Q, redSys.E: system matrix.
%
% f      :  the computed L-inf norm from the final iteration.
%
% info   : struct which contains the following information:
%    info.time        : computing time.
%    info.iterations  : number of iterations needed.
%    info.error       :  = 0: without error
%                        = 1: the maximum number of iteration exceeded. 
% 

    % input parsing  
    narginchk(2, 4);
    [b,opts] = parseInputs(sys,w,varargin{:});

    % start timing
    tic; 

    % prepare options for the reducePH
    optsReducePH = opts;
%   optsReducePH = rmfield(optsReducePH,{'maxiter','tolz','tolf'});

    % set parameters
    tolf = opts.tolf;
    tolz = opts.tolz;
    maxiter = opts.maxiter;
% ---------------------------------------------------------------------------------
    if ~sys.isMIMO    % SISO system
% ---------------------------------------------------------------------------------
        % determine the first reduced system with the initial interpolation points
        [redSys,V] = reducePH(sys,w,optsReducePH);

        % deterine the frequency of the peak L-inf norm
        [f,z] = norm(redSys,Inf);

        % print info
        redOrder = size(redSys.J,1);
        fprintf('Current iteration: 0\t Current interpolation frequency: %f\t Model order: %i\n',z,redOrder);

        % main loop
        iter = 0;
        z_old = -1e15;
        f_old = -1e15;
        while (((abs(z-z_old) > tolz*0.5*(abs(z) + abs(z_old))) && ...
            (abs(f-f_old) > tolf*0.5*(abs(f) + abs(f_old)) || ...
             abs(z-z_old) > 1e3*tolz*0.5*(abs(z)+abs(z_old))) && iter < maxiter) || opts.interactive)
            iter = iter + 1;
            z_old = z;
            f_old = f;
            w = [w z_old];
            % determine the new reduced system
            [redSys,V] = reducePH(sys,w,optsReducePH);

            % interactive mode
            if opts.interactive
               
            end

            % determine the frequency for next iteration
            [f,z] = norm(redSys,Inf);

            % print the info
            redOrder = size(redSys.J,1);
            fprintf('Current iteration: %i\t Current interpolation frequency: %f\t Model order: %i\n',iter,z,redOrder);
        end
% -------------------------------------------------------------------------------
    else     % MIMO system
% -------------------------------------------------------------------------------
        % determine the initial reduced system
        [redSys,V] = reducePH(sys,w,b,optsReducePH);

        % plot
        if opts.interactive
            figure();
            hold on
            bode(sys);
            set(findall(gcf,'type','line'),'linewidth',4)
            bode(redSys);
            disp('Press any key to continue')
            pause
        end

        % main loop
        
    end

    % Output
    info.time = toc;
    info.iterations = iter;
    if iter >= maxiter
        info.error = 1;
    else
        info.error = 0;
    end
end
%% Supporting function --------------------------------------------------------
function [b,opts] = parseInputs(sys,w,varargin)
% Parse the input to get information of tangential directions and options
    % check  the dimention of w
    if min(size(w)) ~= 1
        error('w must be a vector.');
    end
    if length(w) >= 0.5*sys.dim
        error('reduced order exceeds the original order.');
    end

    % opts: the options
    if ~isempty(varargin) && isstruct(varargin{end})
        optsProvided = true;
        % check existing fields
        if ~isfield(opts,'tolz')
            opts.tolz = 1e-6;
        end
        if ~isfield(opts,'tolf')
            opts.tolf = 1e-6;
        end
        if ~isfield(opts,'maxiter')
            opts.maxiter = 30;
        end
        if ~isfield(opts,'interactive')
            opts.interactive = false;
        end
    else  % opts not provided
        optsProvided = false;
        opts.tolz = 1e-6;
        opts.tolf = 1e-6;
        opts.maxiter = 30;
        opts.interactive = false;
    end
    % b: the tangential directions
    if (optsProvided && nargin == 3) || (~optsProvided && nargin == 2)
        % b not provided
        % initialize b as a matrix
        if sys.isMIMO % MIMO system
            % initialize b as ones matrix
            redOrder = 0;
            for i = 1:length(w)
                if w(i) == 0
                    redOrder = redOrder + 1;
                else
                    redOrder = redOrder + 2;
                end
            end
            b = zeros(size(sys.B,2),redOrder);
        else % SISO system
            b = []; % no need for b
        end
    elseif (optsProvided && nargin == 4) || (~optsProvided && nargin == 3)
        % b provided
        b = varargin{1};
    else

    end
end