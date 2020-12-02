function sys = setup_MassSpringDamperMIMOSystem(nsys,n,m,k,d)
%   This demo provides a port-Hamiltonian MIMO system derived by the mechanical
%   example system of connected masses, springs and dampers (See picture
%   '2012_Polyuga_Effort-and-flow-constraint-reduction-methods_figure_5.png' [2]). 
%   Here, nsys such systems are are aligned parallel, each having its own
%   input and output.
%
%   Some parameters which can easily be adapted are:
%       nSys - number of parallel systems
%       n - This changes the dimension of the system
%       m, k, d -   These parameters change the mass, spring constant and
%                   damping coefficient respecitvely
%
%   References:
%       [1] Gugercin et al. (2009), in Proceedings of the 48h IEEE 15.12.2009 - 18.12.2009
%               "Interpolation-based H2 model reduction for port-Hamiltonian systems"
%       [2] Polyuga et al. (2012), 
%               "Effort- and flow-constraint reduction methods for
%               structure preserving model reduction of port-Hamiltonian
%               systems"
%
% 
% 
%   Author: Julius Durmann
%   E-Mail: julius.durmann@tum.de
%   Date:   2020/09/12
%

    %% create linear mass-spring-damper-system
    % create matrices of one system (PH-representation)
    sys_single = setup_MassSpringDamperSystem(n,m,k,d);
    
    c = cell(1,nsys);
    for i = 1:nsys
        c{i} = sys_single.J;
    end
    J = diagonalMatrix(c{:});
    for i = 1:nsys
        c{i} = sys_single.R;
    end
    R = diagonalMatrix(c{:});
    for i = 1:nsys
        c{i} = sys_single.Q;
    end
    Q = diagonalMatrix(c{:});
    
    B = zeros(size(sys_single.B,1)*nsys, size(sys_single.B,2)*nsys);
    for i = 1:nsys
        B((i-1)*size(sys_single.B,1)+1:i*size(sys_single.B,1),(i-1)*size(sys_single.B,2)+1:i*size(sys_single.B,2)) = sys_single.B;
    end

    % create phs system
    opts.inputValidation = true;
    sys = phs(J,R,Q,B, opts);

end

function matrix = diagonalMatrix(varargin)
    dim = [0,0];
    for i = 1:length(varargin)
        dim = dim + size(varargin{i});
    end
    matrix = zeros(dim);
    
    idx1 = 1;
    idx2 = 1;
    for i = 1:length(varargin)
        matrix(idx1:idx1+size(varargin{i},1)-1, idx2:idx2+size(varargin{i},2)-1) = varargin{i};
        idx1 = idx1 + size(varargin{i},1);
        idx2 = idx2 + size(varargin{i},2);
    end
end