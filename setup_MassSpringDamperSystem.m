function sys = setup_MassSpringDamperSystem(n,m,k,d)
% Syntax:   sys = setup_MassSpringDamperSystem(n,m,k,d)
%
%   This demo provides a port-Hamiltonian system derived by the mechanical
%   example system of connected masses, springs and dampers (See picture
%   '2012_Polyuga_Effort-and-flow-constraint-reduction-methods_figure_5.png' [2]). 
%   Some parameters which can easily be adapted are:
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
%   Author: Julius Durmann
%   E-Mail: julius.durmann@tum.de
%   Date:   2020/09/12
%

%% if dampers are connected to a mass on one side and to a wall on the other side
    n = n/2;
% create matrices of system (PH-representation)
    J = diag(ones(2*n-1,1), 1) + diag(-1*ones(2*n-1,1), -1);

    diagR = zeros(2*n,1);
    indexR = 1:2:(2*n-1);
    diagR(indexR) = d;
    R = diag(diagR);

    ma = 1:2*(2*n+1):(2*n*2*n);     %for Q
    ks = 2*n+2:2*(2*n+1):(2*n*2*n); %for Q
    Q = zeros(2*n);
    Q(ma) = 1./m;
    Q(ks) = k;

    B = zeros(2*n,1); B(1) = 1;

    % create phs system
    opts.inputValidation = true;
    sys = phs(J,R,Q,B, opts);

%% if dampers are connected to masses on both sides: 
%     if mod(n,2)~= 0
%         error('The system dimension must be even!');
%     end
%     
%     %% Q
%     index_m = 1:2:n-1;
%     index_s = 2:2:n;
%     diagQ = zeros(n,1);
%     diagQ(index_m) = 1./m;
%     diagQ(index_s) = k;
%     Q = diag(diagQ);
% 
%     %% J
%     J = diag(ones(n-1,1),1) + diag(-ones(n-1,1),-1);
% 
%     %% R
% 
%     if length(d) < 2
%         d = repmat(d,n/2,1);
%     end
%     index_upperLowerDiag = 1:2:n-2;
%     upperDiag = zeros(n-2,1);
%     lowerDiag = zeros(n-2,1);
%     diagonal = zeros(n,1);
% 
%     upperDiag(index_upperLowerDiag) = -d(1:end-1);
%     lowerDiag(index_upperLowerDiag) = -d(2:end);
% 
%     index_diag = 1:2:n-1;
%     diagonal(index_diag) = d(1:end) + [0; d(1:end-1)];
% 
%     R = diag(diagonal) + diag(upperDiag,2) + diag(lowerDiag, -2);
% 
%     %% B
%     B = zeros(n,1);
%     B(1) = 1;
% 
%     %% system
%     sys = phs(J,R,Q,B);

end