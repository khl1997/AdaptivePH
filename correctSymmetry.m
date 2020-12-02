function [J_out,R_out,Q_out] = correctSymmetry(J,R,Q)
% correct slightly the matrices which are possibly not symmetric
J_out = 0.5*((J-R)-(J'-R'));
R_out = -0.5*((J-R)+(J'-R'));
if nargin == 3
    Q_out = 0.5*(Q+Q');
end