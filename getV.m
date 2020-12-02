function V = getV(A,B,E,w,opts,flag_SIMI,b) % now only for SISO
% To determine Krylov subspace matrix
% Here we only consider the situation the reduced system interpolates
%      the original system at the frequencies.
% No consideration about the derivatives.
% Arguments:
% Inputs:
%    A,B,E     : corresponding matrix
%    w         : vector containing the interpolation frequencies, all real number
%    opts      : options
%       .orth  : 'none' : without orthogonalization
%                'mgs'  : single use of Gram-Schmidt
%                '2mgs' : double use of Gram-Schmidt
%
%    flag_SIMI : 'SISO' : SISO system
%                'MIMO' : MIMO system
%    b         : matrix containing tangential directions as columns
% Outputs:
%    V         : Krylov subspace matrix

narginchk(6, 7); % b is optional

% preallocate for V
redOrder = 0;
for i = 1:length(w)
    if w(i) == 0
        redOrder = redOrder + 1;
    else
        redOrder = redOrder + 2;
    end
end
V = zeros(size(A,1),redOrder);

indexV = 0;
% calculate Krylov subspace matrix
for i = 1:length(w)
    if w(i) == 0 % real shift
        if strcmp(flag_SIMI,'SISO') % SISO system
            v = -A\B;
        elseif strcmp(flag_SIMI,'MIMO') % MIMO system
            v = -A\B*b(:,i);
        else
            error('Wrong flag');
        end
        % orthogonalization
        switch opts.orth
        case 'none' % without orthogonalization
        case 'mgs' % single use
            v = gramSchmidt(v,V(:,1:indexV));
        case '2mgs' % double use
            v = gramSchmidt(v,V(:,1:indexV));
            v = gramSchmidt(v,V(:,1:indexV));
        end
        % extend V
        V(:,indexV+1) = v;
            indexV = indexV + 1;
    else % imaginary shift
        if strcmp(flag_SIMI,'SISO') % SISO system
            v = (1i*w(i)*E-A)\B;
        elseif strcmp(flag_SIMI,'MIMO') % MIMO system
            v = (1i*w(i)*E-A)\B*b(:,i);
        else
            error('Wrong flag');
        end
        % real part
        switch opts.orth
        case 'none' % without orthogonalization
            v_real = real(v);
        case 'mgs'
            v_real = gramSchmidt(real(v),V(:,1:indexV));
        case '2mgs'
            v_real = gramSchmidt(real(v),V(:,1:indexV));
            v_real = gramSchmidt(v_real,V(:,1:indexV));
        end
        % extend real part
        V(:,indexV+1) = v_real;
        indexV = indexV + 1;

        % imaginary part
        switch opts.orth
        case 'none' % without orthogonalization
            v_imag = imag(v);
        case 'mgs'
            v_imag = gramSchmidt(imag(v),V(:,1:indexV));
        case '2mgs'
            v_imag = gramSchmidt(imag(v),V(:,1:indexV));
            v_imag = gramSchmidt(v_imag,V(:,1:indexV));
        end
        % extend the imaginary part
        V(:,indexV+1) = v_imag;
        indexV = indexV + 1;
    end
end

function v = gramSchmidt(v,V)
% applies Gram Schmidt method to make v orthogonal to
% all columns in V
% copied from the arnoldiPH.m file
% simplified to single direction:  size(v,2) == 1 
    if ~isempty(V)
        % in case v is used as the first column of the matrix V
        for i = 1:size(V,2)
            v = v - (v'*V(:,i)/norm(V(:,i)))*V(:,i);
        end
        v = v/norm(v);
    else
        v = v/norm(v);
    end
end
