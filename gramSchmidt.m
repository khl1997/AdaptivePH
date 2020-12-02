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