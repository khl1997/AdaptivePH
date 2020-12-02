function v = funplot(sys,wmin,wmax,n)
% Input Arguments:
%   sys   : phs-object
%   wmin  : minimum interpolation frequency
%   wmax  : maximum interpolation frequency
%   n     : number of evaluation points
w = linspace(wmin,wmax,n);
% preallocate for v
v = zeros(length(w),1);
for i = 1:n
    H = ((sys.B+sys.P)'*sys.Q) * ((1i*w(i)*sys.E-(sys.J-sys.R)*sys.Q) \ (sys.B-sys.P));
    v(i) = max(svd(H));
end

% plot the max singular values over interval [wmin,wmax]
figure();
plot(w,v);
xlabel('Frequency');
ylabel('Maximun sigular value');

return