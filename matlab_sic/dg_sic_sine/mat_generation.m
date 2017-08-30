function A = mat_generation(x, st, n, k, dim)
% function to generate nonlinear matrix A
% st: start of y; n: y[0]~y[n]; k: half of channel; dim: nonlinear
% dimension

A0 = toeplitz(x(st + k - 1:st + n + k - 1),fliplr(x(st - k :st + k -1)));
% A0 = fliplr(A0);            % base for linear cancellation
A = [];
for i = 1:dim
    A = [A, A0.^i];
end

end
    