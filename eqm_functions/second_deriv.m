function out = second_deriv(etaout, f)
% function out = second_deriv(etaout, f)
%
% Calculate the second derivative using centered differences in the
% interior, forward difference at the left boundary point, and backward
% difference at the right boundary point.
%
% The order and accuracy are both 2. The right boundary point is a
% reflecting boundary, so the optimal handling of the right boundary point
% is to apply a ghost node approach and calculate the second derivative
% consistent with a zero first derivative at the right boundary point.
    N = length(etaout);
	Dxx_coefs  = cell2mat(arrayfun(@(x0) get_fd_coef(2, 2, x0, etaout, 1), etaout(2:end-1), 'UniformOutput', false));
	Dxx_coefs0 = get_fd_coef(2, 2, etaout(1), etaout, 2);
	Dxx_coefsN = get_fd_coef(2, 2, etaout(end), etaout, 2);
	diag_mat1 = [Dxx_coefs(2:end,1); 0; 0; 0];
	diag_mat2 = [0; Dxx_coefs(2:end,2); 0; 0];
	diag_mat3 = [0; 0; Dxx_coefs(2:end,3); 0];
	diag_mat4 = [0; 0; 0; Dxx_coefs(2:end,4)];
	diag_mat = [diag_mat1 diag_mat2 diag_mat3 diag_mat4];
	D_xx = spdiags(diag_mat, -2:1, N, N);
	D_xx(1, 1:length(Dxx_coefs0)) = Dxx_coefs0;
	D_xx(end, end - length(Dxx_coefsN) + 1:end) = Dxx_coefsN;
    out = D_xx * f;
end
