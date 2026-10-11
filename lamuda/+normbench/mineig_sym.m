function value = mineig_sym(matrix)
%MINEIG_SYM Smallest eigenvalue of the symmetric part.
matrix = 0.5*(matrix+matrix.');
value = min(real(eig(full(matrix))));
end
