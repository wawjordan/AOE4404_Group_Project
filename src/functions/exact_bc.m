function U = exact_bc(soln,U)
soln.U(soln.i_high+1:end) = soln.Uex(soln.i_high+1:end);
soln.U(1:soln.i_low-1) = soln.Uex(1:soln.i_low-1);
end