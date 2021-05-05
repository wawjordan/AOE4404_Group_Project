function U = periodic_bc(soln,U)
i_low = soln.i_low:soln.i_low+soln.grid.n_ghost-1;
i_high = soln.i_high-soln.grid.n_ghost+1:soln.i_high;
U(soln.i_high+1:end) = U(i_low);
U(1:soln.i_low-1) = U(i_high);
end