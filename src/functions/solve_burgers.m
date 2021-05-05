function soln = solve_burgers(xi,inputs,flux,limiter,RK,bc)
soln = inviscid_burgers1D(xi,2,inputs);
N = length(inputs.time_range(1):inputs.dt:inputs.time_range(2));
for j =1:N-1
    soln.t = soln.t + soln.dt;
    soln.Uex = soln.calc_exact(soln.grid.xc,soln.t);
    soln.U = RK.eval(@(t,u)burgers(t,u,flux,limiter,soln,bc),soln.U,soln.t,soln.dt);
    soln.E = soln.U(soln.i)-soln.Uex;
end

function dUdt = burgers(~,U,flux,limiter,soln,bc)
    dUdt = zeros(length(U),1);
    [left,right] = MUSCL_extrap(soln,limiter,bc);
    F = flux.calc_flux(left,right);
    R = soln.residual(F);
    dUdt(soln.i) = -(1./soln.grid.dx(soln.i)).*R;
end

end

