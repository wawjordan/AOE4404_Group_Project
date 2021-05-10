function [soln,out] = solve_burgers(xi,inputs,flux,limiter,RK,bc,out_interval)
out = struct();
soln = inviscid_burgers1D(xi,2,inputs);
N = length(inputs.time_range(1):inputs.dt:inputs.time_range(2));
out.t = nan(N,1);
out.tplot = nan(N,1);
out.norm1 = nan(N,1);
out.norm2 = nan(N,1);
out.norminf = nan(N,1);
out.U = cell(N,1);
out.Uex = cell(N,1);
out.E = cell(N,1);


out.t(1) = soln.t;
out.tplot(1) = soln.t;
out.norm1(1) = 0;
out.norm2(1) = 0;
out.norminf(1) = 0;
out.U{1} = soln.U(soln.i);
out.Uex{1} = soln.Uex(soln.i);
out.E{1} = soln.U(soln.i)-soln.Uex(soln.i);

out.x = soln.grid.xc(soln.i);

for j =1:N-1
    soln.t = soln.t + soln.dt;
    soln.Uex = soln.calc_exact(soln.grid.xc,soln.t);
    soln.U = RK.eval(@(t,u)burgers(t,u,flux,limiter,soln,bc),...
        soln.U,soln.t,soln.dt);
    soln.E = soln.U(soln.i)-soln.Uex(soln.i);
    
    if mod(j,out_interval)==0
        out.t(j+1) = soln.t;
        out.norm1(j+1) = sum(abs(soln.E))/(length(xi)-1);
        out.norm2(j+1) = sqrt(sum((soln.E).^2)/(length(xi)-1));
        out.norminf(j+1) = max(abs(soln.E));
        out.tplot(j+1) = soln.t;
        out.U{j+1} = soln.U(soln.i);
        out.Uex{j+1} = soln.Uex(soln.i);
        out.E{j+1} = soln.E;
    end
    
end

out.t = out.t(~isnan(out.t));
out.tplot = out.tplot(~isnan(out.tplot));
out.norm1 = out.norm1(~isnan(out.norm1));
out.norm2 = out.norm2(~isnan(out.norm2));
out.norminf = out.norminf(~isnan(out.norminf));
out.U = out.U(~cellfun('isempty',out.U));
out.Uex = out.Uex(~cellfun('isempty',out.Uex));
out.E = out.E(~cellfun('isempty',out.E));

function dUdt = burgers(~,U,flux,limiter,soln,bc)
    dUdt = zeros(length(U),1);
%     [left,right] = MUSCL_extrap(soln,limiter,bc);
    left1 = U(soln.i_low-1:soln.i_high);
    right1 = U(soln.i_low:soln.i_high+1);
    F1 = flux.calc_flux(left1,right1);
    R1 = soln.residual(F1);
    U1 = U;
    U1(soln.i) = U(soln.i)-(0.5./soln.grid.dx(soln.i)).*R1*soln.dt;
    [left,right] = MUSCL_extrap(soln,limiter,bc);
    left = left - left1 + U1(soln.i_low-1:soln.i_high);
    right = right - right1 + U1(soln.i_low:soln.i_high+1);
% 
    F = flux.calc_flux(left,right);
    R = soln.residual(F);
    dUdt(soln.i) = -(1./soln.grid.dx(soln.i)).*R;
end

end

