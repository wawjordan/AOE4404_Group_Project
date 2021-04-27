%% Test Script
clc; clear; close all;

Ncells = 16;
n_ghost = 2;
xi = linspace(-1,1,Ncells+1);
inputs = struct();
inputs.dt = 0.01;
inputs.time_range = [0,0.5];
inputs.order = 2;
inputs.kappa = 1;
inputs.exact_solution_type = 'initial_disc';

inputs.uLeft = -1;
inputs.uRight = 1;
soln = inviscid_burgers1D(xi,n_ghost,inputs);
flux = fluxes('scheme','godunov');
limiter = limiters('scheme','van_leer');
% limiter = limiters('scheme','beta_lim','beta',1);
RK = RK_Explicit('Method','RK41');

N = length(inputs.time_range(1):inputs.dt:inputs.time_range(2));

for i =1:N
    soln.t = soln.t + soln.dt;
%     left = soln.U(soln.i_low-1:soln.i_high);
%     right = soln.U(soln.i_low:soln.i_high+1);
%     [left,right] = MUSCL_extrap(soln,limiter);
%     soln.F = flux.calc_flux(left,right);
%     soln.R = soln.residual(soln.F);
%     soln.U(soln.i) = soln.U(soln.i)-(1./soln.grid.dx(soln.i)).*soln.R*soln.dt;
    soln.U = RK.eval(@(t,u)burgers(t,u,flux,limiter,soln),soln.U,soln.t,soln.dt);
    soln.E = soln.U(soln.i)-soln.calc_exact(soln.grid.xc(soln.i),soln.t);
end
hold on;
plot(soln.grid.xc(soln.i),soln.U(soln.i),'r')
plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,soln.t),'k')
xlim([soln.grid.xmin,soln.grid.xmax]);
hold off;


function dUdt = burgers(~,U,flux,limiter,soln)
    dUdt = zeros(length(U),1);
    [left,right] = MUSCL_extrap(soln,limiter);
    F = flux.calc_flux(left,right);
    R = soln.residual(F);
    dUdt(soln.i) = -(1./soln.grid.dx(soln.i)).*R;
end