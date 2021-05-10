%% Test Script
clc; clear; close all;

Ncells = 17;
n_ghost = 2;
xi = linspace(0,1,Ncells+1);
inputs = struct();
inputs.dt = 0.001;
inputs.time_range = [0,0.1];
inputs.order = 2;
inputs.kappa = 1;
inputs.exact_solution_type = 'initial_sine';

inputs.uLeft = -1;
inputs.uRight = 1;
soln = inviscid_burgers1D(xi,n_ghost,inputs);
flux = fluxes('scheme','godunov');
% limiter = limiters('scheme','van_leer');
limiter = limiters('scheme','beta_lim','beta',2);
RK = RK_Explicit('Method','Euler');
BC = @periodic_bc;
N = length(inputs.time_range(1):inputs.dt:inputs.time_range(2));

for i =1:N
    clf;
    soln.t = soln.t + soln.dt;
%     left = soln.U(soln.i_low-1:soln.i_high);
%     right = soln.U(soln.i_low:soln.i_high+1);
%     [left,right] = MUSCL_extrap(soln,limiter);
%     soln.F = flux.calc_flux(left,right);
%     soln.R = soln.residual(soln.F);
%     soln.U(soln.i) = soln.U(soln.i)-(1./soln.grid.dx(soln.i)).*soln.R*soln.dt;
    soln.U = RK.eval(@(t,u)burgers(t,u,flux,limiter,soln,BC),soln.U,soln.t,soln.dt);
    soln.E = soln.U(soln.i)-soln.calc_exact(soln.grid.xc(soln.i),soln.t);
%     hold on
%     plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,soln.t),'k')
%     plot(soln.grid.xc(soln.i),soln.U(soln.i),'r')
%     plot(soln.grid.xc(soln.i),soln.E,'r')
%     drawnow
%     hold off
end
%%
hold on;
% plot(soln.grid.xc(soln.i),soln.U(soln.i),'r')
% plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,soln.t),'k')
% xlim([soln.grid.xmin,soln.grid.xmax]);
hold off;
% plot(soln.E)

function dUdt = burgers(~,U,flux,limiter,soln,BC)
    dUdt = zeros(length(U),1);
%     [left,right] = MUSCL_extrap(soln,limiter,BC);
%     F = flux.calc_flux(left,right);
%     R = soln.residual(F);

    left1 = U(soln.i_low-1:soln.i_high);
    right1 = U(soln.i_low:soln.i_high+1);
    F1 = flux.calc_flux(left1,right1);
    R1 = soln.residual(F1);
    U1 = U;
    U1(soln.i) = U(soln.i)-(0.5./soln.grid.dx(soln.i)).*R1*soln.dt;
    [left,right] = MUSCL_extrap(soln,limiter,BC);
    left = left - left1 + U1(soln.i_low-1:soln.i_high);
    right = right - right1 + U1(soln.i_low:soln.i_high+1);
    F = flux.calc_flux(left,right);
    R = soln.residual(F);

    dUdt(soln.i) = -(1./soln.grid.dx(soln.i)).*R;
end