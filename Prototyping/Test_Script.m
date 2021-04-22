%% Test Script
clc; clear; close all;

Ncells = 200;
xi = linspace(-1,2,Ncells+1);
inputs = struct();
inputs.n_ghost = 2;
inputs.dt = 0.01;
inputs.time_range = [0,1.5];
soln = inviscid_burgers1D(xi,inputs);
flux = fluxes('scheme','eo');
RK = RK_Explicit('Method','Euler');

hold on;
N = length(inputs.time_range(1):inputs.dt:inputs.time_range(2));
for i =1:N
    clf;
    soln.t = soln.t + soln.dt;
    left = soln.U(soln.i_low-1:soln.i_high);
    right = soln.U(soln.i_low:soln.i_high+1);
    soln.F = flux.calc_flux(left,right);
    soln.R = soln.residual(soln.F);
%     soln.U(soln.i) = soln.U(soln.i)-(1./soln.grid.dx(soln.i)).*soln.R*soln.dt;
    soln.U = RK.eval(@(t,u)burgers(t,u,flux,soln),soln.U,soln.t,soln.dt);
    hold on;
    plot(soln.grid.xc(soln.i),soln.U(soln.i),'r')
    plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,soln.t),'k')
    
    hold off;
    pause(0.01);
end


function dUdt = burgers(~,U,flux,soln)
    dUdt = zeros(length(U),1);
    left = U(soln.i_low-1:soln.i_high);
    right = U(soln.i_low:soln.i_high+1);
    F = flux.calc_flux(left,right);
    R = soln.residual(F);
    dUdt(soln.i) = -(1./soln.grid.dx(soln.i)).*R;
end