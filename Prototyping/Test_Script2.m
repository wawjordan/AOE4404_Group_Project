%% Test Script
clc; clear; close all;

Ncells = 16;
n_ghost = 2;
xi = linspace(-1,1,Ncells+1);
inputs = struct();
inputs.dt = 0.01;
inputs.time_range = [0,0.5];
inputs.order = 1;
inputs.kappa = -1;
inputs.exact_solution_type = 'initial_disc';

inputs.uLeft = -1;
inputs.uRight = 1;
soln = inviscid_burgers1D(xi,n_ghost,inputs);
flux = fluxes('scheme','roe');
limiter = limiters('scheme','van_leer');
% limiter = limiters('scheme','beta_lim','beta',1);
RK = RK_Explicit('Method','RK21');

N = length(inputs.time_range(1):inputs.dt:inputs.time_range(2));

hfig=figure(1);
clf(hfig);
dim = [7.5 5.5 3.25 2.5];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');

hold on;
plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,soln.t),'k--')
for i =1:N-1
    soln.t = soln.t + soln.dt;
    soln.U = RK.eval(@(t,u)burgers(t,u,flux,limiter,soln),soln.U,soln.t,soln.dt);
    soln.E = soln.U(soln.i)-soln.calc_exact(soln.grid.xc(soln.i),soln.t);
end
plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,soln.t),'k')
plot(soln.grid.xc(soln.i),soln.U(soln.i),'r')
title(sprintf('$N_{cells} = %d$, $t = %0.4f$',Ncells,soln.t),'interpreter','latex')
xlim([soln.grid.xmin,soln.grid.xmax]);
% ylim([min(soln.U)-1,max(soln.U)+1])
xlabel('x','interpreter','latex')
ylabel('u','rotation',0,'interpreter','latex')
hold off;
plen1= legend({'initial condition','exact','solution'},'location','northwest','interpreter','latex');
plen1.ItemTokenSize = [20,10];
filename = sprintf('order1_flux-%s_int-%s_N-%d',flux.scheme,RK.Method,Ncells);

% filename = sprintf('flux=%s_int=%s_dt=%0.4f_N=%d',flux.scheme,RK.Method,soln.dt,Ncells);
disp(filename);
print(['C:\Users\Will\Desktop\',filename],'-dpng','-r600');
function dUdt = burgers(~,U,flux,limiter,soln)
    dUdt = zeros(length(U),1);
    [left,right] = MUSCL_extrap(soln,limiter);
    F = flux.calc_flux(left,right);
    R = soln.residual(F);
    dUdt(soln.i) = -(1./soln.grid.dx(soln.i)).*R;
end