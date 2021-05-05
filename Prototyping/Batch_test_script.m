%% Test Script
clc; clear; close all;

% fluxess = {'godunov','eo','roe'};
fluxess = {'roe'};
limiterss = {'van_leer','van_albada','minmod','beta_lim'};
RK = RK_Explicit('Method','RK41');
BC = @periodic_bc;

Ncells = 2.^(4:10);
M = length(Ncells);
dt = 0.001*2.^(0:-1:-(M-1));


K = length(limiterss);
% K = length(fluxess);

mark = cell(K,1);
marker = {'o-','^-',''};
color = {'k','r','g','b','m','c'};
for ii = 1:K
    mark{ii} = strcat(marker{1},color{ii});
end

inputs = struct();
inputs.time_range = [0,0.1];
inputs.order = 2;
inputs.kappa = -1.0;
input.beta = 1.5;
inputs.exact_solution_type = 'initial_sine';
% inputs.uLeft = -1;
% inputs.uRight = 1;

labels = {'van Leer limiter','van Albada limiter','minmod limiter',...
    sprintf('$\\beta$-limiter ($\\beta=%0.2f$)',input.beta)};
% labels = {'Engquist-Osher Flux','Roe Flux'};
% labels = {'Euler''s method','Heun''s method','RK41'};


E = struct();
for j = 1:K
for i = 1:M
fprintf('Limiter=%s    N = %d     (%d/%d)\n',labels{j},Ncells(i),(j-1)*M+i,M*K);
flux = fluxes('scheme',fluxess{1});
% limiter = limiters('scheme',limiterss{1});
if strcmp(limiterss{j},'beta_lim')
limiter = limiters('scheme',limiterss{j},'beta',input.beta);
else
limiter = limiters('scheme',limiterss{j});
end
inputs.dt = dt(i);
xi = linspace(0,1,Ncells(i)+1);
% xi = linspace(-1,2,Ncells(i)+1);
soln = solve_burgers(xi,inputs,flux,limiter,RK,BC);
E(i,j).x = soln.grid.xc(soln.i);
E(i,j).dx = min(soln.grid.dx);
E(i,j).error = soln.E;
E(i,j).norm1 = sum(abs(soln.E))/Ncells(i);
E(i,j).norm2 = sqrt(sum((soln.E).^2)/Ncells(i));
E(i,j).norminf = max(abs(soln.E));
E(i,j).soln = soln.U(soln.i);
E(i,j).exsoln = soln.Uex;
end
end
%% Plotting
hfig=figure(1);
clf(hfig);
dim = [7.5 5.5 6.25 2.5];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(gca,'Units','Inches');

ax1 = subplot(1,2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax1,'Units','Inches');
set(ax1,'DefaultAxesFontName','Helvetica');
set(ax1,'DefaultTextFontName','Helvetica'); 
set(ax1,'DefaultAxesFontSize',8);
set(ax1,'DefaultTextFontSize',8);

hold on;
for j=1:K
plot([E(:,j).dx],[E(:,j).norm1],mark{j});
end
hold off;
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$\Delta{x}$','interpreter','latex')
ylabel('$\Vert{\epsilon}\Vert_1$','interpreter','latex','rotation',0)
ylim([1e-4,1e-1])
pbaspect([1.2 1 1])


ax2 = subplot(1,2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(ax2,'Units','Inches');
set(ax2,'DefaultAxesFontName','Helvetica');
set(ax2,'DefaultTextFontName','Helvetica'); 
set(ax2,'DefaultAxesFontSize',8);
set(ax2,'DefaultTextFontSize',8);

r = Ncells(2:end)./Ncells(1:end-1);
hold on;
for j=1:K
eex = [E(1:M-1,j).norm1]./[E(2:M,j).norm1];
px = log(eex)./log(r);
plot([E(2:end,j).dx],px,mark{j})
end
hold off;
set(gca,'xscale','log')
ylim([0,3])
xlim(get(ax1,'xlim'));
xlabel('$\Delta{x}$','interpreter','latex')
ylabel('$\hat{p}\quad\quad$','interpreter','latex','rotation',0)
legend(labels,'interpreter','latex','location','northeast');
pbaspect([1.2 1 1])

dirname1 = 'C:\Users\Will Jordan\Documents\MATLAB\Grad_School\AOE4404\Project\Results\Data\';
dirname2 = 'C:\Users\Will Jordan\Documents\MATLAB\Grad_School\AOE4404\Project\Results\Figures\';
% filename = sprintf('sine_test_int=%s',RK.Method);
filename = sprintf('sine_test_int=%s',RK.Method);
% filename = 'flux=EO_limiter=vary_int=Heun';
filename = [dirname1,filename];
save([dirname1,filename],E);
exportgraphics(gcf,[dirname2,filename,'.png'],'Resolution',600)
% print(filename,'-dpng','-r600');
% crop(strcat(filename,'.png'))