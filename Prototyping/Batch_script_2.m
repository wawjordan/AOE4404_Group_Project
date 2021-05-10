%% Test Script
clc; clear; close all;
fluxess = {'godunov'};
limiterss = {'van_leer','van_albada','minmod','beta_lim'};
RK = RK_Explicit('Method','Euler');
BC = @periodic_bc;

Ncells = 2.^(4:10);
M = length(Ncells);
dt = 0.01*2.^(0:-1:-(M-1));


K = length(limiterss);

inputs = struct();
inputs.time_range = [0,0.5];
inputs.order = 2;
inputs.kappa = 1;
input.beta = 1.5;
inputs.exact_solution_type = 'initial_sine';

labels = cell(length(limiterss),1);

filename = sprintf(...
    'test=%s_flux=%s_k=%0.1f_BC=%s_tf=%0.2f',...
    strrep(inputs.exact_solution_type,'initial_',''),...
    fluxess{1},...
    inputs.kappa,...
    func2str(BC),...
    inputs.time_range(2));

dirname1 = ['C:\Users\Will Jordan\Documents\MATLAB\Grad_School\AOE4404',...
    '\Project\Results\Data\'];
dirname2 = ['C:\Users\Will Jordan\Documents\MATLAB\Grad_School\AOE4404',...
    '\Project\Results\Figures\'];

E = struct();
intervals = (inputs.time_range(2)-inputs.time_range(1))./dt;
out_interval = 1;
for i = 1:10
    out_interval = max(out_interval,gcd(i,intervals(1)));
end
intervals = intervals/out_interval;
for j = 1:K
for i = 1:M
flux = fluxes('scheme',fluxess{1});
if strcmp(limiterss{j},'beta_lim')
limiter = limiters('scheme',limiterss{j},'beta',input.beta);
else
limiter = limiters('scheme',limiterss{j});
end
fprintf('Limiter=%s    N = %d     (%d/%d)\n',...
    limiter.label,Ncells(i),(j-1)*M+i,M*K);
inputs.dt = dt(i);
xi = linspace(0,1,Ncells(i)+1);
[soln,out] = solve_burgers(xi,inputs,flux,limiter,RK,BC,intervals(i));
E(i,j).t = out.t;
E(i,j).x = out.x;
E(i,j).dx = min(soln.grid.dx);
E(i,j).dt = dt(i);
E(i,j).flux = flux.scheme;
E(i,j).limiter = limiter.scheme;
E(i,j).tplot = out.tplot;
E(i,j).error = out.E;
E(i,j).norm1 = out.norm1;
E(i,j).norm2 = out.norm2;
E(i,j).norminf = out.norminf;
E(i,j).soln = out.U;
E(i,j).exsoln = out.Uex;
end
labels{j} = limiter.label;
end
%% Plotting

save([dirname1,filename,'.mat'],'E');

% hfig1 = ooa_plot(1,E,K,M,Ncells,labels,1);
% exportgraphics(gcf,[dirname2,filename,'N1','.png'],'Resolution',600)
% 
% hfig2 = ooa_plot(2,E,K,M,Ncells,labels,2);
% exportgraphics(gcf,[dirname2,filename,'N2','.png'],'Resolution',600)
% 
% hfig3 = ooa_plot(3,E,K,M,Ncells,labels,3);
% exportgraphics(gcf,[dirname2,filename,'N3','.png'],'Resolution',600)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function hax = stdsubplot(i,j,k)
hax = subplot(i,j,k);
set(hax,'Units','Inches');
set(hax,'DefaultAxesFontName','Helvetica');
set(hax,'DefaultTextFontName','Helvetica'); 
set(hax,'DefaultAxesFontSize',8);
set(hax,'DefaultTextFontSize',8);
end

function hfig = stdplot(i)
hfig=figure(i);
clf(hfig);
dim = [7.5 5.5 6.25 2.5];
set(hfig,'Units','Inches','Position',dim);
set(hfig,'DefaultAxesFontName','Helvetica');
set(hfig,'DefaultTextFontName','Helvetica'); 
set(hfig,'DefaultAxesFontSize',8);
set(hfig,'DefaultTextFontSize',8);
set(hfig,'PaperUnits',get(gcf,'Units'));
pos = get(hfig,'Position');
set(hfig,'PaperPosition',[0 0 pos(3) pos(4)]);
set(gca,'Units','Inches');
set(hfig,'DefaultLineLineWidth',0.5)
end

function hfig = ooa_plot(i,E,K,M,Ncells,labels,norm)

mark = cell(K,1);
marker = {'o-','^-',''};
color = {'k','r','g','b','m','c'};
for ii = 1:K
    mark{ii} = strcat(marker{1},color{ii});
end


norms = {'norm1','norm2','norminf'};
normlabel = {'1','2','\infty'};

hfig = stdplot(i);
ax1 = stdsubplot(1,2,1);
hold on;
for j=1:K
    e1 = [E(:,j).(sprintf('%s',norms{norm}))];
    plot([E(:,j).dx],e1(end,:),mark{j});
end
hold off;
set(ax1,'xscale','log')
set(ax1,'yscale','log')
xlabel('$\Delta{x}$','interpreter','latex')
ylabel(sprintf('$\\Vert{\\epsilon}\\Vert_{%s}$',normlabel{norm}),...
    'interpreter','latex','rotation',0)
ylim([1e-6,1e-1])
pbaspect([1.2 1 1])

ax2 = stdsubplot(1,2,2);
r = Ncells(2:end)./Ncells(1:end-1);
hold on;
for j=1:K
    e1 = [E(1:M-1,j).(sprintf('%s',norms{norm}))];
    e2 = [E(2:M,j).(sprintf('%s',norms{norm}))];
    eex = e1(end,:)./e2(end,:);
px = log(eex)./log(r);
plot([E(2:end,j).dx],px,mark{j})
end
hold off;
set(ax2,'xscale','log')
ylim([0,5])
xlim(get(ax1,'xlim'));
xlabel('$\Delta{x}$','interpreter','latex')
ylabel('$\hat{p}\quad\quad$','interpreter','latex','rotation',0)
legend(labels,'interpreter','latex','location','northeast');
pbaspect([1.2 1 1])
end