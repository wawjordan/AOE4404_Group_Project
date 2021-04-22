%% Test Script
clc; clear; close all;
% fprintf('Blah blah blah\n')
Ncells = 500;
Lmax = 3;
L = 1;
A = 1;
u0 =0;
xi = linspace(0,Lmax,Ncells-1);
t = 0:0.1/pi:10/pi;
g = @(x) A*sin(pi*x/L).*(x<L) + u0;
gprime = @(x) A*pi/L*cos(pi*x/L).*(x<L);

% u = @(x,t,x0) g(x-g(x0)*t);
u1 = A;
u2 = u0;
x01 = 0.5;
x02 = 1;
ts = 1/pi;

% x = xi + g(xi)*ts;
% u = g(x-g(xi)*ts);
% plot(x,u,'k')
hold on;
for i = 1:length(t)
    clf;
    x = xi + g(xi)*t(i);
    x1 = x01 + g(x01)*t(i);
    x2 = x02 + g(x02)*t(i);
    u1 = g(x1-g(x01)*t(i));
    u2 = g(x2-g(x02)*t(i));
    C = 0.5*(u2+u1);
    u = g(x-g(xi)*t(i));
    if t(i) >= ts
        xs = 0.5*(x01+x02)+C*t(i);
        for j = 1:length(u)
            if x(j) > xs
                u(j) = u0;
%         ind = find(x>xs,1);
%         u(ind:end) = u0;
%         x(x>xs) = xs;
            end
            if x(j) < xi(j)
                x(j) = xs;
            end
        end
    end
    plot(x,u,'k')
    xlim([0,3])
    ylim([0,1]);
    drawnow;
%     pause(0.1)
end
hold off;

% for i = 1:length(t)
%     clf;
%     u = zeros(length(xi),1);
%     x = xi;% + g(xi)*t(i);
%     for j = 1:length(xi)
%         if t(i)<=1
%             if x(j)<t(i)
%                 u(j) = 1;
%             elseif (x(j)>t(i))&&(x(j)<1)
%                 u(j) = (1-x(j))/(1-t(i));
%             else
%                 u(j) = 0;
%             end
%         else
%             if x(j)<(t(i)+1)/2
%                 u(j) = 1;
%             elseif x(j)>(t(i)+1)/2
%                 u(j) = 0;
%             end
%         end
%     end
%         
%     
%     plot(x,u,'k')
%     xlim([0,3])
%     ylim([-1,2]);
%     drawnow;
% %     pause(0.1)
% end
% hold off;


% inputs = struct();
% inputs.n_ghost = 2;
% inputs.dt = 0.1;
% inputs.time_range = [0,10];
% inputs.IC = @(x) sin(pi*x);
% soln = inviscid_burgers1D(xi,inputs);
% plot(soln.grid.xc,soln.U)