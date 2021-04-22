%% Test Script
clc; clear; close all;
% fprintf('Blah blah blah\n')
Ncells = 128;
xi = linspace(-1,2,Ncells+1);
inputs = struct();
inputs.n_ghost = 2;
inputs.dt = 0.1;
inputs.time_range = [0,1];
%inputs.IC = @(x) sin(pi*x);
soln = inviscid_burgers1D(xi,inputs);
hold on;
t = 0;
for i =1:10
plot(soln.grid.xc,soln.calc_exact(soln.grid.xc,t),'.')
t = t + soln.dt;
end
hold off;
xlim([-1,2])
ylim([-1,2])


