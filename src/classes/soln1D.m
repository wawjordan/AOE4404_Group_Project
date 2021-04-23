classdef soln1D < handle
   properties
      grid, IC
      N {mustBeNumeric}
      dt {mustBeNumeric}
      U, F, E, R {mustBeNumeric}
      Rnorm, Rinit, count {mustBeNumeric}
   end
   methods
       function this = soln1D(xi,n_ghost)
           this.N = length(xi)-1;
           p = inputParser;
           addRequired(p,'xi',...
               @(x)validateattributes(x,{'numeric'},...
               {'nonempty','vector','nonnan','finite','increasing'}));
           addRequired(p,'n_ghost',...
               @(x)validateattributes(x,{'numeric'},...
               {'nonempty','scalar','integer'}));
           parse(p,xi,n_ghost);
           xi = p.Results.xi;
           n_ghost = p.Results.n_ghost;
           this.grid = grid1D( xi, n_ghost );
           this.F = zeros(this.N+1,1);
           this.R = zeros(this.N,1);
           this.E = zeros(this.N,1);
           this.Rnorm = 0;
           this.Rinit = 0;
           this.count = 0;
       end
   end
end