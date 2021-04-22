classdef soln1D < handle
   properties
      grid, IC
      N, Lref {mustBeNumeric}
      dt, t0, tf {mustBeNumeric}
      U, F, E, R {mustBeNumeric}
      Rnorm, Rinit, count {mustBeNumeric}
   end
   methods
       function this = soln1D(xi,varargin)
           this.N = length(xi)-1;
           default_ghost = 0;
           default_time_range = [0,1];
           default_dt = 0.01;
           default_Lref = max(xi)-min(xi);
           default_IC = @(x) zeros(length(x),1);
           
           p = inputParser;
           valid_func = @(x) isa(x,'function_handle');
           addRequired(p,'xi',@(x) validateattributes(x,{'numeric'},{'nonempty','vector','nonnan','finite','increasing'}));
           addParameter(p,'n_ghost',default_ghost,@(x) validateattributes(x,{'numeric'},{'nonempty','scalar','integer'}));
           addParameter(p,'IC',default_IC,@(x)(valid_func(x)));
           addParameter(p,'dt',default_dt,@(x) validateattributes(x,{'numeric'},{'nonempty','positive'}));
           addParameter(p,'Lref',default_Lref,@(x) validateattributes(x,{'numeric'},{'nonempty','positive'}));
           addParameter(p,'time_range',default_time_range,@(x) validateattributes(x,{'numeric'},{'nonempty','vector','numel',2,'nonnan','finite'}));
           
           parse(p,xi,varargin{:});
           xi = p.Results.xi;
           n_ghost = p.Results.n_ghost;
           this.grid = grid1D( xi, n_ghost );
           this.Lref = p.Results.Lref;
           this.IC = p.Results.IC;
           this.dt = p.Results.dt;
           this.t0 = p.Results.time_range(1);
           this.tf = p.Results.time_range(2);
           this.U = this.IC(this.grid.xc);
           this.F = zeros(this.N+1,1);
           this.R = zeros(this.N,1);
           this.E = zeros(this.N,1);
           this.Rnorm = 0;
           this.Rinit = 0;
           this.count = 0;
       end
   end
end