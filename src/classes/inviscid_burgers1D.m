classdef inviscid_burgers1D < soln1D
    properties
      exact_solution_type
      i_low, i_high, i
      t {mustBeNumeric}
      CFL, eps {mustBeNumeric}
    end
    properties (Dependent = true, SetAccess = private)
        
    end
    methods
        function this = inviscid_burgers1D(xi,inputs,varargin)
            this = this@soln1D(xi,inputs);
            this.i_low = this.grid.i_low;
            this.i_high = this.grid.i_high;
            this.i = (this.grid.i_low:this.grid.i_high)';
            default_exact_solution_type = 'ramp_shock';
            default_CFL = 0.5;
            expected_solutions = {'ramp_shock'};
            p = inputParser;
%             valid_func = @(x) isa(x,'function_handle');
            addParameter(p,'CFL',default_CFL,@(x) validateattributes(x,{'numeric'},{'nonempty','positive'}));
            addParameter(p,'exact_solution_type',default_exact_solution_type,...
               @(x) any(validatestring(x,expected_solutions)));
            parse(p,varargin{:});
            
            this.eps = 0.001*min(this.dt);
            switch(p.Results.exact_solution_type)
               case 'ramp_shock'
                   this.exact_solution_type = @ramp_shock;
                   this.IC = @(x)this.exact_solution_type(this,x,0);
                   this.U = this.IC(this.grid.xc);
                   this.t0 = 0;
                   this.tf = 1;
            end
            this.t = this.t0;
        end
        function uex = calc_exact(this,x,t)
          uex = this.exact_solution_type(this,x,t);
        end
        function uex = ramp_shock(this,x,t)
            u1 = 1*(x-t<this.eps);
            u2 = (1-x)./(1-t).*(x-t>this.eps&x-1<this.eps);
            u2(isnan(u2))=0;
            u3 = 1*(x-(t+1)/2<this.eps);
            uex = (u1+u2).*(t-1<=this.eps) + u3*(t-1>this.eps);
        end
        function res = residual(this,F)
           res = F(2:this.grid.i_max+1)-F(1:this.grid.i_max);
       end
    end
end