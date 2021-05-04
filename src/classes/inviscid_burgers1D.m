classdef inviscid_burgers1D < soln1D
    properties
      exact_solution_type
      i_low, i_high, i {mustBeNumeric}
      t0, tf, t {mustBeNumeric}
      CFL, eps, Lref {mustBeNumeric}
      order {mustBeNumeric}
      kappa_MUSCL {mustBeNumeric}
      epsilon_MUSCL {mustBeNumeric}
      uLeft, uRight, L {mustBeNumeric}
    end
    properties (Dependent = true, SetAccess = private)
        
    end
    methods
        function this = inviscid_burgers1D(xi,n_ghost,varargin)
            this = this@soln1D(xi,n_ghost);
            this.i_low = this.grid.i_low;
            this.i_high = this.grid.i_high;
            this.i = (this.grid.i_low:this.grid.i_high)';
            
            default_time_range = [0,1];
            default_dt = 0.01;
            default_Lref = max(xi)-min(xi);
            default_IC = @(x) zeros(length(x),1);
            
            default_exact_solution_type = 'initial_disc';
            default_uLeft = 1;
            default_uRight = 0;
            default_L = 0.25*(this.grid.xmax-this.grid.xmin);
            default_order = 1;
            default_kappa = 1;
            default_CFL = 0.5;
            expected_solutions = {'initial_linear','initial_disc',...
                'initial_sine'};
            p = inputParser;
            addParameter(p,'IC',default_IC,@(x)(valid_func(x)));
            addParameter(p,'dt',default_dt,...
                @(x)validateattributes(x,{'numeric'},...
                {'nonempty','positive'}));
            addParameter(p,'Lref',default_Lref,...
                @(x)validateattributes(x,{'numeric'},...
                {'nonempty','positive'}));
            addParameter(p,'time_range',default_time_range,...
                @(x)validateattributes(x,{'numeric'},...
                {'nonempty','vector','numel',2,'nonnan','finite'}));
            addParameter(p,'order',default_order,@(x)(x==1)||(x==2));
            addParameter(p,'kappa',default_kappa,...
                @(x)validateattributes(x,{'numeric'},{'nonempty'}));
            addParameter(p,'CFL',default_CFL,...
                @(x)validateattributes(x,{'numeric'},...
                {'nonempty','positive'}));
            addParameter(p,'exact_solution_type',...
                default_exact_solution_type,...
                @(x)any(validatestring(x,expected_solutions)));
            addParameter(p,'uLeft',default_uLeft,...
                @(x)validateattributes(x,{'numeric'},{'nonempty'}));
            addParameter(p,'uRight',default_uRight,...
                @(x)validateattributes(x,{'numeric'},{'nonempty'}));
            addParameter(p,'L',default_L,...
                @(x)validateattributes(x,{'numeric'},...
                {'nonempty','positive'}));
            parse(p,varargin{:});
            
            this.Lref = p.Results.Lref;
            this.IC = p.Results.IC;
            this.dt = p.Results.dt;
            this.t0 = p.Results.time_range(1);
            this.tf = p.Results.time_range(2);
            this.U = this.IC(this.grid.xc);
            this.eps = 0.001*min(this.dt);
            this.kappa_MUSCL = p.Results.kappa;
            this.order = p.Results.order;
            this.epsilon_MUSCL = p.Results.order - 1;
            this.uLeft = p.Results.uLeft;
            this.uRight = p.Results.uRight;
            this.L = p.Results.L;
            switch(p.Results.exact_solution_type)
               case 'initial_linear'
                   this.exact_solution_type = @initial_linear;
                   this.IC = @(x)this.exact_solution_type(this,x,0);
                   this.U = this.IC(this.grid.xc);
                   this.t0 = 0;
                   this.tf = 1;
               case 'initial_disc'
                   this.exact_solution_type = @initial_disc;
                   this.IC = @(x)this.exact_solution_type(this,x,0);
                   this.U = this.IC(this.grid.xc);
                   this.t0 = 0;
                   this.tf = 1;
               case 'initial_sine'
                   this.exact_solution_type = @initial_sine;
                   this.IC = @(x)this.exact_solution_type(this,x,0);
                   this.U = this.IC(this.grid.xc);
                   this.t0 = 0;
                   this.tf = 0.99/pi;
            end
            this.t = this.t0;
        end
        function uex = calc_exact(this,x,t)
          uex = this.exact_solution_type(this,x,t);
        end
        function uex = initial_linear(this,x,t)
            u1 = 1*(x-t<this.eps);
            u2 = (1-x)./(1-t).*(x-t>this.eps&x-1<this.eps);
            u2(isnan(u2))=0;
            u3 = 1*(x-(t+1)/2<this.eps);
            uex = (u1+u2).*(t-1<=this.eps) + u3*(t-1>this.eps);
        end
        function uex = initial_disc(this,x,t)
            uL = this.uLeft;
            uR = this.uRight;
            xmid = this.grid.xmin + 0.5*(this.grid.xmax-this.grid.xmin);
            us1 = uL*(x<xmid+0.5*(uL+uR)*t);
            us2 = uR*(x>=xmid+0.5*(uL+uR)*t);
            u_shock = us1+us2;
            uf1 = uL*((x-xmid)/t<=uL);
            uf1(isnan(uf1)) = uL;
            uf2 = ((x-xmid)/t).*((x-xmid)/t>uL & (x-xmid)/t<uR);
            uf2(isnan(uf2)) = 0;
            uf3 = uR*((x-xmid)/t>=uR);
            uf3(isnan(uf3)) = uR;
            u_fan = uf1 + uf2 + uf3;
            uex = u_shock*(uL>uR) + u_fan*(uL<uR);
        end
        function uex = initial_sine(this,x,t)
            a = 0.5*(this.grid.xmax-this.grid.xmin);
            tol = 1e-14;
            maxiter = 100;
            init = @(a,x) a - sin(pi*x);
            f = @(eta,x,t,a) eta + (a-sin(pi*eta))*t-x;
            df = @(eta,t) 1-pi*t*cos(pi*eta);
            N = length(x);
            xmax = max(x);
            xmin = min(x);
            x2 = x + init(a,x)*t;
            xnew = [x(x2>xmax);x(x2<=xmax)];
            uex = zeros(N,1);
            for j = 1:N
                left = xmin-a;
                right = x(j);
                [eta,~,~] = newton_safe(@(eta)f(eta,x(j),t,a),...
                    @(eta)df(eta,t),xnew(j),left,right,tol,maxiter);
                uex(j) = init(a,eta);
            end
        end
        function res = residual(this,F)
           res = F(2:this.grid.i_max+1)-F(1:this.grid.i_max);
       end
    end
end