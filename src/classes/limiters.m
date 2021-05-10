classdef limiters < handle
    properties
        scheme
        label
        limiterfun
        beta
        frozen = false
    end
    properties (Dependent = true)%, SetAccess = private)
        limiter
    end
    methods
        function this = limiters(varargin)
            default_scheme = 'none';
            default_beta = 2.0;
            expected_names = ...
                {'none','van_leer','van_albada','minmod','beta_lim'};
            p = inputParser;
            addOptional(p,'scheme',default_scheme,...
                @(x) any(validatestring(x,expected_names)));
            addOptional(p,'beta',default_beta,...
                @(x)validateattributes(x,{'numeric'},...
                {'nonempty','positive'}));
            parse(p,varargin{:});
            this.scheme = p.Results.scheme;
            this.beta = p.Results.beta;
            switch(p.Results.scheme)
                case 'none'
                    this.limiterfun = @no_limiter;
                    this.label = 'no limiter';
                case 'van_leer'
                    this.limiterfun = @van_leer_limiter;
                    this.label = 'van Leer limiter';
                case 'van_albada'
                    this.limiterfun = @van_albada_limiter;
                    this.label = 'van Albada limiter';
                case 'minmod'
                    this.limiterfun = @minmod_limiter;
                    this.label = 'minmod limiter';
                case 'beta_lim'
                    this.limiterfun = @beta_limiter;
                    this.label = ...
                    sprintf('$\\beta$-limiter $(\\beta=%0.1f)$',this.beta);
            end
        end
        function set.limiter(this,r)
            if ~this.frozen
                this.limiter = this.limiterfun(this,r);
            end
        end
        function set.frozen(this,val)
           this.frozen = val;
        end
        function psi = limit(this,r)
          psi = this.limiterfun(this,r);
        end
        function psi = no_limiter(~,r)
            psi = ones(size(r,1),size(r,2));
        end
        function psi = van_leer_limiter(~,r)
            psi = (r + abs(r))./(1+r);
            psi(r<0)=0;
        end
        function psi = van_albada_limiter(~,r)
            psi = (r.^2 + r)./(1 + r.^2);
            psi(r<0)=0;
        end
        function psi = minmod_limiter(~,r)
            psi = min(r,1);
            psi(r<0)=0;
        end
        function psi = beta_limiter(this,r)
            psi = max(0,max(min(this.beta*r,1),min(r,this.beta)));
        end
    end
end