classdef limiters
    properties
        scheme
        limiterfun
        beta
    end
    properties (Dependent = true, SetAccess = private)
        frozen
        unfrozen
    end
    methods
        function this = limiters(varargin)
            default_scheme = 'minmod';
            default_beta = 2.0;
            expected_names = ...
                {'van_leer','van_albada','minmod','beta_lim'};
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
                case 'van_leer'
                    this.limiterfun = @van_leer_limiter;
                case 'van_albada'
                    this.limiterfun = @van_albada_limiter;
                case 'minmod'
                    this.limiterfun = @minmod_limiter;
                case 'beta_lim'
                    this.limiterfun = @beta_limiter;
            end
        end
        function psi = limit(this,r)
          psi = this.limiterfun(this,r);
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