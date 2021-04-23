classdef fluxes
    properties
        scheme
        fluxfun
    end
    methods
        function this = fluxes(varargin)
            default_scheme = 'eo';
            expected_names = ...
                {'central','godunov','eo','roe'};
            p = inputParser;
            addOptional(p,'scheme',default_scheme,...
                @(x) any(validatestring(x,expected_names)));

            parse(p,varargin{:});
            this.scheme = p.Results.scheme;
            switch(p.Results.scheme)
                case 'central'
                    this.fluxfun = @central_flux;
                case 'godunov'
                    this.fluxfun = @godunov_flux;
                case 'eo'
                    this.fluxfun = @eo_flux;
                case 'roe'
                    this.fluxfun = @roe_flux;
            end
        end
        function F = calc_flux(this,left,right)
          F = this.fluxfun(this,left,right);
        end
        function F = central_flux(~,left,right)
            ui = 0.5*(left+right);
            F = 0.5*ui.^2;
        end
        function F = godunov_flux(~,left,right)
            uplusL = max(left,0);
            uminusR = min(right,0);
            F = max(0.5*uminusR.^2,0.5*uplusL.^2);
        end
        function F = eo_flux(~,left,right)
            FplusL = 0.5*(max(left,0).^2);
            FminusR = 0.5*(min(right,0).^2);
            F = FplusL + FminusR;
        end
        function F = roe_flux(~,left,right)
            a_bar = 0.5*(right+left);
            eps = max(0,max((a_bar-left),(right-a_bar)));
            a_bar_mod1 = abs(a_bar);
            a_bar_mod2 = eps;
%             a_bar_mod2 = 0.5*((a_bar.^2)./eps + eps);
%             a_bar_mod2(isnan(a_bar_mod2)|isinf(a_bar_mod2)) = 0;
            a_bar_mod = a_bar_mod1.*(a_bar_mod1>=eps) ...
                + a_bar_mod2.*(a_bar_mod1<eps);
            FL = 0.5*left.^2;
            FR = 0.5*right.^2;
            F = 0.5*(FL + FR) - 0.5*a_bar_mod.*(right-left);
        end
    end
end