classdef grid1D < handle
    properties
        xi, xc, dx {mustBeNumeric}
        n_ghost, i_low, i_high {mustBeNumeric}
        xmin, xmax {mustBeNumeric}
        i_max {mustBeNumeric}
        L {mustBeNumeric}
    end
    methods
        function this = grid1D( xi, n_ghost )
            this.xmin = min(xi);
            this.xmax = max(xi);
            this.i_max = length(xi)-1;
            this.L = max(xi)-min(xi);
            if nargin < 2
                this.n_ghost = 0;
                this.i_low = 1;
                this.i_high = this.i_max;
                this.xi = xi(:);
                xplus = [(2*xi(1)-xi(2));xi(:);(2*xi(end)-xi(end-1))];
                this.dx = 0.5*(diff(xplus(1:end-1))+diff(xplus(2:end)));
            else
                this.n_ghost = n_ghost;
                this.i_low = 1+n_ghost;
                this.i_high = this.i_max+n_ghost;
                dxlow = xi(1)-xi(2);
                dxhigh = xi(end)-xi(end-1);
                xlow = flipud((xi(1)+dxlow:dxlow:xi(1)+(n_ghost)*dxlow)');
                xhigh = (xi(end)+dxhigh:dxhigh:xi(end)+(n_ghost)*dxhigh)';
                this.xi = [xlow;xi(:);xhigh];
                xplus = [(2*this.xi(1)-this.xi(2));this.xi(:);...
                    (2*this.xi(end)-this.xi(end-1))];
                this.dx = 0.5*(diff(xplus(1:end-1))+diff(xplus(2:end)));
            end
            this.xc = 0.5*(this.xi(1:this.i_max-1)+this.xi(2:this.i_max));
        end
    end
end