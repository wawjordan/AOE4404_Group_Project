classdef RK_Explicit
    properties
        c, A, b, s
        Method
    end
    methods
        function this = RK_Explicit(varargin)
            default_Method = 'RK41';
            expected_names = ...
                {'Euler','RK21','RK22','RK31','RK32','RK41','RK42'};
            p = inputParser;
            addOptional(p,'Method',default_Method,...
                @(x) any(validatestring(x,expected_names)));

            parse(p,varargin{:});
            this.Method = p.Results.Method;
            switch(p.Results.Method)
                case 'Euler'
                    this.A = 0;
                    this.b = 1;
                    this.c = 0;
                    this.s = 1;
                case 'RK21'
                    this.A = [0,0;1,0];
                    this.b = [0.5;0.5];
                    this.c = [0;1];
                    this.s = 2;
                case 'RK22'
                    this.A = [0,0;0.5,0];
                    this.b = [0;1];
                    this.c = [0;0.5];
                    this.s = 2;
                case 'RK31'
                    this.A = [0,0,0; 2/3,0,0; 1/3,1/3,0];
                    this.b = [0.25;0;0.75];
                    this.c = [0;2/3;2/3];
                    this.s = 3;
                case 'RK32'
                    this.A = [0,0,0; 0.5,0,0; -1,2,0];
                    this.b = [1/6;2/3;1/6];
                    this.c = [0;0.5;1];
                    this.s = 3;
                case 'RK41'
                    this.A = diag([0.5;0.5;1],-1);
                    this.b = [1/6;1/3;1/3;1/6];
                    this.c = [0;0.5;0.5;1];
                    this.s = 4;
                case 'RK42'
                    this.A = diag([0.25;0.5;2],-1) + ...
                        diag([0;-2],-2) + diag(1,-3);
                    this.b = [1/6;0;2/3;1/6];
                    this.c = [0;0.25;0.5;1];
                    this.s = 4;
            end
        end
        function y = eval(this,dydt,yk,tk,h)
            k = zeros(size(yk,1),this.s);
            y = yk;
            for i = 1:this.s
                kj = zeros(size(yk,1),1);
                for j = 1:i-1
                    kj = kj + this.A(i,j)*k(:,j);
                end
                k(:,i) = dydt(tk+this.c(i)*h,yk+h*kj);
                y = y + h*this.b(i)*k(:,i);
            end
        end
        function [y,yk,t] = solve(this,dydt,y0,t0,tf,h)
            yk = cell(ceil((tf-t0)/h),1);
            t = (t0:h:tf)';
            N = length(t);
            y = y0;
            yk{1} = y0;
            
            for i = 1:N-1
                y = this.eval(dydt,y,t(i),h);
                yk{i+1} = y;
            end
        end
    end
end