function [x,xk,ea] = newton_safe(f,df,x0,a,b,tol,maxiter)
xk = nan(maxiter,1);
ea = nan(maxiter,1);

xk(1) = x0;
fk1 = f(xk(1));
xk(2) = xk(1) - fk1/df(xk(1));
k = 2;
ea(k) = 1;
while (ea(k)>tol)&&(k <= maxiter)
    fk = f(xk(k));
    xk(k+1) = xk(k) - fk/df(xk(k));
    if ( xk(k+1) < a )||( xk(k+1) > b)
        xk(k+1) = a +(b-a)/2;
        if sign(f(a)) == sign(f(xk(k+1)))
            a = xk(k+1);
        else
            b = xk(k+1);
        end
    end
    ea(k+1) = abs(fk)/abs(fk1);
    k = k + 1;
end
x = xk(k);
xk = xk(~isnan(xk));
ea = ea(~isnan(xk));
end