function [x,xk,ea] = newton_safe(f,df,x0,a,b,tol,maxiter)
xk = nan(maxiter,1);
ea = nan(maxiter,1);

xk(1) = x0;
xk(2) = xk(1) - f(xk(1))/df(xk(1));
k = 2;
ea(k) = abs(xk(k)-xk(k-1))/abs(xk(k));
while (ea(k)>tol)&&(k <= maxiter)
    xk(k+1) = xk(k) - f(xk(k))/df(xk(k));
    if ( xk(k+1) < a )||( xk(k+1) > b)
        xk(k+1) = a +(b-a)/2;
        if sign(f(a)) == sign(f(xk(k+1)))
            a = xk(k+1);
        else
            b = xk(k+1);
        end
    end
    ea(k+1) = abs(xk(k+1)-xk(k))/abs(xk(k+1));
    k = k + 1;
end
x = xk(k);
xk = xk(~isnan(xk));
ea = ea(~isnan(xk));
end