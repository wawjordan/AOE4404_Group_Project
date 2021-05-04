function [left,right] = MUSCL_extrap(soln,limiter)
    k = soln.kappa_MUSCL;
    e = soln.epsilon_MUSCL;
    
    delta = 1e-6;
    r_plus = 0.0*soln.U;
    r_minus = 0.0*soln.U;
    i = soln.i_low-1:soln.i_high;
    
    i_low2 = soln.i_low:soln.i_low+soln.grid.n_ghost-1;
    i_high2 = soln.i_high-soln.grid.n_ghost+1:soln.i_high;
    
    soln.U(soln.i_high+1:end) = soln.U(i_low2);
    soln.U(1:soln.i_low-1) = soln.U(i_high2);
    
    den = soln.U(i+1) - soln.U(i);
    den = sign(den).*max(abs(den),delta);
    den(den==0)= delta;
    
    num = soln.U(i+2) - soln.U(i+1);
    r_plus(i) = num./den;
    num = soln.U(i) - soln.U(i-1);
    r_minus(i) = num./den;
    
%     r_plus(soln.i_low-2) = r_plus(soln.i_low-1);
%     r_minus(soln.i_low-2) = r_minus(soln.i_low-1);
% 
%     r_plus(soln.i_high+1) = r_plus(soln.i_high);
%     r_minus(soln.i_high+1) = r_minus(soln.i_high);

    r_plus(soln.i_low-2) = r_plus(soln.i_high);
    r_minus(soln.i_low-2) = r_minus(soln.i_high);

    r_plus(soln.i_high+1) = r_plus(soln.i_low);
    r_minus(soln.i_high+1) = r_minus(soln.i_low);

    
    psi_plus = limiter.limit(r_plus);
    psi_minus = limiter.limit(r_minus);
    
    left = soln.U(i) + (e/4)*(...
        (1-k)*psi_plus(i-1).*(soln.U(i)-soln.U(i-1)) + ...
        (1+k)*psi_minus(i).*(soln.U(i+1)-soln.U(i)));
    right = soln.U(i+1) - (e/4)*(...
        (1-k)*psi_minus(i+1).*(soln.U(i+2)-soln.U(i+1)) + ...
        (1+k)*psi_plus(i).*(soln.U(i+1)-soln.U(i)));
end