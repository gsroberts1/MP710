function [F] = t1_model(theta,ti)
    pd_mod = theta(1);
    r1_mod = theta(2);
    
    mz = pd_mod.*(1-2*exp(-ti.*r1.*mod));
    F = abs(mz);
end

