function [k,m_moment,v_moment] = Train_Adam(m0,sgCurr,m,v,iter,params,optims)
    
    alpha = optims.alpha;
    beta1 = optims.beta1;
    beta2 = optims.beta2;
        
    % Update biased 1st moment estimate
    m_moment = beta1.*m + (1 - beta1).*sgCurr;
    % Update biased 2nd raw moment estimate
    v_moment = beta2.*v + (1 - beta2).*(sgCurr.^2);
        
    % Compute bias-corrected 1st moment estimate
    mHat = m_moment./(1 - beta1^iter);
    % Compute bias-corrected 2nd raw moment estimate
    vHat = v_moment./(1 - beta2^iter);
        
    % Update decision variables
    k = m0 - alpha.*mHat./(sqrt(vHat) + sqrt(eps));