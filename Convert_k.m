function k = Convert_k(k0,params)
    k = k0;
    for i = 1:length(k)
        if params.UB(i) == inf && params.LB(i) ~= inf
            k(i) = log(1+exp(k(i))) + params.LB(i);
        elseif params.UB(i) ~= inf && params.LB(i) == inf
            k(i) = -log(1+exp(k(i))) + params.UB(i);
        elseif params.UB(i) ~= inf && params.LB(i) ~= inf
            if params.UB(i) ~= params.LB(i)
                k(i) = params.LB(i) + (params.UB(i)-params.LB(i))./(1+exp(-params.beta*k(i)));
            elseif params.UB(i) == params.LB(i)
                k(i) = params.LB(i);
            end
        else
            k(i) = k(i);
        end
    end