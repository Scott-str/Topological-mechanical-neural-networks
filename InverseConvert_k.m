function k = InverseConvert_k(k0,params)
    k = k0;
    for i = 1:length(k)
        if params.UB(i) == inf && params.LB(i) ~= inf
            k(i) = log(exp(k(i)-params.LB(i))-1);
        elseif params.UB(i) ~= inf && params.LB(i) == inf
            k(i) = log(exp(-(k(i)-params.UB(i)))-1);
        elseif params.UB(i) ~= inf && params.LB(i) ~= inf
            if params.UB(i) ~= params.LB(i)
                k(i) = -1/params.beta*log((params.UB(i)-params.LB(i))/(k(i)-params.LB(i))-1);
            elseif params.UB(i) == params.LB(i)
                k(i) = 0;
            end
        else
            k(i) = k(i);
        end
    end