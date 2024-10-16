function dydt = linearsolver(ts,y,k,X,params)

    ts
    omega = params.omega;

    flag_node = params.flag_node;
    inx1 = find(flag_node==1);
    F = zeros(2*6*params.N1*params.N2,1);
    F(2*params.ind_input-1) = [X(1)*sin(omega*ts);X(3)*sin(omega*ts)];
    F(2*params.ind_input) = [X(2)*sin(omega*ts);X(4)*sin(omega*ts)];
    F([2*inx1-1,2*inx1]) = [];
    C_cur = params.C(params.flag_prune_bonds==0,:);
    C_cur(:,[2*inx1-1,2*inx1]) = [];
    D = C_cur'*diag(k)*C_cur;
    M = params.M;
    M([2*inx1-1,2*inx1],:) = [];
    M(:,[2*inx1-1,2*inx1]) = [];

    dydt = zeros(length(y),1);

    dydt(1:length(y)/2,1) = y(length(y)/2+1:end);
    dydt(length(y)/2+1:end,1) = -M^(-1)*D*y(1:length(y)/2) + M^(-1)*F;

    % fixed boundary conditions
    flag_node(inx1) = [];
    inx2 = find(flag_node == 2);
    dydt(2*inx2-1,1) = 0;
    dydt(2*inx2,1) = 0;
    dydt(2*inx2-1+size(M,1),1) = 0;
    dydt(2*inx2+size(M,1),1) = 0;

end