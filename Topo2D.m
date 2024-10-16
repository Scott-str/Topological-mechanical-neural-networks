function [fcost,u_sol,e_sol,F_err] = Topo2D(k0,F,u_output,params,optims,condition)
    
    if isequal(condition,'1')
        k = Convert_k(k0,params);
    elseif isequal(condition,'0')
        k = k0;
    end
    flag_node = params.flag_node;
    inx1 = find(flag_node==1);
    F([2*inx1-1,2*inx1]) = [];
    C_cur = params.C(params.flag_prune_bonds==0,:);
    C_cur(:,[2*inx1-1,2*inx1]) = [];
    D = C_cur'*diag(k)*C_cur;
    M = params.M;
    M([2*inx1-1,2*inx1],:) = [];
    M(:,[2*inx1-1,2*inx1]) = [];
    
    % boundary conditions
    flag_node(inx1) = [];
    inx2 = find(flag_node == 2);
    F([2*inx2-1,2*inx2]) = [];
    D(:,[2*inx2-1,2*inx2]) = [];
    D([2*inx2-1,2*inx2],:) = []; 
    M(:,[2*inx2-1,2*inx2]) = [];
    M([2*inx2-1,2*inx2],:) = [];
    u_free = (-params.omega^2*M+D)\F;

    u_cur = zeros(size(C_cur,2),1);
    inx3 = find(flag_node == 0);
    u_cur(2*inx3-1) = u_free(1:2:end);
    u_cur(2*inx3) = u_free(2:2:end);
    e_sol = C_cur*u_cur;
    
    u_sol = zeros(2*6*params.N1*params.N2,1);
    inx4 = find(params.flag_node == 0);
    u_sol(2*inx4-1) = u_free(1:2:end);
    u_sol(2*inx4) = u_free(2:2:end);
    
    u_sol_output_1 = zeros(2*length(params.ind_output_1),1);
    u_sol_output_1(1:2:end) = u_sol(2*params.ind_output_1-1);
    u_sol_output_1(2:2:end) = u_sol(2*params.ind_output_1);

    u_sol_output_2 = zeros(2*length(params.ind_output_2),1);
    u_sol_output_2(1:2:end) = u_sol(2*params.ind_output_2-1);
    u_sol_output_2(2:2:end) = u_sol(2*params.ind_output_2);

    u_sol_output = [u_sol_output_1;u_sol_output_2];

    d_sol_output_1 = u_sol_output_1(1:2:end)'*u_sol_output_1(1:2:end) + u_sol_output_1(2:2:end)'*u_sol_output_1(2:2:end);
    d_sol_output_2 = u_sol_output_2(1:2:end)'*u_sol_output_2(1:2:end) + u_sol_output_2(2:2:end)'*u_sol_output_2(2:2:end);
    d_sol_output_norm = [d_sol_output_1/(d_sol_output_1+d_sol_output_2);d_sol_output_2/(d_sol_output_1+d_sol_output_2)];

    % cross entropy
    fcost = -log(u_output*d_sol_output_norm);

    F_err = zeros(2*6*params.N1*params.N2,1);

    v = sym('v',[1,length(u_output)],'real');
    ux = sym('ux',[length(params.ind_output_1)+length(params.ind_output_2),1],'real');
    uy = sym('uy',[length(params.ind_output_1)+length(params.ind_output_2),1],'real');
    u_sym = sym(zeros(2*(length(params.ind_output_1)+length(params.ind_output_2)),1));
    u_sym(1:2:end) = ux;
    u_sym(2:2:end) = uy;
    u1 = ux(1:length(ux)/2)'*ux(1:length(ux)/2)+uy(1:length(uy)/2)'*uy(1:length(uy)/2);
    u2 = ux(length(ux)/2+1:end)'*ux(length(ux)/2+1:end)+uy(length(uy)/2+1:end)'*uy(length(uy)/2+1:end);
    u_norm = [u1/(u1+u2);u2/(u1+u2)];

    J = jacobian(-log(v*u_norm),u_sym);
    
    F_sym = -J;
    F_sym = subs(F_sym,u_sym,u_sol_output);
    F_sym = subs(F_sym,v,u_output);
    F_sym = double(F_sym);
    F_err(2*params.ind_output_1-1) = F_sym(1:2:length(F_sym)/2);
    F_err(2*params.ind_output_1) = F_sym(2:2:length(F_sym)/2);
    F_err(2*params.ind_output_2-1) = F_sym(length(F_sym)/2+1:2:end);
    F_err(2*params.ind_output_2) = F_sym(length(F_sym)/2+2:2:end);
