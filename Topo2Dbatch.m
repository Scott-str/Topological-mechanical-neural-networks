function [fcost,sgCurr] = Topo2Dbatch(k0,X_batch,y_batch,params,optims,condition)
    fGrad = [];fcost = 0;
    parfor i = 1:size(X_batch,1)
        F = zeros(2*6*params.N1*params.N2,1);
        F(2*params.ind_input-1) = X_batch(i,[1,3]);
        F(2*params.ind_input) = X_batch(i,[2,4]);
        u_output = y_batch(i,:);
        [fcost_temp,~,eori,F_err] = Topo2D(k0,F,u_output,params,optims,condition);
        [~,~,eadj,~] = Topo2D(k0,F_err,u_output,params,optims,condition);
        fGrad(:,i) = eadj.*eori;
        fcost = fcost + fcost_temp;
    end
    fcost = fcost/optims.batch_size;
    sgCurr = mean(fGrad,2);
