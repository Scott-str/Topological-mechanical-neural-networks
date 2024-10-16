function [pred,acc] = prediction(k,params,optims,X,y)
    pred = [];
    parfor i = 1:size(X,1)
        F = zeros(2*6*params.N1*params.N2,1);
        F(2*params.ind_input-1) = X(i,[1,3]);
        F(2*params.ind_input) = X(i,[2,4]);
        u_output = y(i,:);
        [~,pred_temp,~,~] = Topo2D(k,F,u_output,params,optims,'0');
        pred_temp_x_1 = pred_temp(2*params.ind_output_1-1);
        pred_temp_y_1 = pred_temp(2*params.ind_output_1);
        pred_temp_x_2 = pred_temp(2*params.ind_output_2-1);
        pred_temp_y_2 = pred_temp(2*params.ind_output_2);
        pred(i,:) = [pred_temp_x_1'*pred_temp_x_1+pred_temp_y_1'*pred_temp_y_1,pred_temp_x_2'*pred_temp_x_2+pred_temp_y_2'*pred_temp_y_2];
    end
    [~,ind_pred] = max(pred,[],2);
    [~,ind_y] = max(y,[],2);
    acc = sum(ind_pred == ind_y)/length(ind_pred);