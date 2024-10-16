function [flag_prune_bonds,flag_node] = pruneBott(pos,bonds,params)
    % prune some bonds and nodes
    flag_node = zeros(size(pos,1),1);
    % remove nodes and bonds for the boundary
    ind_remove_node = [];
    row_remove_node = 1:params.N2;col_remove_node = [1:params.N1/4,3*params.N1/4+1:params.N1];
    for i = 1:length(row_remove_node)
        for j = 1:length(col_remove_node)
            ind_remove_node = [ind_remove_node,[1:6]+6*(row_remove_node(i)-1)*params.N1+6*(col_remove_node(j)-1)];
        end
    end
    row_remove_node = [3*params.N2/4+1:params.N2];col_remove_node = [params.N1/4:3*params.N1/4];
    for i = 1:length(row_remove_node)
        for j = 1:length(col_remove_node)
            ind_remove_node = [ind_remove_node,[1:6]+6*(row_remove_node(i)-1)*params.N1+6*(col_remove_node(j)-1)];
        end
    end
    row_remove_node = [1];col_remove_node = [params.N1/4:3*params.N1/4];
    for i = 1:length(row_remove_node)
        for j = 1:length(col_remove_node)
            ind_remove_node = [ind_remove_node,[1:6]+6*(row_remove_node(i)-1)*params.N1+6*(col_remove_node(j)-1)];
        end
    end
    ind_remove_node = unique(ind_remove_node);
    % set boundary
    row_boundary_node = 1:3*params.N2/4;col_boundary_node = [params.N1/4,3*params.N1/4+1];
    ind_boundary_node = [];
    for i = 1:length(row_boundary_node)
        for j = 1:length(col_boundary_node)
            if col_boundary_node(j) == params.N1/4
                if row_boundary_node(i) == 2
                    ind_boundary_node = [ind_boundary_node,[1]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                elseif row_boundary_node(i) == 1
                    ind_boundary_node = [ind_boundary_node,[]];
                else
                    ind_boundary_node = [ind_boundary_node,[1,6]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                end
            else
                if row_boundary_node(i) == 3*params.N2/4
                    ind_boundary_node = [ind_boundary_node,[4]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                elseif row_boundary_node(i) == 1
                    ind_boundary_node = [ind_boundary_node,[]];
                else
                    ind_boundary_node = [ind_boundary_node,[3,4]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                end
            end
        end
    end
    row_boundary_node = [1,3*params.N2/4+1];col_boundary_node = params.N1/4:3*params.N1/4;
    for i = 1:length(row_boundary_node)
        for j = 1:length(col_boundary_node)
            if row_boundary_node(i) == 1
                if col_boundary_node(j) == params.N1/4
                    ind_boundary_node = [ind_boundary_node,[]];
                elseif col_boundary_node(j) == params.N1/4+1
                    ind_boundary_node = [ind_boundary_node,[2]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                elseif col_boundary_node(j) == 3*params.N1/4
                    ind_boundary_node = [ind_boundary_node,[2,3]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1),...
                        [3]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1+1)];
                else
                    ind_boundary_node = [ind_boundary_node,[2,3]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                end
            else
                if col_boundary_node(j) == params.N1/4
                    ind_boundary_node = [ind_boundary_node,[6]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                elseif col_boundary_node(j) == 3*params.N1/4
                    ind_boundary_node = [ind_boundary_node,[5]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                else
                    ind_boundary_node = [ind_boundary_node,[5,6]+6*(row_boundary_node(i)-1)*params.N1+6*(col_boundary_node(j)-1)];
                end
            end
        end
    end
    ind_boundary_node = unique(ind_boundary_node);
    inx = [];
    for i = 1:length(ind_boundary_node)
        inx = [inx,find(ind_boundary_node(i) == ind_remove_node)];
    end
    ind_remove_node(inx) = [];
    % prune corresponding bonds
    flag_prune_bonds = zeros(size(bonds,1),1);
    ind_remove_bond = [];
    for i = 1:length(ind_remove_node)
        ind_remove_bond = [ind_remove_bond;find(bonds(:,1) == ind_remove_node(i));find(bonds(:,2) == ind_remove_node(i))];
    end
    for i = 1:length(ind_boundary_node)
        for j = 1:length(ind_boundary_node)
            if ind_boundary_node(i) ~= ind_boundary_node(j)
                [q1,inx1] = ismember([ind_boundary_node(i),ind_boundary_node(j)],bonds,'rows');
                [q2,inx2] = ismember([ind_boundary_node(j),ind_boundary_node(i)],bonds,'rows');
                if q1 == 1
                    ind_remove_bond = [ind_remove_bond;inx1];
                elseif q2 == 1
                    ind_remove_bond = [ind_remove_bond;inx2];
                end
            end
        end
    end
    flag_prune_bonds(ind_remove_bond) = 1;
    flag_node(ind_remove_node) = 1;
    flag_node(ind_boundary_node) = 2;
end