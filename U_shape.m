function [pos,bonds,flag_bonds] = U_shape(params)

    a1 = [3*params.a,0];
    a2 = [3*params.a*cos(pi/3),3*params.a*sin(pi/3)];
    pos_base = [params.a*cos(0:pi/3:2*pi-pi/3)',params.a*sin(0:pi/3:2*pi-pi/3)'];
    intra_bonds_base = [1,2;2,3;3,4;4,5;5,6;6,1];
    pos = [];intra_bonds = [];inter_bonds = [];
    flag_intra = [];flag_inter = [];
    for i = 1:params.N2
        for j = 1:params.N1
            pos = [pos;pos_base+(j-1)*a1+(i-1)*a2];
            inter_ind1 = 1+6*(i-1)*params.N1+6*(j-1);
            inter_ind2 = 4+6*(i-1)*params.N1+6*(j-1+1);
            inter_ind3 = 2+6*(i-1)*params.N1+6*(j-1);
            inter_ind4 = 5+6*(i-1+1)*params.N1+6*(j-1);
            inter_ind5 = 3+6*(i-1)*params.N1+6*(j-1);
            inter_ind6 = 6+6*(i-1+1)*params.N1+6*(j-1-1);
            intra_bonds = [intra_bonds;intra_bonds_base+6*(i-1)*params.N1+6*(j-1)];
            if j <= params.N1/4
                flag_intra = [flag_intra;ones(size(intra_bonds_base+6*(i-1)*params.N1+6*(j-1),1),1)];
            elseif j > params.N1/4 && j <= 3*params.N1/4
                if i <= 3*params.N2/4
                    flag_intra = [flag_intra;zeros(size(intra_bonds_base+6*(i-1)*params.N1+6*(j-1),1),1)];
                else
                    flag_intra = [flag_intra;ones(size(intra_bonds_base+6*(i-1)*params.N1+6*(j-1),1),1)];
                end
            else
                flag_intra = [flag_intra;ones(size(intra_bonds_base+6*(i-1)*params.N1+6*(j-1),1),1)];
            end
            if i == params.N2
                if j ~= params.N1
                    inter_bonds = [inter_bonds;[inter_ind1,inter_ind2]];
                    flag_inter = [flag_inter;0];
                end
            else
                if j == 1
                    inter_bonds = [inter_bonds;[inter_ind1,inter_ind2];[inter_ind3,inter_ind4]];
                    flag_inter = [flag_inter;[0;0]];
                elseif j == params.N1
                    inter_bonds = [inter_bonds;[inter_ind3,inter_ind4];[inter_ind5,inter_ind6]];
                    flag_inter = [flag_inter;[0;0]];
                else
                    inter_bonds = [inter_bonds;[inter_ind1,inter_ind2];[inter_ind3,inter_ind4];[inter_ind5,inter_ind6]];
                    if i < 3*params.N2/4
                        if j < params.N1/4
                            flag_inter = [flag_inter;[0;0;0]];
                        elseif j == params.N1/4
                            flag_inter = [flag_inter;[2;0;0]];
                        elseif j == params.N1/4+1
                            flag_inter = [flag_inter;[1;1;2]];
                        elseif j > params.N1/4+1 && j < 3*params.N1/4
                            flag_inter = [flag_inter;[1;1;1]];
                        elseif j == 3*params.N1/4
                            flag_inter = [flag_inter;[2;1;1]];
                        elseif j == 3*params.N1/4+1
                            flag_inter = [flag_inter;[0;0;2]];
                        else
                            flag_inter = [flag_inter;[0;0;0]];
                        end
                    elseif i == 3*params.N2/4
                        if j < params.N1/4
                            flag_inter = [flag_inter;[0;0;0]];
                        elseif j == params.N1/4
                            flag_inter = [flag_inter;[2;0;0]];
                        elseif j >= params.N1/4+1 && j < 3*params.N1/4
                            flag_inter = [flag_inter;[1;2;2]];
                        elseif j == 3*params.N1/4
                            flag_inter = [flag_inter;[2;2;2]];
                        else
                            flag_inter = [flag_inter;[0;0;0]];
                        end
                    else
                        flag_inter = [flag_inter;[0;0;0]];
                    end
                end
            end
        end
    end
    bonds = [intra_bonds;inter_bonds];
    flag_bonds = [flag_intra;flag_inter];
end

