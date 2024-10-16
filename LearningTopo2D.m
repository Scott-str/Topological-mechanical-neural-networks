clear;clc;

%% data
load fisheriris;
X = meas;
y = species;
X(1:50,:) = [];y(1:50) = [];
X = 1./(max(X)-min(X)).*X-1./(max(X)-min(X)).*max(X)+1;
[unique_labels, ~, label_indices] = unique(y);
num_labels = numel(unique_labels);
one_hot_labels = eye(num_labels);
one_hot_labels = one_hot_labels(label_indices,:);
split_ratio = 0.7;
num_train_examples = round(size(X, 1) * split_ratio);
rand_indices = randperm(size(X, 1));
X_train = X(rand_indices(1:num_train_examples), :);
y_train = one_hot_labels(rand_indices(1:num_train_examples), :);
X_test = X(rand_indices(num_train_examples+1:end), :);
y_test = one_hot_labels(rand_indices(num_train_examples+1:end), :);

figure;
scatter(X_train(:,1),X_train(:,3),20,y_train(:,1),'filled');
axis equal;

%% parameters
params = struct('a',1,'N1',4*6,'N2',4*4,'m',1,'omega',370,'beta',5);
% visualization
[pos,bonds,flag_bonds] = U_shape(params);
[flag_prune_bonds,flag_node] = prune(pos,bonds,params);
[flag_prune_bonds_Bott,flag_node_Bott] = pruneBott(pos,bonds,params);

params.pos = pos;
params.bonds = bonds;
params.flag_bonds = flag_bonds;

% input
params.ind_input = [2+6*(3*params.N2/4-1)*params.N1+6*(params.N1/2-1),...
    5+6*(3*params.N2/4-1+1)*params.N1+6*(params.N1/2-1)];
% fixed
params.ind_fix = find(flag_node == 2);
% free
params.ind_free = find(flag_node == 0);
% two output ports
row_probe_1 = [2:3*params.N2/4];col_probe_1 = [params.N1/4,params.N1/4+1];
ind_output_1 = [];
for i = 1:length(row_probe_1)
    for j = 1:length(col_probe_1)
        ind_output_1 = [ind_output_1,[1:6]+6*(row_probe_1(i)-1)*params.N1+6*(col_probe_1(j)-1)];
    end
end
params.ind_output_1 = unique(ind_output_1);

row_probe_2 = [2:3*params.N2/4];col_probe_2 = [3*params.N1/4,3*params.N1/4+1];
ind_output_2 = [];
for i = 1:length(row_probe_2)
    for j = 1:length(col_probe_2)
        ind_output_2 = [ind_output_2,[1:6]+6*(row_probe_2(i)-1)*params.N1+6*(col_probe_2(j)-1)];
    end
end
params.ind_output_2 = unique(ind_output_2);

optims = struct('alpha',5e-2,'beta1',0.9,'beta2',0.999,...
    'batch_size',16,'epochs',30);

figure('Name','Topological Metamaterials','color','k');
for i = 1:size(params.bonds,1)
    if params.flag_bonds(i) == 1
        plot([pos(params.bonds(i,1),1),pos(params.bonds(i,2),1)],[pos(params.bonds(i,1),2),pos(params.bonds(i,2),2)],'r-');
    elseif params.flag_bonds(i) == 2
        plot([pos(params.bonds(i,1),1),pos(params.bonds(i,2),1)],[pos(params.bonds(i,1),2),pos(params.bonds(i,2),2)],'c-');
    else
        plot([pos(params.bonds(i,1),1),pos(params.bonds(i,2),1)],[pos(params.bonds(i,1),2),pos(params.bonds(i,2),2)],'y-');
    end
    hold on;
end
axis equal;axis off;

figure('Name','Topological Metamaterials','color','k');
bonds_left = bonds(flag_prune_bonds == 0,:);
flag_bonds_left = flag_bonds(flag_prune_bonds == 0);
for i = 1:size(bonds_left,1)
    if flag_bonds_left(i) == 1
        plot([pos(bonds_left(i,1),1),pos(bonds_left(i,2),1)],[pos(bonds_left(i,1),2),pos(bonds_left(i,2),2)],'r-');
    elseif flag_bonds_left(i) == 2
        plot([pos(bonds_left(i,1),1),pos(bonds_left(i,2),1)],[pos(bonds_left(i,1),2),pos(bonds_left(i,2),2)],'c-');
    else
        plot([pos(bonds_left(i,1),1),pos(bonds_left(i,2),1)],[pos(bonds_left(i,1),2),pos(bonds_left(i,2),2)],'y-');
    end
    hold on;
end
% plot(pos(:,1),pos(:,2),'k.','MarkerSize',5);
plot(pos(params.ind_fix,1),pos(params.ind_fix,2),'g.','MarkerSize',5);
plot(pos(params.ind_input,1),pos(params.ind_input,2),'wx','MarkerSize',5);
plot(pos(params.ind_output_1,1),pos(params.ind_output_1,2),'yx','MarkerSize',5);
plot(pos(params.ind_output_2,1),pos(params.ind_output_2,2),'cx','MarkerSize',5);
% plot(pos(flag_node == 0,1),pos(flag_node == 0,2),'k.','MarkerSize',2);
% plot(pos(flag_node == 2,1),pos(flag_node == 2,2),'cx','MarkerSize',2);
axis equal;axis off;

params.bonds_new = params.bonds(flag_prune_bonds == 0,:);
params.flag_bonds_new = params.flag_bonds(flag_prune_bonds == 0);

figure('Name','NonTrivial Phase','color','k');
bonds_left = bonds(flag_prune_bonds_Bott == 0,:);
flag_bonds_left = flag_bonds(flag_prune_bonds_Bott == 0);
for i = 1:size(bonds_left,1)
    if flag_bonds_left(i) == 1
        plot([pos(bonds_left(i,1),1),pos(bonds_left(i,2),1)],[pos(bonds_left(i,1),2),pos(bonds_left(i,2),2)],'r-');
    elseif flag_bonds_left(i) == 2
        plot([pos(bonds_left(i,1),1),pos(bonds_left(i,2),1)],[pos(bonds_left(i,1),2),pos(bonds_left(i,2),2)],'c-');
    else
        plot([pos(bonds_left(i,1),1),pos(bonds_left(i,2),1)],[pos(bonds_left(i,1),2),pos(bonds_left(i,2),2)],'y-');
    end
    hold on;
end
plot(pos(flag_node_Bott == 2,1),pos(flag_node_Bott == 2,2),'g.','MarkerSize',5);
axis equal;axis off;

%%
params.C = compatibility_matrix(params.pos,params.bonds);

%%
% Bound
ind_softer = find(params.flag_bonds_new == 0);
ind_stiffer = find(params.flag_bonds_new == 1);
ind_interface = find(params.flag_bonds_new == 2);
LB1 = 0.75e5*ones(length(ind_softer),1);
UB1 = 0.85e5*ones(length(ind_softer),1);
LB2 = 1.15e5*ones(length(ind_stiffer),1);
UB2 = 1.25e5*ones(length(ind_stiffer),1);
LB3 = 1e5*ones(length(ind_interface),1);
UB3 = 1e5*ones(length(ind_interface),1);
params.LB(ind_softer) = LB1;
params.LB(ind_stiffer) = LB2;
params.LB(ind_interface) = LB3;
params.UB(ind_softer) = UB1;
params.UB(ind_stiffer) = UB2;
params.UB(ind_interface) = UB3;
k0_vec = zeros(size(params.bonds_new,1),1);
m0_vec = params.m*ones(2*6*params.N1*params.N2,1);
params.M = diag(m0_vec);
k0_vec(ind_softer) = 0.8e5;
k0_vec(ind_stiffer) = 1.2e5;
k0_vec(ind_interface) = 1e5;
k0_vec = InverseConvert_k(k0_vec,params);

params.flag_prune_bonds = flag_prune_bonds;
params.flag_node = flag_node;
params.flag_prune_bonds_Bott = flag_prune_bonds_Bott;
params.flag_node_Bott = flag_node_Bott;

%% optimization
results.kHist = []; % optimized k
results.fHist = []; % cost function
results.pred_train = []; % prediction of train
results.pred_test = []; % prediction of test
results.acc_train = []; % accuracy of training data
results.acc_test = []; % accuracy of testing data

% Initialize moment estimates
m_moment = zeros(size(params.bonds_new,1), 1);
v_moment = zeros(size(params.bonds_new,1), 1);

iter = 1;
for epoch = 1:optims.epochs
    ind_shuffle = randperm(size(X_train,1));
    X_train_shuffle = X_train(ind_shuffle,:);
    y_train_shuffle = y_train(ind_shuffle,:);
    count = 1;
    for i = 1:optims.batch_size:size(X_train_shuffle,1)
        idx = i:min(i+optims.batch_size-1,size(X_train_shuffle,1));
        X_batch = X_train_shuffle(idx,:);
        y_batch = y_train_shuffle(idx,:);
        % Get gradients
        [results.fHist{epoch}(count),sgCurr] = Topo2Dbatch(k0_vec,X_batch,y_batch,params,optims,'1');
        results.kHist{epoch}(:,count) = Convert_k(k0_vec,params);
        % Adam optimizer
        [k0_vec,m_moment,v_moment] = Train_Adam(k0_vec,sgCurr,m_moment,v_moment,iter,params,optims);
        iter = iter + 1;
        count = count + 1;
    end
    [results.pred_train{epoch},results.acc_train(epoch)] = prediction(results.kHist{epoch}(:,end),params,optims,X_train,y_train);
    [results.pred_test{epoch},results.acc_test(epoch)] = prediction(results.kHist{epoch}(:,end),params,optims,X_test,y_test);
    formatSpec = 'Epoch:%d/%d  Cost:%e  Accuracy(Train):%f  Accuracy(Test):%f\n';
    fprintf(formatSpec,epoch,optims.epochs,mean(results.fHist{epoch}),results.acc_train(epoch),results.acc_test(epoch));
end

%% display
figure;
fcost = [];
for i = 1:length(results.fHist)
    fcost(i) = mean(results.fHist{i});
end
subplot(1,2,1);
plot(fcost,'r-');
xlabel('epoch');ylabel('cost function');
subplot(1,2,2);
plot(results.acc_train,'DisplayName','Train');
hold on;
plot(results.acc_test,'DisplayName','Test');
xlabel('epoch');ylabel('accuracy');
legend();

k_cur = results.kHist{end}(:,end);

k_cur_ori = zeros(size(params.bonds,1),1);
k_cur_ori(params.flag_prune_bonds == 0) = k_cur;
figure('Name','Topological metamaterials','Color','k');
for i = 1:size(params.bonds,1)
    if params.flag_prune_bonds(i) == 0
        [c,~] = colorline(0.75e5,1.25e5,k_cur_ori(i),'RdYlBu');
        plot([params.pos(bonds(i,1),1),params.pos(bonds(i,2),1)],[params.pos(bonds(i,1),2),params.pos(bonds(i,2),2)],'-',...
        'linewidth',1,'color',c);
        hold on;
    end
end
plot(params.pos(params.ind_input,1),params.pos(params.ind_input,2),'wx','MarkerSize',5);
plot(params.pos(params.flag_node == 2,1),params.pos(params.flag_node == 2,2),'g.','MarkerSize',5);
axis equal;axis off;

k_ori = zeros(size(params.bonds,1),1);
k_ori(params.flag_prune_bonds == 0) = results.kHist{1}(:,1);
figure('Name','Difference','Color','k');
for i = 1:size(params.bonds,1)
    if params.flag_prune_bonds(i) == 0
        [c,~] = colorline(-0.05e5,0.05e5,k_cur_ori(i)-k_ori(i),'RdYlBu');
        plot([params.pos(bonds(i,1),1),params.pos(bonds(i,2),1)],[params.pos(bonds(i,1),2),params.pos(bonds(i,2),2)],'-',...
        'linewidth',1,'color',c);
        hold on;
    end
end
plot(params.pos(params.flag_node == 2,1),params.pos(params.flag_node == 2,2),'g.','MarkerSize',5);
axis equal;axis off;

k_ori = zeros(size(params.bonds,1),1);
k_ori(params.flag_prune_bonds == 0) = results.kHist{1}(:,1);
k_cur_ori = zeros(size(params.bonds,1),1);
k_cur_ori(params.flag_prune_bonds == 0) = k_cur;
figure('Name','NonTrivial Phase','Color','w');
for i = 1:size(params.bonds,1)
    if params.flag_prune_bonds_Bott(i) == 0
        [c,~] = colorline(-0.05e5,0.05e5,k_cur_ori(i)-k_ori(i),'RdYlBu');
        plot([params.pos(bonds(i,1),1),params.pos(bonds(i,2),1)],[params.pos(bonds(i,1),2),params.pos(bonds(i,2),2)],'-',...
        'linewidth',1,'color',c);
        hold on;
    end
end
plot(params.pos(params.flag_node_Bott == 2,1),params.pos(params.flag_node_Bott == 2,2),'g.','MarkerSize',5);
axis equal;axis off;

%% Time domain
k_cur = results.kHist{end}(:,end);

num = 20;
fps = 1000;
dt = 1/fps;

x0 = zeros(2*2*sum(params.flag_node~=1),1);

options = odeset('RelTol',1e-8,'AbsTol',1e-10);
TIME = 0:dt:1;
[Ts,X] = ode45(@(t,x) linearsolver(t,x,k_cur,X_train(num,:),params),TIME,x0,options);

ux = X(:,1:2:size(X,2)/2);
uy = X(:,2:2:size(X,2)/2);
vx = X(:,size(X,2)/2+1:2:end);
vy = X(:,size(X,2)/2+2:2:end);

inx1 = find(params.flag_node==1);
pos_eff = params.pos;
pos_eff(inx1,:) = [];

figure;
for i = 1:length(TIME)
    scatter(pos_eff(:,1),pos_eff(:,2),20,vx(i,:),'filled');
    clim([-1e-3,1e-3]);
    colormap jet;colorbar;
    axis equal;axis off;
    title(['t=',num2str(TIME(i)),'s']);
    drawnow;
end
