% Test sperimentally the correctness of iterative formulation
% with two sets of datas: 
% 1. J1 f.r.r. and J2 f.r.r. (1.a J2 l.i. to J1, 1.b J2 not necessarily
%    l.i. to J1)
%       J1 active
%       J2 variable activation

% 2. J1 f.r.r. and J2 f.r.r.
%       J1 active
%       J2 variable activation

%% TEST #1: Define two f.r.r. Jacobians such that the first is l.i. to the second one
clear all; close all; clc;
n_trials = 1000;
n = 6;

for i = 1 : n_trials
    m = randi([3 n],1,1);
    m1 = randi([1 m-1],1,1);
    m2 = m - m1;
    
    %Find m1 l.i. rows
    J1tmp = rand(1,n);
    while size(J1tmp,1) < m1
        J1_ = rand(1,n);
        if rank([J1tmp; J1_])  > rank(J1tmp)
            J1tmp = [J1tmp; J1_];
        end
    end   
    J1{i} = J1tmp;
    
    %Find m2 rows l.i. with J1 and with each other
    J2tmp = [];
    while isempty(J2tmp) || size(J2tmp,1) < m2
        J2_ = rand(1,n);
        if rank([J1tmp; [J2tmp; J2_]])  > rank([J1tmp; J2tmp])
            J2tmp = [J2tmp; J2_];
        end
    end    
    J2{i} = J2tmp;
    li2{i} = 1:size(J2tmp,1);
end

%clear garbage
clear J1tmp J2tmp m1 m2 m J1_ J2_ i

%Control your script
for i = 1 : n_trials
    J1_ = J1{i};
    J2_ = J2{i};
    if (size(J1_,1) > rank(J1_))
        disp('Error: J1 not f.r.r')
    end
    
    if (size(J2_,1) > rank(J2_))
        disp('Error: J2 not f.r.r.')
    end
    
    if (rank([J1_;J2_]) < size([J1_;J2_],1))
        disp('Error: J2 not l.i. to J1')
    end
end

%clear garbage
clear J1_ J2_ i

return

%% TEST #2: Define two f.r.r. Jacobians removing the hypothesis of J1 l.i. to J2
clear all; clc;
n_trials = 1000;
n = 6;

for i = 1 : n_trials
    m = randi([3 n],1,1);
    m1 = randi([1 m-1],1,1);
    m2 = m - m1;
    
    %Find m1 l.i. rows
    J1tmp = rand(1,n);
    while size(J1tmp,1) < m1
        J1_ = rand(1,n);
        if rank([J1tmp; J1_])  > rank(J1tmp)
            J1tmp = [J1tmp; J1_];
        end
    end
    J1{i} = J1tmp;
    
    %%% Find m2 rows l.i. with each other %%%
    li2_ = [];              %indices of rows l.i. to J1
%     J2tmp = rand(1,n);
%     if rank([J2tmp; J1tmp])  > rank(J1tmp)
%         li2_ = 1;
%     end
    J2tmp = J1tmp(1,:); %force to have that
    count = size(J2tmp,1);
    while size(J2tmp,1) < m2
        J2_ = rand(1,n);
        if rank([J2tmp; J2_])  > rank(J2tmp)
            count = count + 1;
            if rank([[J2tmp; J2_]; J1tmp])  > rank([J2tmp; J1tmp])
                li2_ = [li2_; count];
            end
            J2tmp = [J2tmp; J2_];
        end
    end    
    li2{i} = li2_;
    J2{i} = J2tmp;
end


%clear garbage
clear J1tmp J2tmp li2_ m1 m2 m J1_ J2_ i

%Control your script
for i = 1 : n_trials
    J1_ = J1{i};
    J2_ = J2{i};
    if (size(J1_,1) > rank(J1_))
        disp('Error: J1 not f.r.r')
    end
    
    if (size(J2_,1) > rank(J2_))
        disp('Error: J2 not f.r.r.')
    end
    
    if isempty(li2{i}) && rank([J1_; J2_]) > rank(J1_)
       disp('Error: wrong set of l.i. rows. It should not be empty.')
    else
        if rank([J1_;J2_(li2{i},:)]) < size([J1_;J2_(li2{i},:)],1)
            disp('Error: wrong set of l.i. rows.')
        end
    end
end

%clear garbage
clear J1_ J2_ i

return

%% Define activation matrices and damping

%set the type of activation for the second Jacobian
%     1 "disactive"
%     2 "transient"
%     3 "active"
%     4 "li_active", just for A2
set_test_type = [2, 2];

kvo = 0;
lambda = 1; 
epsilon = 1;

% task 1 full active, task 2 with random activations
for i = 1 : n_trials
    J{1} = J1{i};
    J{2} = J2{i};

    %no one active
    for ii = 1:2
        m = size(J{ii},1); 
        switch set_test_type(ii)
            case 1 %disactive
                Atmp = zeros(m,1);
            case 2 %transient
                Atmp = randi([0,1],m,1);
                ind_disactive = find(Atmp == 0);
                if ~isempty(ind_disactive)
                    ind_transient = datasample(ind_disactive,randi(length(ind_disactive)));
                    Atmp(ind_transient) = rand(length(ind_transient),1);
                end
            case 3     %all active
                Atmp = ones(m,1);
            case 4  %active only the linearly independent ones
                Atmp = zeros(m,1);
                Atmp(li2{i}) = ones(length(li2{i}),1);
            default    %all active
                Atmp = ones(m,1);
        end
        A{ii} = Atmp;
    end
    
    A1{i} = A{1};
    A2{i} = A{2};
end

%% iCAT conservative

for i = 1 : n_trials
    J1_ = J1{i};
    A1_ = eye(size(J1_,1));
    J2_ = J2{i};
    A2_ = A2{i};
    
    %correct: use just the l.i. raws
    if isempty(li2{i}) 
        J2_li = [];
        Ptilda_c{i} = eye(n);
    else
        J2_li = J2_(li2{i},:);
        Ptilda_c{i} = eye(n)-damped_pinv(J2_li,epsilon,lambda)*diag(A2_(li2{i}))*J2_li;
    end   
    
    %iterative formulation
    JA_k = [];
    PAinv = eye(n);
    
    Jk = J2_;
    H = diag(A2_);
    JA_k = [Jk;JA_k];
    
    JAinv = damped_pinv(JA_k, epsilon, lambda);
    Tk = JAinv(:,1:size(Jk,1));

    %Compute null space and projector
    TkTkinv = zeros(size(Tk,1));
    [u,~,~] = svd(Tk);
    for j = 1 : rank(Tk)
        TkTkinv = TkTkinv + u(:,j)*u(:,j)';
    end

    %iterative law to speed up the computation (if tk
    %computed as rank update)
    Ptilda = (eye(n)-TkTkinv)*PAinv+TkTkinv;
    PAinv = (eye(n)-Tk*Jk)*PAinv+Tk*(eye(size(H,1))-H)*Jk;
    
    Jk = J1_;
    H = diag(A1_);
    JA_k = [Jk;JA_k];
    
    JAinv = pinv(JA_k);
    Tk = JAinv(:,1:size(Jk,1));

    %Compute null space and projector
    TkTkinv = zeros(size(Tk,1));
    [u,~,~] = svd(Tk);
    for j = 1 : rank(Tk)
        TkTkinv = TkTkinv + u(:,j)*u(:,j)';
    end

    Ptilda = (eye(n)-TkTkinv)*PAinv+TkTkinv;
    PAinv = (eye(n)-Tk*Jk)*PAinv+Tk*(eye(size(H,1))-H)*Jk;
    
    Ptilda_i{i} = Ptilda;
    
    %control rank
    rank_ec{i} = abs(rank(J1_*Ptilda_c{i})-rank(J1_));
    rank_ei{i} = abs(rank(J1_*Ptilda_i{i})-rank(J1_));
    
    %Frobenius norm
    fro{i} = abs(norm(Ptilda_c{i},'fro')-norm(Ptilda_i{i},'fro'));
    
    %Direct comparison
    spa{i} = norm(Ptilda_c{i} - Ptilda_i{i});
    %spa{i} = norm(J1_*Ptilda_c{i} - J1_*Ptilda_i{i});
end

%anche con damping, il rango resta costante.
rank_ec_ = cell2mat(rank_ec);
rank_ei_ = cell2mat(rank_ei);
figure('Name', 'Rank invariance')
plot(rank_ec_)
hold on
plot(rank_ei_)

fro_ = cell2mat(fro);
figure('Name', 'Error in frobenius norm')
plot(fro_)

max(fro_)
std(fro_)
mean(fro_)

spa_ = cell2mat(spa);
figure('Name', 'Error in span')
plot(spa_)

max(spa_)
std(spa_)
mean(spa_)

clear fro_ rank_ec_ rank_ei_ Jk Ptilda Tk TkTkinv u H JA_k JAinv Painv m2 A2tmp i j J1_ J2_ J2li 

%% iCAT RP

%ATTENZIONE Senza damping l'errore commesso tra calcolo corretto e formulazioni iterative è praticamente nullo,
%con damping il calcolo corretto e quello dampato differiscono.
%La frequenza con cui si ha quest'errore è funzione degli autovalori di J'AJ
%quindi sia dell'attivazione che degli autovalori di J.

for i = 1 : n_trials
    
    J1_ = J1{i};
    A1_ = A1{i};
    J2_ = J2{i};
    A2_ = A2{i};
    up_ = [];
    
    %correct: use just the l.i. raws
    if isempty(li2{i}) 
        J2_li = [];
        Ptilda_c{i} = eye(n);
    else
        J2_li = J2_(li2{i},:);
        Ptilda_c{i} = eye(n)-pinvSC(J2_li,diag(A2_(li2{i})),eye(n),kvo,epsilon,lambda)*J2_li;
    end   
    
    %iterative formulation
    % I - pinvSC(J_{i}) * J_{i} = (I - TkTkinv) * (I - pinvSC(J_{i+1}) *J_{i+1}) + TkTkinv
    JA_k = [];
    Hbuff = [];
    
    Jk = J2_;
    H = diag(A2_);
    PAinvc = eye(n);
    PAinv = eye(n);
    JA_k = [Jk;JA_k];
    Hbuff = blkdiag(H,Hbuff);
    
    JAinv = pinvSC(JA_k,Hbuff,eye(n),kvo,epsilon,lambda);
    JAinv_ = JAinv*pinv(Hbuff);
    Tk = JAinv_(:,1:size(Jk,1));

    %Compute null space and projector
    TkTkinv = zeros(size(Tk,1));
    [u,~,~] = svd(Tk);
    for j = 1 : rank(Tk)
        TkTkinv = TkTkinv + u(:,j)*u(:,j)';
    end

    %iterative law to speed up the computation (if tk
    %computed as rank update)
    %Ptilda = (eye(n)-TkTkinv)*(eye(n)-JAinv*JA_k)+TkTkinv;
    Ptilda = (eye(n)-TkTkinv)*PAinvc+TkTkinv;
    PAinvc = (eye(n)-Tk*Jk)*PAinv+Tk*(eye(size(H,1))-H)*Jk;
    PAinv = eye(n)-JAinv*JA_k;
    up_ = [up_, norm(PAinvc-PAinv)];
    
    Jk = J1_;
    H = diag(A1_);
    JA_k = [Jk;JA_k];
    Hbuff = blkdiag(H,Hbuff);
    
    JAinv = pinvSC(JA_k,Hbuff,eye(n),kvo,epsilon,lambda);
    JAinv_ = JAinv*pinv(Hbuff);
    Tk = JAinv_(:,1:size(Jk,1));

    %Compute null space and projector
    TkTkinv = zeros(size(Tk,1));
    [u,~,~] = svd(Tk);
    for j = 1 : rank(Tk)
        TkTkinv = TkTkinv + u(:,j)*u(:,j)';
    end

    %computed as rank update)
    %Ptilda = (eye(n)-TkTkinv)*(eye(n)-JAinv*JA_k)+TkTkinv;
    Ptilda = (eye(n)-TkTkinv)*PAinvc+TkTkinv;
    PAinvc = (eye(n)-Tk*Jk)*PAinv+Tk*(eye(size(H,1))-H)*Jk;
    PAinv = eye(n)-JAinv*JA_k;
    up_ = [up_, norm(PAinvc-PAinv)];

    Ptilda_i{i} = Ptilda;
    up{i} = up_;
    
    %control rank
    rank_ec{i} = abs(rank(J1_*Ptilda_c{i})-rank(J1_));
    rank_ei{i} = abs(rank(J1_*Ptilda_i{i})-rank(J1_));
    
    %Frobenius norm
    fro{i} = abs((norm(Ptilda_c{i},'fro')-norm(Ptilda_i{i},'fro'))/norm(Ptilda_c{i},'fro'));
    
    %Direct comparison
    spa{i} = norm((Ptilda_c{i}-Ptilda_i{i}))/norm(Ptilda_c{i});
    %spa{i} = norm(J1_*Ptilda_c{i} - J1_*Ptilda_i{i});
end

%sempre nullo, indipendentemente dal valore di lambda usato (quindi
%sembrerebbe sempre verificato
rank_ec_ = cell2mat(rank_ec);
rank_ei_ = cell2mat(rank_ei);
figure('Name', 'Rank invariance')
plot(rank_ec_)
hold on
plot(rank_ei_)

fro_ = cell2mat(fro);
figure('Name', 'Error in frobenius norm')
plot(fro_)

max(fro_)
std(fro_)
mean(fro_)

spa_ = cell2mat(spa);
figure('Name', 'Error in span')
plot(spa_)

max(spa_)
std(spa_)
mean(spa_)

clear up_
for i = 1 : n_trials
    up_(i,:) = up{i};
end
figure('Name', 'Error in updating')
plot(up_(:,1))
hold on
plot(up_(:,2))

clear fro_ rank_ec_ rank_ei_ Jk Ptilda Tk TkTkinv u H JA_k JAinv Painv m2 A2tmp i j J1_ J2_ J2li 