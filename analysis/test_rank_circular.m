clear all
n_agents = 4;
q = sym('q',[2*n_agents,1],'real');
center = sym(zeros(2,1));
for i = 1 : n_agents
    center(1) = center(1)+q(2*i-1);
    center(2) = center(2)+q(2*i);
end
center = 1/n_agents*center;
pr = sym(zeros(2,1));

J6 = zeros(2,2*n_agents);
for i = 1 : n_agents
    J6(:,(i-1)*2+1:2*i) = 1/n_agents*eye(2);                   
end

%Non garantisce errore nullo neanche se fatto da solo
% for i = 1 : n_agents 
%     pr(1:2) = q(2*i-1:2*i)-qc(1:2);
%     J9(i, 2*i-1:2*i) = pr';
%     xdot9(i) = pr(1:2)'*qcdot(1:2) + .4 * (r^2/2 -.5*(pr(1:2))'*(pr(1:2)));
% end
    for i = 1 : n_agents
        pr = q(2*i-1:2*i)-center;
        for j = 1 : n_agents
            if (j ~= i)
                J9(i, 2*j-1:2*j) = (-1/n_agents)*pr';
            else
                J9(i, 2*j-1:2*j) = (1-1/n_agents)*pr';
            end
        end
    end
    
%Task #10: PERIMETER
for i = 1 : n_agents
    if i > 1 && i < n_agents
        J10(1,2*i-1:2*i) = (q(2*i-1:2*i)-q(2*(i-2)+1:2*(i-1)))'+(q(2*i-1:2*i)-q(2*i+1:2*i+2))';
    else
        if (i == 1)
            J10(1,2*i-1:2*i) = (q(2*i-1:2*i)-q(2*(n_agents-1)+1:2*(n_agents-1)+2))'+(q(2*i-1:2*i)-q(2*i+1:2*i+2))';
        else
            J10(1,2*i-1:2*i) = (q(2*i-1:2*i)-q(2*(i-2)+1:2*(i-1)))'+(q(2*i-1:2*i)-q(1:2))';
        end
    end
end

J9fun = matlabFunction(J9,'var',{q});
J10fun = matlabFunction(J10,'var',{q});
centerFun = matlabFunction(center,'var',{q});
q = zeros(2*n_agents,1);
for i = 1 : n_agents
    q(2*i-1) = cos(2*(i-1)*pi/n_agents);
    q(2*i) = sin(2*(i-1)*pi/n_agents);
end
rank(J9fun(q))
rank([J9fun(q);J10fun(q)])
svd([J9fun(q);J10fun(q)])
J10sym = J10fun(q);
J9sym = J9fun(q);
rank(null(J6))
rank(null(J9sym))
rank(null([J6;J9sym]))
M33 = J10sym*(eye(2*n_agents)-pinv(J6)*J6)*(eye(2*n_agents)-pinv(J9sym)*J9sym)*pinv(J10sym)
nll=null([J6;J9sym]);
rank(nll)
oth=orth((eye(2*n_agents)-pinv(J6)*J6)*(eye(2*n_agents)-pinv(J9sym)*J9sym));
rank(oth)
rank([nll,oth])

q = rand(2*n_agents,1);
rank(J9fun(q))
rank([J9fun(q);J10fun(q)])
J10sym = J10fun(q);
J9sym = J9fun(q);
M33 = J10sym*(eye(2*n_agents)-pinv(J6)*J6)*(eye(2*n_agents)-pinv(J9sym)*J9sym)*pinv(J10sym)
nll=null([J6;J9sym]);
rank(nll)
oth=orth((eye(2*n_agents)-pinv(J6)*J6)*(eye(2*n_agents)-pinv(J9sym)*J9sym));
rank(oth)
rank([nll,oth])
%%
clear all
n_agents = 4;
q = sym('q',[n_agents,1],'real');
center = sym(0);
for i = 1 : n_agents
    center(1) = center(1)+q(i);
end
center = 1/n_agents*center;
pr = sym(0);
%Non garantisce errore nullo neanche se fatto da solo
% for i = 1 : n_agents 
%     pr(1:2) = q(2*i-1:2*i)-qc(1:2);
%     J9(i, 2*i-1:2*i) = pr';
%     xdot9(i) = pr(1:2)'*qcdot(1:2) + .4 * (r^2/2 -.5*(pr(1:2))'*(pr(1:2)));
% end
    for i = 1 : n_agents
        pr = q(i)-center;
        for j = 1 : n_agents
            if (j ~= i)
                J9(i,j) = (-1/n_agents)*pr';
            else
                J9(i,j) = (1-1/n_agents)*pr';
            end
        end
    end
    
%Task #10: PERIMETER
for i = 1 : n_agents
    if i > 1 && i < n_agents
        J10(1,i) = (q(i)-q(i-1))'+(q(i)-q(i+1))';
    else
        if (i == 1)
            J10(1,i) = (q(1)-q(n_agents))'+(q(1)-q(2))';
        else
            J10(1,i) = (q(n_agents)-q(n_agents-1))'+(q(n_agents)-q(1))';
        end
    end
end
