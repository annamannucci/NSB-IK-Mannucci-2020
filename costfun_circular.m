function [f,gradf] = costfun_circular(r, q, w)
%COSTFUN_CIRCULAR Summary of this function goes here
%   Detailed explanation goes here

    n_agents = length(q)/2;    
    if (nargin < 3)
        w = ones(n_agents, 1); %[n_agents,1]
    end
    
    center = compute_centroid(q);
    f = 0;
    gradf = zeros(size(q));
    for i = 1 : n_agents
        pr = q(2*i-1:2*i)-center;
        gradf(2*i-1:2*i) = gradf(2*i-1:2*i) + w(i)*(r^2-pr'*pr)*pr;
        f = f + w(i)*(r^2-pr'*pr)^2;
    end
    gradf = -4*(1-(1/n_agents))*gradf;
end

