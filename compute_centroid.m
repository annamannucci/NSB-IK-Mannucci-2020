function center = compute_centroid(q)
%CENTER compute the centroid of the formation
    center = zeros(2,1);
    n_agents = length(q)/2;
    for i = 1 : 1 : n_agents
        center = center + q(2*i-1:2*i);
    end
    center = center/n_agents;
end

