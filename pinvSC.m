function [ Xinv, damp, smin ] = pinvSC(X, A, Q, k_mcv, epsilon_svo, lambda_svo)
%PINVSC compute the smoothed pseudo-inverse for swithching activation tasks
%   Implementation of the new pseudo-inverse operator introfuced by
%Simetti Casalino in the paper "A Novel Practical Technique To Integrate
%Inequality Control Objectives and Task In Transitions in Priority Based
%Control."

n = size(X,2);

damp = 0;

B = X'*A*X+k_mcv*(eye(n)-Q)'*(eye(n)-Q);
[~,S,V] = svd(B);

P = zeros(size(S));

r = rank(S);
smin = 0;

if r > 0
   smin = S(r,r);
   lambdaq = 0;
    if sqrt(S(r,r)) < epsilon_svo
        lambdaq = lambda_svo*(1-(sqrt(S(r,r))/epsilon_svo)^2);
        damp = 1;
    end
    for ii=1:rank(S)   
        if S(ii,ii) < epsilon_svo
             P(ii,ii) = lambdaq;
        end
    end
%    %Gaussian to preserve continuity 
%     stdev = epsilon/4;
%     if S(r,r) < epsilon
%         damp = 1;
%     end
%     for i = 1 : r
%         if (S(i,i) < epsilon)
%             P(i,i) = lambda_svo*exp(-1/2*(S(i,i)/stdev)^2);
%         end
%     end
    inv = pinv(B + V*P*V');
else
    inv = zeros(size(S));
end

Xinv = inv*X'*A*A;

end

