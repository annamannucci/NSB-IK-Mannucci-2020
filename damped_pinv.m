function [ Minv, damped ] = damped_pinv(M, lambdaq_max, epsilon)
%DUMPED_PINV Compute dumped pseudo-inverse of matrix M
%   Parameters:
%   - lambdaq_max := dumping factor
%   - epsilon := amplitude of the singular value region
%   Output:
%   Minv dumped pseudo-inverse
%   dumped flag, 1 if dumped.
% keyboard
damped = 0;
% [U,S,V] = svd(M'*M);
% 
% epsi = epsilon;
% epsilon = epsilon^2;
% lambdaq_max = 1e-1;
% lambdaq = 0;
% 
% P = zeros(size(S));
% r = rank(S,1e-12);
% 
% if r > 0
%     if  S(r,r) < epsilon
%         lambdaq = (1-(S(r,r)/epsilon))*lambdaq_max;
%         dumped = 1;
%     end
%     %stdev = epsilon/4;
%     for i = 1 : r
%         if S(i,i) < epsilon
%             P(i,i) = lambdaq;%_max*exp(-1/2*(S(i,i)/stdev)^2);
%         end
%     end
%     inv = pinv(M'*M+V*P*V');
% else
%     inv = zeros(size(S));
% end
% 
% Minv = inv*M';

%%
% epsilon = epsi;

if epsilon == 0
    Minv = pinv(M);
    damped = 0;
else
    lambdaq = 0;
    [U,Sig,V] = svd(M);
    P = zeros(size(Sig,1));
    r = rank(Sig,1e-12);
    %stdev = epsilon/4;

    A = zeros(size(M'));

    if r > 0
        %lambdaq = lambdaq_max*exp(-1/2*(Sig(r,r)/stdev)^2);%(1-(Sig(r,r)/epsilon)^2)*lambdaq_max;
        if Sig(r,r) < epsilon
            lambdaq = (1-(Sig(r,r)/epsilon)^2)*lambdaq_max;
            %lambdaq=lambdaq_max;
            damped = 1;
        end

        for i = 1 : r
            %if Sig(i,i) < epsilon
                P(i,i) = lambdaq;%_max*exp(-1/2*(S(i,i)/stdev)^2);
            %end
        end

        A = zeros(size(M'));
        for i = 1 : r
            v = V(:,i);
            u = U(:,i);
            A = A + Sig(i,i)/(Sig(i,i)^2 + P(i,i)) * v * u';
        end
    end

    Minv = A;
end

