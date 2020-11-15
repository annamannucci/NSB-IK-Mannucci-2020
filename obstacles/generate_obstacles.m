clear all;

trials = 10;
n_obs = 7;

pobs = [
      140 2.5;
      145 -5;
      150 5;
      167.5 2.5;
      167.5 -2.5;
      170 5;
      170 -5
        ];
save("obs0.mat", 'pobs');

for i = 1 : trials
    x = randi([140 170],1,n_obs);
    y = randi([-7 7],1,n_obs);
    pobs = zeros(n_obs,2);
    for ii = 1 : n_obs
        pobs(ii,:) = [x(ii) y(ii)];
    end

    save("obs"+i+".mat", 'pobs');
end

clear x y trials n_obs
