function [rand2] = Sobol_QMC(nSamples)  %nSamples=1000


rng default  % For reproducibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sobol quasirandom point set
p = sobolset(nSamples,'Skip',1e3,'Leap',1e2);

p = scramble(p,'MatousekAffineOwen'); %Apply a random linear scramble combined with a random digital shift by using scramble.

% Use net to generate the first 1000 points.
rand2 = net(p,1);  % This is equivalent to X0 = p(1:500,:);