function rhat= csgelrub(nu)
% CSGELRUB Gelman and Rubin convergence diagnostic.
%
%   RHAT = CSGELRUB(NU) Returns the Gelman and Rubin convergence 
%   diagnostic for Markov chain Monte Carlo. The input to the
%   function is a mtarix of scalar summaries NU. Each row of the
%   matrix is a sequence of scalar summaries from the chain.
%   For example, the scalar summaries might be the mean, median, etc.


[k,n,size_q]=size(nu);
% First find the Between Sequence variance
museq = mean(nu,2);	% mean of each row
mutot = mean(nu(:,:,:),1:2);
% square each element and sum, between-sequence variance
B = n/(k-1)*sum((museq - mutot).^2);

% Find the within-sequence variance
varseq = var(nu,0,2);	% variance of each row
% varseq_mean = mean(varseq,3);	% mean of flux variances

W = mean(varseq,1); % within-sequence variance of each flux

varhat = (n-1)/n*W + B/n;

rhat = varhat./W;

