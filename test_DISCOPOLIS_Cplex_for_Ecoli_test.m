clear all
changeCobraSolver ('ibm_cplex', 'all');
rng(0)
nSamples=1e3;  % One parameter to tune!!!
nGrid=5e0;  % just 5 grids?? % One parameter to tune!!! 
load('Ecoli_core_model')
model = changeRxnBounds(model,'EX_glc(e)',-18.5,'l');
model.c=0*model.c;  %???

% ATTENTION : bounds modification !!!

model.lb(model.lb==-1e3)=-20*ones;
model.ub(model.ub==1e3)=20*ones;

[A, b, v0, Z] = A_b_from_Cobra_model_Cplex(model);
t0=cputime;
%[q, w_q, w, Q_min, Q_max, errors]=DISCOPOLIS_Cplex(A,b,nSamples,nGrid);
[q, w_q, w, Q_min, Q_max]=DISCOPOLIS_Cplex_test(A,b,nSamples,nGrid);
tf=cputime-t0

%var_q=var(q')';  std(q,0,2), variance of q for one run of nSamples
%We can use bootstrapping k times of nSamples to get variance of q_mean?
%var_v_mean=Z*std(q_mean)^2*Z'


q_mean=sum(w_q,2);  % sum of each row (here, 1000 columns of q based 1000 samples)
v_mean=v0+Z*q_mean;

load v_round_ecoli_ub_20
figure; plot(v_mean,'x'); hold; plot(v_round,'or');

steps_q=1*1e1;  % One parameter to tune!!!  is it a variable impacting the modelling?
mpdf_q=marginal_distributions(q,w,Q_min,Q_max,steps_q);
figure;
for i=1:24
    subplot(4,6,i); histogram('BinEdges',linspace(Q_min(i),Q_max(i),steps_q+1),'BinCounts',mpdf_q(i,:));
    xlabel(strcat('q',num2str(i))); ylabel(strcat('mpdf(q',num2str(i),')'))
    axis([Q_min(i) Q_max(i) 0 max(mpdf_q(i,:))])
end
