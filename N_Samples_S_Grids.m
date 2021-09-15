% Test of different N_Samples and S_Grids
clear all


P.A=[-1 0 ; 0 -1 ; 0 1 ; 1 1];
P.b=[0 ; 0 ; 1 ; 2];

size_q = size(P.A,2);

%rng(0)
%nSamples=1e3;
%nGrid;

% Set up preliminaries.
%sig = 5;
% We will generate 5000 iterations of the chain.
n = 5000; % here we define "n =nSamples" is a parameter of N_sample to optimize!!!!
nSamples = n;
nGrid=1e1;

numchain =[1 2 3 4 5 6 7 8 9 10]; %[1 2 3]; %[2 5 20 100]; % here we define "numchain =nGrid" is another parameter of number of grid to optimize!!!!
%rng(0)
 
% Set up the vectors to store the samples.
% This is 4 chains, 5000 samples.
%X = zeros(length(numchain),n);

q_mean= zeros(length(numchain),n,size_q);

% This is 3 sequences (rows) of summaries.
nu = zeros(length(numchain),n,size_q);		%Note: difference between q_mean and nu

rhat = zeros(1,n,size_q);	%?????	% Track the rhat for each iteration.
% Get the starting values for the chain.
% Use over-dispersed starting points.
% X(1,1) = -10;
% X(2,1) = 10;
% X(3,1) = -5;
% X(4,1) = 5;
% Run the chain.


for j = 2:n	 % j=nSamples
    
    for i = 1:length(numchain)  % numchain=nGrid

%       % Generate variate from proposal distribution.
%       y = randn(1)*sig + X(i,j-1);
%       % Generate variate from uniform.
%       u = rand(1);
%       % Calculate alpha.
%       alpha = normpdf(y,0,1)/normpdf(X(i,j-1),0,1);
%       if u <= alpha
%          % Then set the chain to the y.
%          X(i,j) = y;
%       else
%          X(i,j) = X(i,j-1);
%          
%       end
%       
rng(i)
% nGrid=numchain(i);
%       [q, w_q, w, Q_min, Q_max]=DISCOPOLIS_Cplex(P.A,P.b,j,nGrid);
      [q, w_q, w, Q_min, Q_max]=DISCOPOLIS_Cplex(P.A,P.b,j,nGrid);
      
      
      q_mean(i,j,:)=sum(w_q,2);
     

   end
   
% Get the scalar summary - means of each row.
   nu(:,j,:) = mean(q_mean(:,1:j,:),2);
   rhat(:,j,:) = csgelrub(nu(:,1:j,:));   %?????????? 
 
   min_rhat(j)=min(rhat(:,j,:),[],3);
   max_rhat(j)=max(rhat(:,j,:),[],3);
   
    if max_rhat(j) < 1.2 & j >=50  % diagonize if it converges
    return
    end
%    
%   if max_rhat(j) < 1.2
%       break
%      if max_rhat(j) > 1.2
%       continue
%       end
%    end
      
end


text = 'actual number of samples=';   
n_sample = j;
fprintf('%s %d .\n',text,n_sample);

figure; 
% subplot(2,2,1)
plot(q(1,:),q(2,:),'.m'); 
% xlim([-1 3])
% ylim([-1 2])
hold on; 
plot(q_mean(1,j,1),q_mean(1,j,2),'xk','LineWidth',2,'MarkerSize',15); 
hold on; 
plot(q_mean(2,j,1),q_mean(2,j,2),'b*','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(3,j,1),q_mean(3,j,2),'r*','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(4,j,1),q_mean(4,j,2),'y*','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(5,j,1),q_mean(5,j,2),'k*','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(6,j,1),q_mean(6,j,2),'w*','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(7,j,1),q_mean(7,j,2),'bk','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(8,j,1),q_mean(8,j,2),'c*','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(9,j,1),q_mean(9,j,2),'rk','LineWidth',2,'MarkerSize',10); 
hold on;
plot(q_mean(10,j,1),q_mean(10,j,2),'b*','LineWidth',2,'MarkerSize',10); 

xlabel('v_1'); ylabel('v_2')
 xlim([-0.01 2.01])
ylim([-0.01 1.01])
print('toy_example_Gird=10+nSample5000','-dpdf')

figure; 
plot(max_rhat(:,:),'r.')
hold on
plot(min_rhat(:,:),'g.')
title('monitor convergence')
xlabel('nSamples'); ylabel('min and max R-hat')
xlim([0 5000])
hold off
print('toy_example_Gird=10+nSample5000_convergetrend','-dpdf')

% mu_rhat = mean(rhat,3);
% var_rhat = var(rhat,0,3);
% 
% min_rhat=min(rhat,[],3);
% max_rhat=max(rhat,[],3);

% figure; 
% plot(rhat(:,:,1))
% figure; 
% plot(var_rhat(:,:))

figure; 
subplot(3,2,1)
% plot(nu(1,:,1))
plot(q_mean(1,1:j,1))
title('rng(1)')
xlabel('nSamples'); ylabel('v_1')

subplot(3,2,2)
% plot(nu(1,:,2))
plot(q_mean(1,1:j,2))
%plot(nu(2,:))
title('rng(1)')
xlabel('nSamples'); ylabel('v_2')

subplot(3,2,3)
% plot(nu(1,:,1))
plot(q_mean(2,1:j,1))
title('rng(2)')
xlabel('nSamples'); ylabel('v_1')

subplot(3,2,4)
% plot(nu(1,:,2))
plot(q_mean(2,1:j,2))
title('rng(2)')
xlabel('nSamples'); ylabel('v_2')

subplot(3,2,5)
% plot(nu(1,:,1))
plot(q_mean(3,1:j,1))
title('rng3)')
xlabel('nSamples'); ylabel('v_1')

subplot(3,2,6)
% plot(nu(1,:,2))
plot(q_mean(3,1:j,2))
title('rng(3)')
xlabel('nSamples'); ylabel('v_2')
print('toy_example_Gird=10+nSample5000-c_converge','-dpdf')


%plot(rhat)

% percentage of samples whose sum of weights represents 99.9% of the total
% sum of weights
% sum_weight=sum(w_q,2)
