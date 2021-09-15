% function [ output_args ] =Sobol_sensitivity(k, N)
%%%%%k----the number of variables
%%%%%N----the number of QMC samples
    tic
%     k=5;N=256;   
    k=1;N=8;   
   %%%% Plot the convergence of measures
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%Sobol index analysis
    geninit=1;   %%%%%%%%%%%QMC sampling
    xx = []; 
    genfirst = (geninit-1)*N+1;
    genlast = geninit*N;
        if (size(xx,1) < genlast)
            dum = zeros(genlast-size(xx,1)-1,k);
            for i = (size(xx,1)+1):genlast
                dum(i-size(xx,1),:) = LPTAU51(i,k);
            end    
            xx = [xx; dum]; 
            dum = [];
        end
    X = xx(genfirst:genlast,1:k);%%%%%Samples of uniform(0 1)
    %%% T=unifrnd(0,1,N,2*k);%%%%% MC sampling
%     T=norminv(T,0,1); %%%%%Normal variables
%     
%    %PREPARATION OF THE RADIAL SAMPLE MATRIX X(k+2,k) 
%    for i=1:size(X,1)
%        y(i,:)=user_defined(X(i,:)');%%%%get all the function value,size=(2k+2)*N
%    end
%    for i=1:9
%        regression(X,y(:,1));
%        load HDMR.mat;
%        S(i,:)=Si;
%        ST(i,:)=Si_tot;
%    end
%    subplot(1,2,1)   
%    plot(S)
%    legend('S_1','S_2','S_3','S_4','S_5')
%    subplot(1,2,2)
%    plot(ST)
%    legend('S^{tot}_1','S^{tot}_2','S^{tot}_3','S^{tot}_4','S^{tot}_5')
%    toc  
%    
% % end
