function [mpdf_q]=marginal_distributions(q,w,Q_min,Q_max,steps_q)
w_tot=sum(w);
mpdf_q=zeros(size(q,1),steps_q);
for i=1:size(q,1)
    int_q=(Q_max(i)-Q_min(i))/steps_q;
    for j=1:steps_q
        I_int_q=find(q(i,:)>=Q_min(i)+(j-1)*int_q & q(i,:)<Q_min(i)+j*int_q);
        mpdf_q(i,j)=sum(w(I_int_q)/w_tot);
    end
end