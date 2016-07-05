function [L_new,R_new,O_new,S_new,A_new]=bar_update(R,L,O,S,A,p,E_L,E_R,dt)
% the T,H,I is the parameter for all the cells(vector),which is embeded in
% the data structure of L, R, O, S, A
% T,H,I each has two column,T,H first column is vertical, second column is
% horizontal; for I, first column is from vertical competition, second column is from
% horizontal competition
% L is the position of the neuron, which is a unit vector(in the vectorized function, it should be a identity matrix)
n=size(L.T,1);
% n is the total number of the neurons
% tau, g , sigma are the parameters, E_R E_L is the input
%% compute the activation term P_R, P_L , which is the normalized true input of the neuron, this complicated implementation is for traveling in ring
% using the circulant matrix for simulation
cir_mat=repmat(1:n,n,1);
% x is the relfative distance of the neuron to others
x=abs(cir_mat'-cir_mat);
% P_R,P_L should truned out to be two column matrix, the first column is
% the vertical, the second column is horizontal
% g0 is the new factor for the
I_L=[sum(O.I_L,2) sum(O.I_L,2)];
I_R=[sum(O.I_R,2) sum(O.I_R,2)];
% p.omiga=1 is the parameter for the attention factor
P_R=max(E_R-p.gI*exp(-x.^5/p.sigma^5)*I_R+p.g1*(R.T'*exp(-x.^5/(2*p.sigma)^5))'+p.g0*(L.T'*exp(-x.^5/(2*p.sigma)^5))' ,0).*(1+p.omiga*exp(-x.^5/p.sigma^5)*A);
P_L=max(E_L-p.gI*exp(-x.^5/p.sigma^5)*I_L+p.g2*(L.T'*exp(-x.^5/(2*p.sigma)^5))'+p.g0*(R.T'*exp(-x.^5/(2*p.sigma)^5))' ,0).*(1+p.omiga*exp(-x.^5/p.sigma^5)*A);

%% updating using Eular
%% updating the monocular layer
% update the firing rate
P_normalize=sum([P_R.^2 P_L.^2],2);
L_new.T=L.T+1/p.tau*(-L.T+100*P_L.^2./((10+L.H).^2+repmat(P_normalize,1,2)))*dt;
L_new.H=L.H+1/p.tau_H*(-L.H+2*L.T)*dt;
R_new.T=R.T+1/p.tau*(-R.T+100*P_R.^2./((10+R.H).^2+repmat(P_normalize,1,2)))*dt;
R_new.H=R.H+1/p.tau_H*(-R.H+2*R.T)*dt;
%% update the binocular layer
% the binocular layer also adapted over time, so it has the parameter H, T
% is the firing rate, the structure is the same to the monocular, except
% the input is the firing rate of the two monocular layers
S_new.T=S.T+1/p.tau*( -S.T+100*(R.T+L.T).^2./((10+S.H).^2+repmat( sum((R.T+L.T).^2,2) ,1,2)) )*dt;
% first set the parameter to be the same as monocular layers
S_new.H=S.H+1/p.tau_H*(-S.H+S.T)*dt;
%% update all the opponency neurons
O_normalize=sum([max((R.T-L.T),0).^2 max((L.T-R.T),0).^2] ,2);
O_new.I_L=O.I_L+1/p.tau_I*(-O.I_L+30*max((R.T-L.T),0).^2./(p.alpha^2+repmat(O_normalize,1,2)))*dt;
O_new.I_R=O.I_R+1/p.tau_I*(-O.I_R+30*max((L.T-R.T),0).^2./(p.alpha^2+repmat(O_normalize,1,2)))*dt;
%% update the attention layer
% A is the attention layer activation, the first column is vertical, the
% second column is horizontal
% p.alpha_A is the normalization parameter for attention, p.tau_A is the
% parameter for attention
Dif=[S.T(:,1)-S.T(:,2) S.T(:,2)-S.T(:,1)];
A_new=A+1/p.tau_A*(-A+( sign(Dif).*(Dif.^2) )./(repmat(sum(Dif.^2,2),1,2)+p.alpha_A^2))*dt;



end