
clc; close all; clear all;
n=30; p.tau=20; p.tau_I=10; p.tau_H=900; p.tau_A=50; % p.tau_H is also used in the binocular layer
p.gI=0.27;p.g0=0; p.g1=0.02; p.g2=0.02; p.sigma=2; p.alpha=10;p.alpha_A=10;p.omiga=0.0;dt=10;% the p.omiga should be set to 
E_L=[repmat(70,n,1) zeros(n,1)]; E_R=[zeros(n,1) repmat(70,n,1)];
% initialize the model 
L.T=rand(n,2)/100;L.H=zeros(n,2);
R.T=rand(n,2)/100;R.H=zeros(n,2);
O.I_L=zeros(n,2);O.I_R=zeros(n,2);
S.T=zeros(n,2);S.H=zeros(n,2);
A=zeros(n,2);
itr1=200; % to get the steady state
itr2=1000; % the simulation of 
for i=1:12
    plot_data{i}=zeros(itr1+itr2,n);
end
for i=1:itr1
    [L,R,O,S,A]=bar_update(R,L,O,S,A,p,E_L,E_R,dt);
    plot_data{1}(i,:)=L.T(:,1)'; plot_data{2}(i,:)=L.T(:,2)'; plot_data{3}(i,:)=R.T(:,1)'; plot_data{4}(i,:)=R.T(:,2)';
    plot_data{5}(i,:)=O.I_L(:,1)'; plot_data{6}(i,:)=O.I_L(:,2)'; plot_data{7}(i,:)=O.I_R(:,1)';plot_data{8}(i,:)=O.I_R(:,2)';
    plot_data{9}(i,:)=S.T(:,1)'; plot_data{10}(i,:)=S.T(:,2)';
    plot_data{11}(i,:)=A(:,1)'; plot_data{12}(i,:)=A(:,2)';
end
R.T(1:3,2)=R.T(1:3,2)+100;
[L,R,O,S,A]=bar_update(R,L,O,S,A,p,E_L,E_R,dt);
for i=(itr1+1):(itr2+itr1)
    [L,R,O,S,A]=bar_update(R,L,O,S,A,p,E_L,E_R,dt);
    plot_data{1}(i,:)=L.T(:,1)'; plot_data{2}(i,:)=L.T(:,2)'; plot_data{3}(i,:)=R.T(:,1)'; plot_data{4}(i,:)=R.T(:,2)';
    plot_data{5}(i,:)=O.I_L(:,1)'; plot_data{6}(i,:)=O.I_L(:,2)'; plot_data{7}(i,:)=O.I_R(:,1)';plot_data{8}(i,:)=O.I_R(:,2)';
    plot_data{9}(i,:)=S.T(:,1)'; plot_data{10}(i,:)=S.T(:,2)';
    plot_data{11}(i,:)=A(:,1)'; plot_data{12}(i,:)=A(:,2)';
end

c_mon=[0,100];
subplot(3,4,1);
imagesc(plot_data{1}',c_mon);
title('Left eye vertical response');
hold on;
subplot(3,4,2);
imagesc(plot_data{2}',c_mon);
title('Left eye horizontal response');
hold on;
subplot(3,4,5);
imagesc(plot_data{3}',c_mon);
title('right eye vertical response');
hold on;
subplot(3,4,6);
imagesc(plot_data{4}',c_mon);
title('Right eye horizontal response');
hold on;
c_opo=[0,30];
subplot(3,4,3);
imagesc(plot_data{5}',c_opo);
title('oppency neuron supress left with vertical input')
hold on;
subplot(3,4,4);
imagesc(plot_data{6}',c_opo);
title('oppency neuron supress left with horizontal input')
hold on;
subplot(3,4,7);
imagesc(plot_data{7}',c_opo);
title('oppency neuron supress right with vertical input')
hold on;
subplot(3,4,8);
imagesc(plot_data{8}',c_opo);
title('oppency neuron supress right with horizontal input')
c_bin=[0,100];
subplot(3,4,9);
imagesc(plot_data{9}',c_bin);
title('binocular neuron vertical')
hold on;
subplot(3,4,10);
imagesc(plot_data{10}',c_bin);
title('binocular neuron horizontal')
hold on;

c_att=[-1,1];
subplot(3,4,11);
imagesc(plot_data{11}',c_att);
title('attention neuron vertical')
hold on;
subplot(3,4,12);
imagesc(plot_data{12}',c_att);
title('attention neuron horizontal')
