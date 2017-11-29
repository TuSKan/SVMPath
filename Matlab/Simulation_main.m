clear;
InputFile = 'newdata\wdbc.data';
data= load(InputFile);
data(:,1)=[];



mode=2;       %1:linear  2:RBF   3:Polynomial
switch (mode)
    case 1    
         opt_kernel='linear';
    case 2
        opt_kernel='rbf';
end



rand('seed',0)
%%%% add duplicate data points
percent=0.1;
ind_duplicate=[];
total_points=length(data(:,1));
no_duplicate=round(total_points*percent);

for ii=1:no_duplicate
    ind_duplicate=[ind_duplicate,randperm(total_points,1)];
end
data=[data;data(ind_duplicate,:)];

label=data(:,1);
data=data(:,2:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  data normalize
label=2*label-1;
datamean=mean(data);
datastd=std(data);
for i=1:size(data,2)
    if datastd(i)<10^(-8)
        datastd(i)=1;
    end
end
data= (data - ones(size(data,1),1)*datamean)./ (ones(size(data,1),1)*datastd) ;




InputFile=data_reformulate(data,label);
%%%%%%%%%%%%%%%%%%%%%%---------------------------------------------%%%%%%%%%%%%%%


[alpha,lambda,a0,E_set, R_set, L_set,K,m_R,f] = svmpath_initialization_2016(data, label, opt_kernel,0.1);
[cost_our,lambd_all,max_nb,max_Inactive,mun_non] = svmpath_subsectmethod_5(K, label,alpha,lambda,a0,E_set, R_set, L_set,m_R,f); 

 plot(lambd_all,cost_our,'o')


% tic
% update_QR=1;      %1: Update QR; 0: inv
% Optimality_check=0;
% [cost_old,lambda_old,Length_E,run_time]=SVMpath_our(data',label,mode,update_QR,Optimality_check);
% toc




