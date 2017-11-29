function [alpha,lambda,a0,E_set, R_set, L_set,K_original,m_R,f] = svmpath_initialization_2016(data, label, ktype, params, rho) 



if (nargin <= 4)
    rho = .1;
end

K = zeros(length(label),length(label));

switch (ktype)
    case 'rbf'    
        lambda = params(1);
        for i = 1:length(label),
             K(i,:) = exp(-lambda * sum((data(i,:)'*ones(1,size(data',2)) - data').^2,1));
        end
    case 'linear'
        K = (data*data');
    otherwise
        error('Unrecognized kernel type');
end

K_original=K;

pex = sum(label==1);
mex = sum(label==-1);
idxaux = [];
E_set_all=zeros(pex+mex,1);

if (pex ~= mex),
    idxaux = length(label)+1:length(label) + abs(pex - mex);
    g = K*label;
    if (pex > mex),
        K = [K, g * ones(1,pex-mex) * -rho; 
             ones(pex-mex,1) * g' * -rho, rho^2 * sum(g.*label) * ones(pex-mex,pex-mex)];
        label = [label; -ones(pex-mex,1)];
    elseif (mex > pex)        
        K = [K, g * ones(1,mex-pex) * rho; 
             ones(mex-pex,1) * g' * rho, rho^2 * sum(g.*label) * ones(mex-pex,mex-pex)];
        label = [label; ones(mex-pex,1)];
    end
end




Q = K.*(label*label');

g = K*label;
[bp, ip] = max(g(label==1));
[bn, in] = min(g(label==-1));
ipa = find(label==1);
ina = find(label==-1);
lambda = (bp - bn)/2;
b0 = -(bp + bn)/(bp - bn);
a0 = b0 * lambda;
ip = ipa(ip);
in = ina(in);

num_nb = 2;
E_set = zeros(1,num_nb);

L_set = ones(length(label),1);
R_set = zeros(length(label),1);
E_set(1) = ip;
E_set(2) = in;

L_set(ip) = 0;
L_set(in) = 0;

alpha = ones(length(label),1);
f = 1/lambda * (g + a0);

m_R = 0;
m_R = UpdateCholesky(m_R, Q(ip,ip), label(ip));
m_R = UpdateCholesky(m_R, Q([ip in],[ip in]), label([ip in])');



if pex~=mex 
    
    iter=0;
    A_set=[];
    P_set=[];
   duplicate_ind=0;   
    

while (lambda > 1e-2 && ~isempty(L_set)),
    
    A_set(alpha(A_set)<1e-5)=[];
    ba = solveSub(m_R, Q(E_set,E_set), label(E_set)');
    b = ba(1:end-1);
    b0 = ba(end);
    
   
    lambdat=  zeros(num_nb,1);
    idxnbt = E_set;
    e = 0;
    lambdat(b < -e) = (1 +  lambda* b(b < -e) - alpha(idxnbt(b < -e))) ./(b(b < -e))  ;    
    lambdat(b > e) = ( lambda * b(b > e)- alpha(idxnbt(b > e)))./(b(b > e))  ; 
    
   % lambdat(lambdat>lambda)=0;
     
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
    h = K(:,E_set)*(b.*label(E_set)) + b0;
    lambdaf = zeros(length(L_set),1);
    T_h=label - h;
    T_active=find(abs(T_h)> 1e-6);
    
    lambdaf(T_active)=lambda * ((f(T_active) - label(T_active))./(T_h(T_active))) +lambda ; 
      
    lambdaf(idxnbt)=0;
    logl= (1 - label.*h) < -1e-6 & L_set;
    logr= (1 - label.*h) > 1e-6  & R_set;
    log=~( logl | logr);
    lambdaf=lambdaf.*log;
    
    lambdaf((lambdaf - lambda)/lambda > 1e-6) = 0;
    
    [lambdatmax, idxa] = max(lambdat);    
    [lambdafmax, idxf] = max(lambdaf);
    
    lambdap = max(lambdatmax, lambdafmax);
    

    if (lambdap > 1e-2), 
       
        % Perform incremental updates
        alpha(E_set) = alpha(E_set) - (lambda - lambdap) * b;
        a0 = a0 - (lambda - lambdap) * b0;
        f = lambda/lambdap * (f - h) + h;            
        
        if (lambdatmax > lambdafmax || ...
            (abs(lambdatmax - lambdafmax)/lambdafmax < 1e-9))% && idx_nb(idxa) < idxf)),

        
           %%%%%%%%%%%%%%
              removed_ind=E_set(idxa); 
              if( alpha(removed_ind) < 1e-6),
                    R_set(removed_ind) = 1;
              else
                    L_set(removed_ind) = 1;
              end
             
              m_R = DownDateCholesky(m_R, label(E_set)', idxa);
              E_set(idxa)=[];
              num_nb = num_nb - 1;
              
  
        
              if sum(P_set==removed_ind)
                 remove_iiii=[];
                 if ~isempty(A_set) 
                      [~,a_ind]=unique(sum( Q(A_set,E_set), 2));
                     
                             for ii=1:length(a_ind)
                                     idx_nb_temp=[E_set,A_set(a_ind(ii))];
                                     m_R_temp=UpdateCholesky(m_R, Q(idx_nb_temp, idx_nb_temp), label(idx_nb_temp)');
                                     if (min(abs(diag(m_R_temp))) > 1e-5)
                                                num_nb = num_nb + 1;
                                                E_set(num_nb) = A_set(a_ind(ii));
                                                if (L_set(A_set(a_ind(ii))) == 1), 
                                                    L_set(A_set(a_ind(ii))) = 0;
                                                else
                                                    R_set(A_set(a_ind(ii))) = 0;
                                                end
                                               m_R=m_R_temp;
                                            %   inActive_set(ii)=0;
                                               remove_iiii=[remove_iiii,a_ind(ii)];
                                                break;
                                              
                                     end
                             end
                            
                 end 
                  A_set(remove_iiii)=[];
%                   A_set(alpha(A_set)<1e-5)=[];
              end                
    
                
                
                
                
                
                
 
        else
%             num_nb = num_nb + 1;
            
            temp1=abs(lambdaf-lambdafmax);
            idx_addf=find (temp1 < 1e-12);   
            
            if length(idx_addf)==1
                    num_nb = num_nb + 1;
                    E_set(num_nb) = idxf;
                    if (L_set(idxf) == 1), 
                        L_set(idxf) = 0;
                    else
                        R_set(idxf) = 0;
                    end
                    [m_R] = UpdateCholesky(m_R, Q(E_set, E_set), label(E_set)');
                    A_set(A_set==idxf)=[];
            else
      %%%%%%%%%%%%%%%%%  for mutiple indexes
%                     inActive_set=idx_addf;
%                     inActive_set_all=inActive_set;
                    
                    A_set=unique([A_set;idx_addf]);
                    P_set=unique([P_set;idx_addf]);
                    
                    % check if they are duplicate
                    duplicate_ind=0;
                    ttt=sum(Q(idx_addf,idx_addf),2);
                    if norm(ttt-ttt(end:-1:1,:))< 1e-6
                       duplicate_ind=1;
                    end
                    
                    if duplicate_ind==1   % if duplicate
                             num_nb = num_nb + 1; 
                             E_set(num_nb) = idx_addf(1);
                             if (L_set(idx_addf(1)) == 1), 
                                  L_set(idx_addf(1)) = 0;
                             else
                                  R_set(idx_addf(1)) = 0;
                             end
                             [m_R] = UpdateCholesky(m_R, Q(E_set, E_set), label(E_set)');
                             A_set(A_set==idx_addf(1))=[];
                    else       % if not duplicate
                             for ii=1:length(idx_addf)
                                     idx_nb_temp=[E_set,idx_addf(ii)];
                                     m_R_temp=UpdateCholesky(m_R, Q(idx_nb_temp, idx_nb_temp), label(idx_nb_temp)');
                                     if (min(abs(diag(m_R_temp))) > 1e-5)
                                                num_nb = num_nb + 1;
                                                E_set(num_nb) = idx_addf(ii);
                                                if (L_set(idx_addf(ii)) == 1), 
                                                    L_set(idx_addf(ii)) = 0;
                                                else
                                                    R_set(idx_addf(ii)) = 0;
                                                end
                                               m_R=m_R_temp;
                                               A_set(A_set==idx_addf(ii))=[];
                                      end
                                      
                             end
                          
                    end
                    
                    
                      
              end
        %%%%%%%%%%%%%
        end

        
    end 
  
        iter=iter+1;
        lambda = lambdap;
        if (max(alpha(idxaux)) < 1e-6 &   ~isempty(E_set(E_set<= pex+ mex ) )   )
           add_no=abs(pex- mex);
           alpha(end-add_no+1:end)=[];
           R_set(end-add_no+1:end)=[];
           L_set(end-add_no+1:end)=[];
           f(end-add_no+1:end)=[];
           E_set(E_set>length(alpha))=[];

            break;
        end

    
end

end







function out = solveSub(m_R, Q, y)
   if (length(y) > 1)
        
        Z = [-y(1)*y(2:end); eye(length(y)-1)];
        
        rhs = Z'*ones(length(y),1);
        
        hz = m_R'\rhs;  %2.
        hz = m_R\hz;
                
        h = Z*hz;       %3.                 
        g = y(1)*(1 - Q(1,:)*h);  %g is the beta_0 in (17)
        
        out = [h; g];
    else
        % Just solve the following
        % 1. h = y(1)*r
        % 2. g = y(1)*(q - Q*h)
        %         

        out = [0; y(1)];
    end














