function [cost_our,lambdavec,max_nb,max_Inactive,mun_non] = svmpath_subsectmethod_5(K, label,alpha,lambda,a0, E_set, R_set, L_set,m_R,f) 

% [cost_our,lambdavec,alpha,lambda,E_set, R_set, L_set] = svmpath_subsectmethod_5(K, label,alpha,lambda,a0, E_set, R_set, L_set,m_R,f) 

Q = K.*(label*label');
num_nb=length(E_set);

iter=0;

A_set=[];
P_set=[];
nrepeatiter=0;
no_algorithm1=0;
count_ind=1;
lambdavec=[];
max_nb=0;
max_Inactive=0; num_Inactive=0;
cost_our=[];
duplicate_ind=0;

mun_non=0;


while (lambda > 1e-2 && ~isempty(L_set)  ),
%      [iter,lambda]
%     
%     if iter>=4641
%         111;
%     end
%     
%     if (sum(A_set==3224))
%         111;
%     end
% %     

    if (~checkKKT(alpha, label, f, L_set, R_set)),
        disp('algorithm failure');
 %         Find the offending point and move to idx_nb, etc.
    end
    
    A_set(alpha(A_set)<1e-5)=[];
    num_Inactive=length(A_set);  
    
  max_Inactive = max(max_Inactive, num_Inactive); 
    
  if (count_ind==1)
      iter=iter+1;
      max_nb = max(max_nb, num_nb);
      lambdavec(iter) = lambda;
      cost_our(iter)= -(1/(2*lambda)*alpha'*Q*alpha -sum(alpha));
      
      if ~isempty(A_set)
          mun_non=mun_non+1;
      end
      
      if iter>1 && (abs( lambdavec(iter-1) - lambda)/lambda < 1e-9)
         nrepeatiter=nrepeatiter+1; 
      end
      
  end
   
   count_ind=1;
   

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
                      [~,a_ind]=unique(sum( Q(A_set,E_set), 2)); %%%%%%%%%
                      
                     %%%%%%%%%%+++++++
                     if length(a_ind)>1
                       no_algorithm1=no_algorithm1+1;
                     end
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
                                                count_ind=0;
%                                                break;
                                              
                                     end
                             end
                            
                 end 
                  A_set(remove_iiii)=[];
%                   A_set(alpha(A_set)<1e-5)=[];
                  num_Inactive=length(A_set);
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

            %       idx_addf(alpha(idx_addf)<1e-5)=[];
                    A_set=unique([A_set;idx_addf]);
                    P_set=unique([P_set;idx_addf]);
                    
                    % check if they are duplicate
                    duplicate_ind=0;
                    ttt=sum(Q(idx_addf,idx_addf),2);  %%%%%%%%%
                    if norm(ttt-ttt(end:-1:1,:))< 1e-6
                       duplicate_ind=1;
                    end
                    
                      %%%%%%%%%%+++++++
                  
                    
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
                             num_Inactive=length(A_set);
                    else       % if not duplicate
                      
                        if length(idx_addf)>1
                           no_algorithm1=no_algorithm1+1;
                        end
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
                                               num_Inactive=length(A_set);
                                               count_ind=0;
                                      end
                                      
                             end
                            
                    end
                    
                    
                      
              end
        %%%%%%%%%%%%%
        end

        
    end 
  
    
    lambda = lambdap;

    
end


no_algorithm1
nrepeatiter


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

