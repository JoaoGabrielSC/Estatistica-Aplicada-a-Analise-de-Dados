function ctr = contrib_pc(x, P, a, N, alfa, L)
% Contribution using PCA for each sample x
% Inputs:
% x = sample dimension mx1
% P = matrix of eigenvectors associated to the eigenvalues L
% L = vector highest retained eigenvalues
% a = number of principal components (scalar)
% N = number of samples used to compute the covariance matrix (scalar)
% alfa =confidence level (Example: 0.95)
% 
T=x*P; % one sample of x
[m,c]=size(P); % m = number of variables
ctr=zeros(1,m);
idx=[];
t2=(a*(N-1)*(N+1)/(N*(N-a)))*finv(alfa,a,N-a); % T2 threshold
for j=1:c % c scores from c retained principal components 
    if (((T(j)/sqrt(L(j)))^2)>(1/a)*t2)
        idx=[idx j];  % scores that violate threshold
    end
end;
cont=[];
c=length(idx); % c selected scores  (threshold violated)
if c>0 % If at least one score was violated
    for i=1:c      % computation for each score
        for j=1:m  % Contribution of m variables to score ti
            tn=idx(i);
            ti=T(tn);
            pij=P(j,tn);
            aux=(ti/L(tn))*pij*x(j);
            if aux>0
                cont(i,j)=aux;
            else
                cont(i,j)=0;
            end
        end
    end
    if c>1 cont= sum(cont); end; % Add contribution of m variables 
    ctr= cont/sum(cont); 

end 