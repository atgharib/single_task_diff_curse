function [msdk,MSDdiffusion,Wdiffusion] = Diffusion_sigmaInsensitiveCorrentropy(sig,X,d,LMSiterations,N,M,A,C,Wo,muN)

%% Author: Atieh Gharib

Wdiffusion(:,1)=rand(M,1);
say = zeros(M,N,LMSiterations);

for n = 1 : LMSiterations
    for k = 1 : N
        sumdiff = zeros(M,1);
        for l=1:N
            sumdiff=sumdiff+A(l,k).*say(:,l,n); % Combination
        end
        phi(:,k,n) = sumdiff;
        %%  Adaptation Step for node k
        e_tmp = (d(k,n)-X(k,:,n)*phi(:,k,n));
        G = exp(-lambda*e_tmp.^2 /(2*sig^2));
        say(:,k,n+1)=phi(:,k,n)+ (1/(sig^2))*(muN)* e_tmp *G *exp(lambda*(1-G))*(1-lambda*G)*X(k,:,n)';

        if isnan(say(:,k,n+1))
            dispay(Err:NAN)
        end
    end
    
    %% Error analysis of Network
    Temp=zeros(M,1);
    for k=1:N
        Temp=Temp+(1/N).*say(:,k,n); % Averaging for computing MSD of total network.
    end
    Wdiffusion(:,n)=Temp; % Average over all the nodes in the network at time n
    Saybar(:,n)=Wo(1:M,n)-Wdiffusion(:,n);    % Error Vector        
    MSDdiffusion(n)=norm(Saybar(:,n))^2; % Scalar of error
    
    for k = 1 : N
        msdk(k,n)=norm(Wo(1:M,n)-say(:,k,n))^2;
    end

end
