%%% AUTHOR:  STUART LANE
%%% DATE:    07/06/2023 
%%% PAPER:   OVERIDENTIFICATION TESTING WITH WEAK INSTRUMENTS AND HETERO-
%%%          SKEDASTICITY 
%%% CONTENT: FIGURE 4.2

clear;

spmd
    rng(234435,'combRecursive');
end
                
tStart = tic;

n = 120;
Nrep = 20000;

b0 = 0;
b1 = 0;

Kvals = [3 5]; % K includes constant

alpha = 1;

muuv = [0;0];

rhoa = [-0.990 -0.99:0.05:0.99,0.999];
rhok = 1:length(rhoa);
mu2a = [1;8;16];

tic
sb = zeros(Nrep,length(rhoa),length(mu2a),4);
        
for Kj = 1:length(Kvals)
    K = Kvals(Kj);
    
    for rhoj = 1:length(rhoa)
        rho = rhoa(rhoj);
        sigmauv = [1 rho; rho 1];
        
        for mu2j = 1:length(mu2a)
            mu2 = mu2a(mu2j);
            pp = sqrt(mu2/(n*(K-1)));
             
            parfor rep = 1:Nrep
                stream = RandStream.getGlobalStream();
                stream.Substream = rep;
                
                warning('off');
                
                z = randn(n,1);

                uv = mvnrnd(muuv,sigmauv,n);              
                u = z.*uv(:,1);
                v = z.*uv(:,2);     
                
                z = [z randn(n,K-2)];
                
                PP = pp*ones(K-1,1);

                x = z*PP+v;
                y = b0+x*b1+u;

                z = [ones(n,1) z];
                x = [ones(n,1) x];
               
                kx = size(x,2);
                kz = size(z,2);
                z2 = z(:,kx+1:end);

                izz = inv(z'*z);


                b2sls = (z*(z\x))\y;
                u2sls = y - x*b2sls;
                pi2sls = z\x;
                v2sls = x - z*(z\x);

                x2sls =  z*pi2sls;

                [bliml,piliml,uliml,vliml] = fun_liml(y,x,z,izz,n,kx);

                 kpsample = fun_robust_score(z,z2,piliml,uliml);
                 
                 score2sls = fun_robust_score(z,z2,pi2sls,u2sls);      

                 sb(rep,rhoj,mu2j,:) = [b2sls(2) bliml(2) score2sls kpsample];
            end

        end
        
    end
        
    graph = figure
    p31 = plot(rhoa,mean(chi2cdf(sb(:,rhok,1,3),K-2,'upper')<0.05),'b:')
    hold on
    p32 = plot(rhoa,mean(chi2cdf(sb(:,rhok,2,3),K-2,'upper')<0.05),'b--')
    hold on
    p33 = plot(rhoa,mean(chi2cdf(sb(:,rhok,3,3),K-2,'upper')<0.05),'b')
    hold on
    p34 = plot(rhoa,mean(chi2cdf(sb(:,rhok,1,4),K-2,'upper')<0.05),'r:')
    hold on
    p35 = plot(rhoa,mean(chi2cdf(sb(:,rhok,2,4),K-2,'upper')<0.05),'r--')
    hold on
    p36 = plot(rhoa,mean(chi2cdf(sb(:,rhok,3,4),K-2,'upper')<0.05),'r')
    hold on
    hline = refline(0,0.05) 
    hline.Color = 'k'
    ylim([0 0.5])
    xlim([-1 1])
    xlabel('$\rho$','Interpreter','Latex') 
    ylabel('Rejection Frequency')
    legend('$J$ ($\mu^2=1$)','$J$ ($\mu^2=8$)','$J$ ($\mu^2=16$)','$KP$ ($\mu^2=1$)','$KP$ ($\mu^2=8$)','$KP$ ($\mu^2=16$)','Location','northwest','Interpreter','Latex')

end

toc