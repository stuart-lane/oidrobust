%%% AUTHOR:  STUART LANE
%%% DATE:    07/06/2023 
%%% PAPER:   OVERIDENTIFICATION TESTING WITH WEAK INSTRUMENTS AND HETERO-
%%%          SKEDASTICITY 
%%% CONTENT: FIGURE 4.1

clear;

spmd
    rng(234435,'combRecursive');
end
                
tStart = tic;

n = 120;
Nrep = 20000;

b0 = 0;
b1 = 0;

Kvals = [3 5];
rhovals = [0.2 0.5 0.95];

alphaa = [0,1,2];
alphfac = [0.8,1,3];
muuv = [0;0];

mu2a = [0.1,0.2,0.4,0.6,0.8,1:32];
mu2k = 1:length(mu2a);

tic

for Kj = 1:length(Kvals)
    
    K = Kvals(Kj);
    
    for rhoj = 1:length(rhovals)
        
        rho = rhovals(rhoj);
        sigmauv = [1 rho; rho 1];
        sb = zeros(Nrep,length(mu2a),length(alphaa),4);
        
        for mu2j = 1:length(mu2a)
            
            mu2 = mu2a(mu2j);
        
            for alphaj = 1:length(alphaa)
                
                alpha = alphaa(alphaj);
                fac = alphfac(alphaj);
                pp = sqrt(fac*mu2/(n*(K-1)));
             
                parfor rep = 1:Nrep
                    
                       stream = RandStream.getGlobalStream();
                       stream.Substream = rep;

                       warning('off');

                       z = randn(n,1);
                       uv = mvnrnd(muuv,sigmauv,n);
                       
                       if alpha == 0
                           
                            u = (abs(z).^(0.5)).*uv(:,1);
                            v = (abs(z).^(0.5)).*uv(:,2);
                            
                       else                
                           
                            u = (z.^alpha).*uv(:,1);  
                            v = (z.^alpha).*uv(:,2);
                            
                       end
                       
                       z = [z randn(n,K-2)];
                
                       PP = pp*ones(K-1,1);
                
                       x = z*PP+v;
                       y = b0+x*b1+u;

                       z = [ones(n,1) z];
                       x = [ones(n,1) x];

                       kx = size(x,2);
                       kz = size(z,2);
                       z2 = z(:,kx+1:end);

                       izz = inv(z'*z)

                       b2sls = (z*(z\x))\y;
                       u2sls = y - x*b2sls;
                       pi2sls = z\x;
                       v2sls = x - z*(z\x);

                       x2sls =  z*pi2sls;
                       
                       [bliml,piliml,uliml,vliml] = fun_liml(y,x,z,izz,n,kx);

                       xhat = z*piliml;
                       mxz = z2 - xhat*(xhat\z2);
                       test = mxz'*uliml
                       test2 = z2'*uliml;

                       kpsample = fun_robust_score(z,z2,piliml,uliml);

                       score2sls = fun_robust_score(z,z2,pi2sls,u2sls);           

                       sb(rep,mu2j,alphaj,:) = [b2sls(2) bliml(2) score2sls kpsample];
                end
                
            end
            
        end 
   
        graph = figure
        p1 = plot(mu2a,mean(chi2cdf(sb(:,mu2k,1,3),K-2,'upper')<0.05),'b')
        hold on
        p2 = plot(mu2a,mean(chi2cdf(sb(:,mu2k,2,3),K-2,'upper')<0.05),'b--')
        hold on
        p3 = plot(mu2a,mean(chi2cdf(sb(:,mu2k,3,3),K-2,'upper')<0.05),'b:')
        hold on
        p4 = plot(mu2a,mean(chi2cdf(sb(:,mu2k,1,4),K-2,'upper')<0.05),'r')
        hold on
        p5 = plot(mu2a,mean(chi2cdf(sb(:,mu2k,2,4),K-2,'upper')<0.05),'r--')
        hold on
        p6 = plot(mu2a,mean(chi2cdf(sb(:,mu2k,3,4),K-2,'upper')<0.05),'r:')
        hold on
        hline = refline(0,0.05) 
        hline.Color = 'k'
        ylim([0 0.5])
        xlim([0 32])
        xlabel('$\mu^2$','Interpreter','Latex') 
        ylabel('Rejection Frequency')
        legend('$J$ $(\alpha = 0.5)$','$J$ $(\alpha = 1)$','$J$ $(\alpha = 2)$','$KP$ $(\alpha = 0.5)$','$KP$ $(\alpha = 1)$','$KP$ $(\alpha = 2)$','Location','northeast','Interpreter','Latex')

    end

end