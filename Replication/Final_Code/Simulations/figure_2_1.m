%%% AUTHOR:  STUART LANE
%%% DATE:    07/06/2023 
%%% PAPER:   OVERIDENTIFICATION TESTING WITH WEAK INSTRUMENTS AND HETERO-
%%%          SKEDASTICITY 
%%% CONTENT: FIGURE 2.1

clear;

spmd
    rng(234435,'combRecursive');
end
                
tStart = tic;

n = 120;
Nrep = 20000;

b0 = 0;
b1 = 0;

Kvals = 3;

rhoa = [0.2 0.5 0.95];
alpha = 0.5;
deltaa = [-30:1:-3,-2.5:0.5:2.5,3:30];

mu2vals = 48;

tic
sb = zeros(Nrep,length(deltaa),length(rhoa),2);

for Kj = 1:length(Kvals)

    K = Kvals(Kj);

    for mu2j = 1:length(mu2vals)

        mu2 = mu2vals(mu2j);
        pp = sqrt(mu2/(n*(K-1)));
    
     for deltaj = 1:length(deltaa)

         delta = deltaa(deltaj);
         dd = delta/sqrt(n);
            
         for rhoj = 1:length(rhoa)

             rho = rhoa(rhoj);  
                
             parfor rep = 1:Nrep

                stream = RandStream.getGlobalStream();
                stream.Substream = rep;
                
                warning('off');
                
                z = randn(n,1);
                 
                muuv = [0;0];
                sigmauv = [1 rho; rho 1];
                eps = mvnrnd(muuv,sigmauv,n);
                
                u = dd*z + (abs(z).^(alpha)).*eps(:,1);
                v = (abs(z).^(alpha)).*eps(:,2);  
                
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

                sb(rep,deltaj,rhoj,:) = [score2sls kpsample];

            end

         end

     end
       
     deltan = deltaa/sqrt(n);
    
     deltak = 1:length(deltaa);
    
     graph = figure
     p11 = plot(deltan,mean(chi2cdf(sb(:,deltak,1,1),K-2,'upper')<0.05),'b')
     hold on
     p13 = plot(deltan,mean(chi2cdf(sb(:,deltak,2,1),K-2,'upper')<0.05),'b--')
     hold on
     p15 = plot(deltan,mean(chi2cdf(sb(:,deltak,3,1),K-2,'upper')<0.05),'b-.')
     hold on
     p12 = plot(deltan,mean(chi2cdf(sb(:,deltak,1,2),K-2,'upper')<0.05),'r')
     hold on
     p14 = plot(deltan,mean(chi2cdf(sb(:,deltak,2,2),K-2,'upper')<0.05),'r--')
     hold on
     p16 = plot(deltan,mean(chi2cdf(sb(:,deltak,3,2),K-2,'upper')<0.05),'r-.')
     hline = refline(0,0.05) 
     hline.Color = 'k'
     xlim([-0.6 0.6])
     ylim([0 1])
     xlabel('$\omega$','Interpreter','Latex') 
     ylabel('Power')
     legend('$J (\rho = 0.2)$','$J (\rho = 0.5)$','$J (\rho = 0.95)$','$KP (\rho = 0.2)$','$KP (\rho = 0.5)$','$KP (\rho = 0.95)$','Location','north','Interpreter','Latex')
    
     end

end

toc 