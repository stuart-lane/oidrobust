%%% AUTHOR:  STUART LANE
%%% DATE:    07/06/2023 
%%% PAPER:   OVERIDENTIFICATION TESTING WITH WEAK INSTRUMENTS AND HETERO-
%%%          SKEDASTICITY 
%%% CONTENT: TABLE E.1

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
alphfac = [0.8,1,1.3];
muuv = [0;0];

mu2a = [1 4 8 16 32];

tic

opts = optimset('Display','off');

printresults = zeros(1,12);

for Kj = 1:length(Kvals)
    K = Kvals(Kj);
    
    for rhoj = 1:length(rhovals)
        rho = rhovals(rhoj);
        sigmauv = [1 rho; rho 1];

        for alphaj = 1:length(alphaa)
            alpha = alphaa(alphaj);
            fac = alphfac(alphaj);

             for mu2j = 1:length(mu2a)
                 mu2 = mu2a(mu2j);
                 pp = sqrt(fac*mu2/(n*(K-1)));

                 sb = zeros(Nrep,4);

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
                         u = (abs(z).^alpha).*uv(:,1);  
                         v = (abs(z).^alpha).*uv(:,2);
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

                     izz = inv(z'*z);

                     b2sls = (z*(z\x))\y;
                     u2sls = y - x*b2sls;
                     pi2sls = z\x;
                     v2sls = x - z*(z\x);

                     x2sls =  z*pi2sls;

                     [bliml,piliml,uliml,vliml] = fun_liml(y,x,z,izz,n,kx);

                     kpsample = fun_robust_score(z,z2,piliml,uliml);

                     score2sls = fun_robust_score(z,z2,pi2sls,u2sls);

                     sb(rep,:) = [b2sls(2) bliml(2) score2sls kpsample];
                
                 end

                 meds = median(sb(:,1:2));
                 nines = prctile(sb(:,1:2),90) - prctile(sb(:,1:2),10);
                 pct10 = mean(chi2cdf(sb(:,3:4),K-2,'upper')<0.1);
                 pct1 = mean(chi2cdf(sb(:,3:4),K-2,'upper')<0.01);
                 printresults(end+1,:) = [K rho alpha mu2 meds nines pct10 pct1]

             end
             
        end
        
    end
    
end

toc 

% printtable = printresults(2:end,:)
% T = round(printtable(:,4:end),3)
% 
% % unblock line below to save table of results
% csvwrite('tailvalues.csv',T(:,4:end),0,1)
% 
% T1 = array2table(T)
% table2latex(T1)
% 
% tableT = table2latex(T,'table.tex')