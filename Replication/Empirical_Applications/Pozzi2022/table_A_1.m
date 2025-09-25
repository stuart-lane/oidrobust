%%% AUTHOR:  STUART LANE
%%% DATE:    07/06/2023 
%%% PAPER:   OVERIDENTIFICATION TESTING WITH WEAK INSTRUMENTS AND HETERO-
%%%          SKEDASTICITY 
%%% CONTENT: TABLE A.1
%%% DATASET: POZZI (2022)

clear; 

%% DATA PREPARATION

tic

rng(234435,'combRecursive');

dataset = importdata('dataset.txt');

dataset = dataset.data;

Country = {'AUS','DEN','FIN','FRA','GER','ITA','JAP','NTH','NOR','POR','SPA','SWD','SWT','UK','USA'}';

year = dataset(:,1);
numyears = max(year) - min(year);
totequityreturn = dataset(:,2);
dc = dataset(:,end);
rhr = dataset(:,5);

fulldataset = zeros(66,6,15);
for j = 1:15
    index = (((j-1)*66)+1):(j*66);
    fulldataset(:,:,j) = dataset(index,:);
end

%% RESULTS

lag_var = [1,2];
for ll = 1:2

    lagt = lag_var(ll);
    
    comp = zeros(1,6);
    
    for cn = 1:15

        country_data = fulldataset(:,:,cn);

        a1 = isnan(country_data);
        a2 = any(a1,2);
        a3 = ~a2;

        country_data = country_data(a3,:);

        dc = country_data(2:end,end); % change in consumption
        rhr = country_data(2:end,5); % real housing return
        averagerate = nanmean(fulldataset(a3,5,[1:cn-1 cn+1:end]),3);

        lagrate = rhr(1:end-lagt);
        x = rhr((lagt+1):end);
        y = dc((lagt+1):end);
        z = [averagerate((lagt+2):end) lagrate];
        n = length(x);

        M = eye(n) - ones(n,n)/n;
        x = M*x;
        y = M*y;
        z = M*z;

        % this loop orthonomalises the instruments
        for orthonormalise_instruments = 1
            
            z1 = z(:,1);
            z2 = z(:,2);
            Z1 = z1/norm(z1);
            s2 = (z2'*Z1)*Z1;
            e2 = z2 - s2;
            Z2 = e2/norm(e2);
            Zall = [Z1 Z2];
            z = Zall*sqrt(n);
            
        end
        
        izz = inv(z'*z);

        b2sls = (z*(z\x))\y;
        u2sls = y - x*b2sls;
        pi2sls = z\x;
        v2sls = x - z*(z\x);
        x2sls =  z*pi2sls;
        [bliml,piliml,uliml,vliml] = fun_liml(y,x,z,izz,n,1);

        z2 = z(:,2:end);
        xhat = z*piliml;
        mxz = z2 - xhat*(xhat\z2);

        W = NeweyWest(v2sls,z,3,0);
        trW = trace(W);
        effF = x'*z*z'*x/(trace(W));
        effDOF = (trW^2)*(1+2*10)/(trace(W'*W) + 2*10*trW*max(eig(W)));
        cv = fun_critval(effDOF,10);

        kpsample = fun_robust_score_nw(z,z2,piliml,uliml,3);

        score2sls = fun_robust_score_nw(z,z2,pi2sls,u2sls,3); 

        comp(end+1,:) = [effF cv b2sls bliml score2sls kpsample];

    end
    
    if ll == 1 

        printtable = comp(2:end,:);
        T1 = round(printtable,2)
        
    else
        
        printtable = comp(2:end,:);
        T2 = round(printtable,2)

    end
    
end

toc