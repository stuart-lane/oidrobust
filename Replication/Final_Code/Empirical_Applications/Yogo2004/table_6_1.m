%%% AUTHOR:  STUART LANE
%%% DATE:    07/06/2023 
%%% PAPER:   OVERIDENTIFICATION TESTING WITH WEAK INSTRUMENTS AND HETERO-
%%%          SKEDASTICITY 
%%% CONTENT: TABLE 6.1
%%% DATASET: YOGO (2004)

clear; 

%%% SAVE RESULTS AT END? SET == 1 IF YES, == 0 IF NO:
SaveResults = 1;

%%

tic

rng(234435,'combRecursive');
     
print = zeros(1,8);
Table6_1a = zeros(1,6);
Table6_1b = zeros(1,6);

num_sim = 10;

for country_number = 1:11

    % SELECT COUNTRY
    if country_number == 1
       dataset = importdata('Data\AULQ.txt');
    elseif country_number ==2
       dataset = importdata('Data\CANQ.txt');
    elseif country_number ==3
       dataset = importdata('Data\FRQ.txt');
    elseif country_number ==4
       dataset = importdata('Data\GERQ.txt');
    elseif country_number ==5
       dataset = importdata('Data\ITAQ.txt');
    elseif country_number ==6
       dataset = importdata('Data\JAPQ.txt');
    elseif country_number ==7
       dataset = importdata('Data\NTHQ.txt');
    elseif country_number ==8
       dataset = importdata('Data\SWDQ.txt');
    elseif country_number ==9
       dataset = importdata('Data\SWTQ.txt');
    elseif country_number ==10
       dataset = importdata('Data\UKQ.txt');
    else
        dataset = importdata('Data\USAQ.txt');
    end
    
    % SELECT NUMBER OF LAGS DEPENDING ON COUNTRY
    if country_number <=10
        L = 4;
    else
        L = 6;
    end
    
    dataset = dataset.data;

    dat = dataset(:,1);
    dc = dataset(:,6);
    r = dataset(:,7);
    ir = dataset(:,8);
    Z = dataset(:,9:12);

    obs = length(dataset);
    onevec = ones(obs,1);

    M = eye(obs) - onevec*onevec'/obs;
    dc = M*dc;
    r = M*r;
    ir = M*ir;
    Z = M*Z;
    Z = fun_orthonormalise(Z,obs);
    izz = inv(Z'*Z);
    
    %%%%%%%%%%%%%%
    % TABLE 6.1a %
    %%%%%%%%%%%%%%
    
    pi_tsls_ir = Z\dc;
    tsls_ir = (Z*pi_tsls_ir)\ir; 
    u_tsls_ir = ir - dc*tsls_ir;      
    v_tsls_ir = dc - Z*pi_tsls_ir;
    
    rho_tsls_ir = corr(u_tsls_ir,v_tsls_ir);
    cov_tsls_ir = cov(u_tsls_ir,v_tsls_ir);
    cov_tsls_ir = cov_tsls_ir(1,2);
    
    [liml_ir,pi_liml_ir,u_liml_ir,v_liml_ir] = fun_liml(ir,dc,Z,izz,obs,1);
    rho_liml_ir = corr(u_liml_ir,v_liml_ir);
    
    % COMPUTE EFFECTIVE F STATISTIC
    Wir = NeweyWest(v_tsls_ir,Z,L,0);
    trWir = trace(Wir);
    effFir = dc'*Z*Z'*dc/trace(Wir);
    effDOFir = (trWir^2)*(1+2*10)/(trace(Wir'*Wir) + 2*10*trWir*max(eig(Wir)));
    cvir = fun_critval(effDOFir,10);

    J_tsls_ir = fun_robust_score_nw(Z,Z(:,2:4),pi_tsls_ir,u_tsls_ir,L);
    KP_liml_ir = fun_robust_score_nw(Z,Z(:,2:4),pi_liml_ir,u_liml_ir,L);

    %%%%%%%%%%%%%%
    % TABLE 6.1b %
    %%%%%%%%%%%%%%
    
    pi_tsls_dc = Z\ir;
    tsls_dc = (Z*pi_tsls_dc)\dc;
    u_tsls_dc = dc - ir*tsls_dc;
    v_tsls_dc = ir - Z*pi_tsls_dc;
    
    rho_tsls_dc = corr(u_tsls_dc,v_tsls_dc);
    cov_tsls_dc = cov(u_tsls_dc,v_tsls_dc);
    cov_tsls_dc = cov_tsls_dc(1,2);

    [liml_dc,pi_liml_dc,u_liml_dc,v_liml_dc] = fun_liml(dc,ir,Z,izz,obs,1);
    rho_liml_dc = corr(u_liml_dc,v_liml_dc);

    J_tsls_dc = fun_robust_score_nw(Z,Z(:,2:4),pi_tsls_dc,u_tsls_dc,L);
    KP_liml_dc = fun_robust_score_nw(Z,Z(:,2:4),pi_liml_dc,u_liml_dc,L);
    
    % COMPUTE EFFECTIVE F STATISTIC
    Wdc = NeweyWest(v_tsls_dc,Z,L,0);
    trWdc = trace(Wdc);
    effFdc = ir'*Z*Z'*ir/trace(Wdc);
    effDOFdc = (trWdc^2)*(1+2*10)/(trace(Wdc'*Wdc) + 2*10*trWdc*max(eig(Wdc)));
    cvdc = fun_critval(effDOFdc,10);
    
    % PUT VALUES INTO TABLES
      
    Table6_1a(end+1,:) = [effFdc cvdc tsls_dc liml_dc  J_tsls_dc KP_liml_dc]
    Table6_1b(end+1,:) = [effFir cvir tsls_ir liml_ir J_tsls_ir KP_liml_ir]
 
    while country_number < 11
        country_number = country_number + 1; 
    end
    
end

Table6_1a_print = round(Table6_1a(2:end,:),2)
Table6_1b_print = round(Table6_1b(2:end,:),2)

toc