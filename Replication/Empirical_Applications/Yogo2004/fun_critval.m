function cv = fun_critval(Keff,tau)
    rng(234435,'combRecursive');
    vals = ncx2rnd(Keff,tau*Keff,100,1)/Keff; % SET TO 10000000 FOR ACTUAL RESULTES
    cv = prctile(vals,95);
end
    
