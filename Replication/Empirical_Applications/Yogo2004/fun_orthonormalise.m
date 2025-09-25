function orthoZ = fun_orthonormalise(Z,obs)
    z1 = Z(:,1);
    z2 = Z(:,2);
    z3 = Z(:,3);
    z4 = Z(:,4);
    Z1 = z1/norm(z1);
    s2 = (z2'*Z1)*Z1;
    e2 = z2 - s2;
    Z2 = e2/norm(e2);
    s3 = (z3'*Z2)*Z2 + (z3'*Z1)*Z1; 
    e3 = z3 - s3;
    Z3 = e3/norm(e3);
    s4 = (z4'*Z3)*Z3 + (z4'*Z2)*Z2 + (z4'*Z1)*Z1; 
    e4 = z4 - s4;
    Z4 = e4/norm(e4);
    Zall = [Z1 Z2 Z3 Z4];
    orthoZ = Zall*sqrt(obs);
end