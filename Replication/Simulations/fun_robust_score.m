function kp = fun_kpcal(z,z2,pihat,uhat)
    xhat = z*pihat;
    mxz = z2 - xhat*(xhat\z2);
    kp = uhat'*mxz*(((mxz.*(uhat.^2))'*mxz)\(mxz'*uhat));
end