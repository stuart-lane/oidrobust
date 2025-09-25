function kp = fun_kpcal_nw(z,z2,pihat,uhat,L)
    xhat = z*pihat;
    mxz = z2 - xhat*(xhat\z2);
    varnw = NeweyWest(uhat,mxz,L,0);
    kp = uhat'*mxz*inv(varnw)*mxz'*uhat;
end