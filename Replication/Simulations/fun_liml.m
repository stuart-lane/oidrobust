function [bhat,pihat,uhat,vhat] = fun_liml(y,x,z,izz,n,kx)
    w = [y x];
    sigw = (w'*w/n);
    wpzw = (w'*z)*(z\w);
    mine = min(eig(sigw\(wpzw)))/n;

    bhat = (wpzw(2:kx+1,2:kx+1)-mine*(x'*x))\(wpzw(2:kx+1,1)-mine*(x'*y));
    uhat = y-x*bhat;
    H = z - uhat*(uhat\z);
    pihat = H\x;
    vhat = x - z*pihat;
end