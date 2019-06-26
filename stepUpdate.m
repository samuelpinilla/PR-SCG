function [dk1,uk1,xk1,gk1] = stepUpdate(xk1,xk,gk,uk,dk,Params,y,A,At)
    
    aux = compute_grad(xk1, uk, y, Params, A, At);
    r = Params.r;
    
    if norm(aux,'fro')>= Params.y*uk
        uk1 = uk;
    else
        uk1 = Params.y1*uk;
    end
    
    gk1 = compute_grad(xk1,uk1,y,Params,A,At);
    yk = gk1-gk;
    sk = xk1-xk;
    zk = yk+(Params.e*norm(gk1,'fro')^r+max(0,-real(sk'*yk/norm(sk,'fro')^2)))*sk;
    
    dk1 = -gk1+(real(gk1'*zk)/real(dk'*zk)-(2*norm(zk,'fro')^2*real(gk1'*dk)/real(dk'*zk)^2))*dk+(real(gk1'*dk)/real(dk'*zk))*zk;
    
end