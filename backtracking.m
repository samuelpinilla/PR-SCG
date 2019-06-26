function alpha = backtracking(z,gk,dk,uk,f)
    alpha = 1;
    c = 0.1;
    p = 0.4;
    while f(z+alpha*dk,uk) > f(z,uk)+c*alpha*real(gk'*dk)
       alpha = p*alpha; 
    end
end