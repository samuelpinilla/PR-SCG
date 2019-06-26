%%  Compute the truncated gradient based on the Poisson log-likelihood function

function grad = compute_grad(z, u, y, Params, A, At)
    m = Params.m;
    
    yz = sqrt(abs(A(z)).^2+u^2);
    grad = 2/m*At(((yz-y)./yz).*A(z));
end
      
