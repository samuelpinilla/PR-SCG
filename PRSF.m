%% Implementation of the Smoothing Conjugate Gradient Phase Retrieval Method

function [z0,xk1,Relerrs] = PRSF(x,y, Params, A, At,Amatrix,f)    
    %% Initialization
    y       = y.^2;
    z0      = randn(Params.n, 1);
    z       = z0 / norm(z0, 'fro');    % Initial guess

    normest = sqrt(sum(y(:)) / numel(y(:)));    % Estimate norm to scale eigenvector
    m       = Params.m;
    ymag    = sqrt(y);

    [ysort, ~] = sort(y, 'ascend');
    ythresh    = ysort(round(m / 1.2));
    ind        = y >= ythresh;

    Aselect    = Amatrix(ind, :);
    weights    = (ymag(ind)).^(Params.alpha); % weights w_i

    for tt = 1:Params.npower_iter                   % Power iterations
        z  = Aselect' * (weights .* (Aselect * z));
        z  = z / norm(z, 'fro');
    end

    z = normest * z;
    Relerrs = norm(x - exp(-1i * angle(trace(x' * z))) * z, 'fro') / norm(x, 'fro'); % Initial rel. error   

    %% t=0
    xk = z;
    uk = Params.u0;
    gk = compute_grad(xk, uk, ymag, Params, A, At);
    dk = -gk;
    
    %% Loop
    for t = 1: Params.T
        alpha = backtracking(xk,gk,dk,uk,f);
        xk1 = xk + alpha*dk;                                                            % Gradient update
        
        Relerrs = [Relerrs, norm(x - exp(-1i*angle(trace(x'*xk1))) * xk1, 'fro')/norm(x,'fro')]; % Relative error
        
        [dk,uk,xk,gk] = stepUpdate(xk1,xk,gk,uk,dk,Params,ymag,A,At);                      % update
        
        if min(Relerrs)<=1e-5
            break;
        end
    end
