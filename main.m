clc;
clear;
close all;

%% variables
cont = 0;

for t=1:100
    x = randn(1000,1) + 1i*randn(1000,1);
    n = length(x);
    
    if exist('Params')                == 0,  Params.n           = n;        end
    if isfield(Params, 'L')           == 0,  Params.L           = 2.8;      end
    if isfield(Params, 'T')           == 0,  Params.T           = 500;      end
    if isfield(Params, 'r')           == 0,  Params.r           = 2;        end
    if isfield(Params, 'e')           == 0,  Params.e           = 10^-10;   end
    if isfield(Params, 'y1')          == 0,  Params.y1          = 0.5;      end
    if isfield(Params, 'u0')          == 0,  Params.u0          = 45;       end
    if isfield(Params, 'y')           == 0,  Params.y           = 0.01;     end
    if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 250;      end
    if isfield(Params, 'alpha')       == 0,  Params.alpha       = 0.5;      end
    
    L           = Params.L;
    m           = n*L;
    Params.m    = m;
    
%     display(Params)
    
    Amatrix = (randn(m,n) + 1i*randn(m,n))/sqrt(2);
    
    % Make linear operators;
    A = @(I)  Amatrix*I;
    At = @(I) Amatrix'*I;
    
    %% Make signal and data (noiseless)
    y = abs(A(x));
    
    f = @(I,u) (1/m)*sum((sqrt(abs(Amatrix*I).^2+u^2)-y).^2);
    
    %% main loop
    tic
    [z0,z,Relerrs] = PRSF(x,y,Params, A, At,Amatrix,f);
    toc
    
    if min(Relerrs) <= 1e-5
        cont = cont  + 1;
    end
    T = length(Relerrs);
    fprintf('Percentage: %f, iter: %d, error: %f \n',cont/t,t, Relerrs(1));
end

%% results

fprintf('Relative error after initialization: %f\n', Relerrs(1))
fprintf('Relative error after %d iterations: %f\n', T, Relerrs(end))

figure, semilogy(0:T-1,Relerrs)
xlabel('Iteration'), ylabel('Relative error (log10)'), ...
title('Relative error vs. iteration count')
