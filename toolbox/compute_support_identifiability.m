function Q = compute_support_identifiability(n,p,k,sigma,gamma_list, options)

% compute_support_identifiability - test for L1 support identifiabiilty
%
%   Q = compute_support_identifiability(n,p,k,sigma,gamma, options);
%
%   Q=1 if sign(x*)=sign(x) where
%   
%   x* = argmin_{u} 1/2*|A*u-y|^2 + gamma*|u|_1  where y=A*x+w
%
%   and:
%       A is randn(n,p)/sqrt(n);
%       w is randnom with |w|=sigma/sqrt(n)
%       x(I)=sign(randn(k,1))
%
%   gamma can be a vector of value, and then Q(i) corresponds to the identifiability for gamma(i).
%
%   Set options.computation_mode='fast' to use the fast code.
%   In this setting, set options.maxgauss_precision as the number of
%       samples used to estimated the proba of the max of (p-k) Gaussian
%
%   Copyright (c) 2010 Gabriel Peyre, Jalal Fadili and Charles Dossal


options.null = 0;
noise_type = getoptions(options, 'noise_type', 'gaussian');
computation_mode = getoptions(options, 'computation_mode', 'fast');
maxgauss_precision = getoptions(options, 'maxgauss_precision', 1000);

%%
% Generate noise.

if strcmp(noise_type, 'gaussian')
    w = randn(n,1);
elseif strmpc(noise_type, 'uniform')
    w = rand(n,1) - .5;
end
w = sigma*sqrt(n) * w/norm(w);


%% 
% Generate signal

xI = sign(randn(k,1));
sx = xI;

%%
% Generate matrices.

AI = randn(n,k)/sqrt(n);
BI = AI'*AI;

%%
% Precompute some important vectors.

x1 = BI\sx;
w1 = BI\( AI'*w );

%%
% Check first condition.

ntests = length(gamma_list);
Q = zeros(ntests,1); 
r = zeros(ntests,1);
for i=1:ntests
    gamma = gamma_list(i);
    % the tentative solution should have same sign
    xstar = xI - gamma*x1 + w1;
    Q(i) = norm(sign(xstar) - sx)==0;
    % norm of redisual
    r(i) = norm( AI*xI+w-AI*xstar );
end


%%
% Check second condition.

if strcmp(computation_mode, 'slow')   
    AIc = randn(n,p-k)/sqrt(n); 
    for i=find(Q)'
        gamma = gamma_list(i);
        % the correlation with other atoms should be smaller than gamma
        gamma = gamma_list(i);
        cor = AIc' * (gamma*AI*x1 + w-AI*w1);
        Q(i) = max(abs(cor))<=gamma;
    end
else    
    % cor is a Gaussian variable of std  r(i)/sqrt(n).
    % test wether max(p-k unit gaussian)>sqrt(n)/r(i)
    Q(Q==1) = compute_maxgauss_overproba( p-k, sqrt(n).*gamma_list(Q==1)./r(Q==1), maxgauss_precision );
end

