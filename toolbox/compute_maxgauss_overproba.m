function p  = compute_maxgauss_overproba( q, s, N )

% compute_maxgauss_overproba - empirical cumulative probability of maximum of Gaussians
%
%   p  = compute_maxgauss_overproba( q, s, N );
%
%   p(i) is an estimate (using N samples, so large N imply good precision)
%   that the maximum of q absolute value of unit Gaussian variable are
%   smaller than s(i).
%
%   Copyright (c) 2010 Gabriel Peyré

global gaussmax_P;
global gaussmax_V;
global gaussmax_Vmax;
global gaussmax_Vmin;

if isempty(gaussmax_P)
    %% DO NOT USE LOOKUP %%
    
    U = max(abs(randn(q,N)));
    p = zeros(length(s),1);
    for i=1:length(s)
        p(i) = sum(U<=s(i))/N;
    end
    
else
    %% USE LOOKUP %%
    p = -ones(length(s),1);
    p(s<gaussmax_Vmin(q)) = 0;
    p(s>gaussmax_Vmax(q)) = 1;    
    p(p<0) = interp1( gaussmax_V{q}, gaussmax_P{q}, s(p<0) );
end