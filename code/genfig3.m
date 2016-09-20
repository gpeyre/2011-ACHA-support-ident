%%% GENERATE FIGURE 3 OF THE PAPER %%%
% Probability of support inclusion as a function of T/gamma0 for k=kbeta
% Copyright (c) 2011 Gabriel Peyré

%%
% test for the support recovery when the signal (T) is very large, 
% and for a varying gamma/sigma
% where sigma=|w|/sqrt(n)

addpath('toolbox/');

% dimensions
n = 1000*8;
p = 4*n;

a = .9; b = .9;
a = .8; b = .8;

k = round( a*b*n/(2*log(p)) );


%%
% critical value for gamma
% sigma=1

ntest = 100;
gamma = sqrt( 2*log(p)/(1-a) );
Tcrit = gamma;
Tlist = linspace(.5, 3, ntest) * Tcrit;


renew_rate = 100;

Q = zeros(ntest,1);
i = 0;
while true
    i = i+1;    
    progressbar(mod(i,500)+1,500);    
    % random matrix
    AI = randn(n,k)/sqrt(n);
    BI = AI'*AI;
    % random signs
    xI = sign(randn(k,1));
    x1 = BI\xI;        
    % random noise and its projection
    w = randn(n,1); w = sqrt(n) * w/norm(w);
    w1 = BI\( AI'*w );
    % solutions for various T
    xstar = xI*Tlist + repmat( - gamma*x1 + w1, [1 ntest]);
    % test condition (C1)
    q = sum( abs( sign(xstar) - repmat(xI, [1 ntest]) ) )==0;            
    % increment counter
    Q = Q + q';
    if mod(i,500)==0
        % disp(['** Saved - ' num2str(i)]);
        Qratio = Q/i;
        % save result into file
        save(['backup-noisedep-n' num2str(n) '-p' num2str(p)], 'Qratio', 'gammas_list', 'a', 'b', 'k');
    end
    
end 
