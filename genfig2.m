%%% GENERATE FIGURE 2 OF THE PAPER %%%
% Probability of support recovery for large T as a function of gamma/gamma0 for k=kbeta.
% Copyright (c) 2011 Gabriel Peyré

%%
% test for the support recovery when the signal (T) is very large, 
% and for a varying gamma/sigma
% where sigma=|w|/sqrt(n)

addpath('toolbox/');

%%
% Setup paramters.

n = 1000*8;
p = 4*n;

a = .9; b = .9;
a = .8; b = .8;

k = round( a*b*n/(2*log(p)) );


%%
% critical value for gamma/sigma

ntest = 20;
gammas_crit = sqrt( 2*log(p)/(1-a) );
gammas_list = linspace(.5,1.5, ntest)*gammas_crit;

renew_rate = 100;

Q = zeros(ntest,1);
i = 0;
while true
    i = i+1;
    
    progressbar(mod(i,500)+1,500);
    
    % random matrix
    AI = randn(n,k)/sqrt(n);
    BI = AI'*AI;
    % random noise
    w = randn(n,1); w = sqrt(n) * w/norm(w);
    % random signs
    sx = sign(randn(k,1));
    d0 = AI*(BI\sx);
    % w1 = BI\( AI'*w ); % P_{A_I}(w)
    Pw = w-AI*( BI\( AI'*w ) );
    if mod(i,renew_rate)==1
        AIc = randn(n,p-k)/sqrt(n);
    end
    
    % test condition (C2)
    % forall j \notin I,  |<a_j, gamma*d0 + P_{A_I^bot}(w)>| < gamma
    C = AIc' * ( d0*gammas_list + repmat(Pw, [1 ntest]) );
    q = max(abs(C)) <= gammas_list;
        
    % increment counter
    Q = Q + q';
    if mod(i,500)==0
        % disp(['** Saved - ' num2str(i)]);
        Qratio = Q/i;
        % save result into file
        save(['backup-noisedep-n' num2str(n) '-p' num2str(p)], 'Qratio', 'gammas_list', 'a', 'b', 'k');
    end
    
end 

