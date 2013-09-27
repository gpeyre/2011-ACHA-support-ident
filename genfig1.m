%%% GENERATE FIGURE 1 OF THE PAPER %%%
% Probability of sparsistency as a function of k (the sparsity)
% Copyright (c) 2011 Gabriel Peyré


path(path, 'toolbox/');

noise_type = 'uniform';
noise_type = 'gaussian';

if not(exist('processeur'))
    processeur = 1;
end

%% 
% dimensions

% setup 1, (n,p) = (8000,32000)
n = 1000*8;
p = 4*n;
% setup 2, (n,p) = (4000,36000)
n = 3000;
p = 12*n;

options.noise_type = noise_type;
options.computation_mode = 'slow';
options.computation_mode = 'fast';
options.maxgauss_precision = 1000;

if not(exist('processeur'))
    processeur = 1;
end

% gamma values
gamma_list = 1./[1.8 2 3 4 6 8 10]';
gamma_list = 1./[4 6]';

% a and b values
a_list = [.6 .8 .9 .95];
a_list = [.9 .95];
bmin = .5; bmax = 1.8;

% sigma values
sigma_list = 1/6 * sqrt( (1-a_list)/(2*log(p)) );

% k values
kmin_list = a_list*bmin/2*n/log(p);
kmax_list = a_list*bmax/2*n/log(p);
kmin = floor(min(kmin_list));
kmax = ceil(max(kmax_list));

if n<=1000
    k_list = kmin:1:kmax;
elseif n<=2000
    k_list = kmin:2:kmax;
else
    k_list = kmin:3:kmax;
end
    
if not(exist('Q'))
    Q = [];
    k_svg = [];
    sigma_svg = [];
    i = 0;
else
    q = min([size(Q,2) length(k_svg) length(sigma_svg)]);
    k_svg = k_svg(1:q);
    sigma_svg = sigma_svg(1:q);
    Q = Q(:,1:q);
    i = q;
end

if strcmp(options.computation_mode, 'fast') && not(exist('gaussmax_P'))
    compute_gaussmax_lookup;
end

[Sigma_list,K_list] = meshgrid(sigma_list,k_list);

% infinite loop
while true
    i = i+1;
    
    % choose sigma and k
    i1 = mod(i,length(Sigma_list(:)))+1;
    sigma = Sigma_list(i1);
    k = K_list(i1);
    
    % check for identifiability
    Q(:,end+1) = compute_support_identifiability(n,p,k,sigma,gamma_list, options);
    sigma_svg(end+1) = sigma;
    k_svg(end+1) = k;
    
    if mod(i,500)==0
        disp(['** Saved - ' num2str(i)]);
        % save result into file
        save(['backup-n' num2str(n) '-p' num2str(p) '-proc' num2str(processeur)], 'Q', 'sigma_svg', 'k_svg', 'k_list', 'a_list', 'gamma_list', 'sigma_list');
    end
end
