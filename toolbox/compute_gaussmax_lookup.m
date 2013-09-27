%% 
% Compute lookup tables for the max of Gaussians

% number of trials
lookup_number = 5000;
% number of sample of the proba function
lookup_precision = 100;

pmax = 2000;

%% 
% gaussmax_P(i) is the proba of the max being <= gaussmax_V(i)
global gaussmax_P;
global gaussmax_V;
global gaussmax_Vmax;
global gaussmax_Vmin;
gaussmax_P = {};
gaussmax_V = {};
gaussmax_Vmin = [];
gaussmax_Vmax = [];

disp('#### Building lookup tables for the probability of max of Gaussian, might take a while ###');

for ik=1:length(k_list)
    k = k_list(ik);
    progressbar(ik,length(k_list));
    
    q = p-k;
    % Maximum of abs(qgaussians), by slices of pmax values
    U = zeros(1,lookup_number);
    np = ceil(q/pmax);
    for ip=1:np
        if ip<np
            nbr = pmax;
        else
            nbr = q-(np-1)*pmax;
        end
        U = max(U, max( abs(randn(pmax,lookup_number)) ) );        
    end
    [h,t] = hist(U,lookup_precision); h = h/sum(h);
    I = find(h>=1e-3);
    V = linspace(t(I(1)), t(I(end)), lookup_precision);
    for i=1:lookup_precision
        P(i) = sum(U<=V(i))/lookup_number;
    end
    gaussmax_P{p-k} = P;
    gaussmax_V{p-k} = V;
    gaussmax_Vmin(p-k) = min(V);
    gaussmax_Vmax(p-k) = max(V);
end