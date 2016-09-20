%%% GENERATE FIGURE 2 OF THE PAPER %%%
% Probability of support recovery for large T as a function of gamma/gamma0 for k=kbeta.
% Copyright (c) 2011 Gabriel Peyré

Qratio = Q/i;

if 0
    str = ['backup-noisedep-n' num2str(n) '-p' num2str(p)];
    load(str);
end

rep = 'results/noise-dep/';
if not(exist(rep))
    mkdir(rep);
end
    
rho = .02;
clf; hold on;
h = plot(gammas_list/gammas_crit, Qratio, 'k');  set(h, 'LineWidth', 2);
h = plot([1 1], [-rho 1+rho], 'k--'); set(h, 'LineWidth', 2); 
axis([min(gammas_list/gammas_crit) max(gammas_list/gammas_crit) -rho 1+rho]);
box on;
set(gca, 'FontSize', 20);
saveas(gcf, [rep 'gamma-dependence' '-n' num2str(n) '-p' num2str(p) '-a' num2str(a) '-b' num2str(b) '.eps'], 'eps');
        
