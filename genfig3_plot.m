%%% GENERATE FIGURE 3 OF THE PAPER %%%
% Probability of support inclusion as a function of T/gamma0 for k=kbeta
% Copyright (c) 2011 Gabriel Peyré

rep = 'results/noise-dep/';
if not(exist(rep))
    mkdir(rep);
end
    
rho = .02;
clf; hold on;
h = plot(Tlist/Tcrit, Q/i, 'k');  set(h, 'LineWidth', 2);
% h = plot([1 1], [-rho 1+rho], 'k--'); set(h, 'LineWidth', 2); 
axis([min(Tlist/Tcrit) max(Tlist/Tcrit) -rho 1+rho]);
box on;
set(gca, 'FontSize', 20);
% xlabel('T/\gamma');
saveas(gcf, [rep 'T-dependence' '-n' num2str(n) '-p' num2str(p) '-a' num2str(a) '-b' num2str(b) '.eps'], 'eps');
        


