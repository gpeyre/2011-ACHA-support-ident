%%% GENERATE FIGURE 1 OF THE PAPER %%%
% Probability of sparsistency as a function of k (the sparsity)
% Copyright (c) 2011 Gabriel Peyré

%%
% Analyze and display the result for support identification.

n = 1000*8;
p = 4*n;

%% number of processor

nproc = 2;

Q1 = [];
k_svg1 = []; sigma_svg1 = [];
for proc=1:nproc
    load(['backup-n' num2str(n) '-p' num2str(p) '-proc' num2str(proc)]);
    Q1 = [Q1 Q];
    k_svg1 = [k_svg1 k_svg];
    sigma_svg1 = [sigma_svg1 sigma_svg];
end

%q = min([size(Q,2) length(k_svg) length(sigma_svg) ]);
%k_svg = k_svg(1:q);
%sigma_svg = sigma_svg(1:q);
%Q = Q(:,1:q);

b_list = [1 .9 .8 .7];

%k_list = unique(k_svg);
%sigma_list = unique(sigma_svg);

rho = .02;

for isigma=1:length(sigma_list)
    sigma = sigma_list(isigma);

    R = [];
    for k=k_list(:)'
        I = find(k_svg1==k & sigma_svg1==sigma);
        if length(I)>0
            R(:,end+1) = sum(Q1(:,I),2)/length(I);
        else
            R(:,end+1) = zeros(size(Q1,1), 1);
        end
    end

    lgd = {};
    for i=1:length(gamma_list)
        lgd{end+1} = ['\gamma=T/' num2str(round(1/gamma_list(i)),2)];
    end

    gamma_disp = [2 4 6 7];
    gamma_disp = [1 2];
    gr = {'-' '--' '.-' '-.'};
    ms = [15 10 10 10]; 
    clf; hold on;
    for ig=1:length(gamma_disp)
        r = R(gamma_disp(ig),:);
        % r = (r + r([1 1:end-1]) + r([2:end end]))/3;
        h = plot(k_list, r, ['k' gr{ig}]);
        set(h, 'MarkerSize', ms(ig));
        set(h, 'LineWidth', 2);
    end    
    legend({lgd{gamma_disp}});        
    kdisp = a_list(isigma)*b_list/2*n/log(p);
    col = {'' '--' '-.'  ':' '*-'};
    for ik=1:length(b_list)       
        h = plot(kdisp(ik)*[1 1], [-rho 1+rho], ['k' col{ik}]);
        set(h, 'LineWidth', 2);
    end
    axis([min(k_list) max(k_list) -rho 1+rho]);
    box on;
    set(gca, 'FontSize', 20);
    
    rep = 'results/';
    if not(exist(rep))
        mkdir(rep);
    end
    saveas(gcf, [rep 'identifiability-n' num2str(n) '-p' num2str(p) '-a' num2str(a_list(isigma),2) '.png'], 'png');
    saveas(gcf, [rep 'identifiability-n' num2str(n) '-p' num2str(p) '-a' num2str(a_list(isigma),2) '.eps'], 'eps');
end
