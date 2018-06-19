k = 0:0.5:4; % Values of kappa = transitions/transversions
u = 0:0.001:1;
for i= 1:length(k)
% pABC(i) = 1- (1-u).^2 - 2.*u.*(1-u) - u.*(k(i).^2 + phi^2 + (1-phi)^2)./(1+k(i)).^2;
pABc(i,:) = u.^2.*(k(i).^2 + phi^2 + (1-phi)^2)./(1+k(i)).^2;
end
pABC = (1-u).^2;
SNPS = 1-pABC;

for(i = 1:length(k))
   plot(SNPS, pABc(i,:)); 
   hold on;
end
hold off;
xlabel("SNPS/Length")
ylabel("Convergent mutations/Length")
title(sprintf("Convergent Mutations/L vs SNPs/L for Kappa = {0, 0.5, 1, ... 4}"))

for(i = 1:length(k))
   labels(i) = sprintf("k = %g", k(i));
end

legend(labels);

% k = 0.01:0.01:101;
% pABc = (k.^2 + 1/2)./(k+1).^2;
% %pABc = 1 - 2./(k+1) - 0.75./(k+1).^2;
% plot(k, pABc);