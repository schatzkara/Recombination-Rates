u = 0:0.001:1; % mutation rate (u = alpha + beta)
k = 0:0.1:100; % ratio of transitions to transversions (k = alpha/beta)
phi = 0.5;
pABC = (1-u).^2; % A=B=C
pABc = u.^2.*(k.^2 + phi^2 + (1-phi)^2)./(1+k).^2; % A = B =/= C
paBC = u.*(1-u); % A =/= B = C
pAbC = u.*(1-u); % B =/= A = C
pabc = 1- (1-u).^2 - 2.*u.*(1-u) - u.*(k.^2 + phi^2 + (1-phi)^2)./(1+k).^2;

[U, K] = meshgrid(u, k);
ratio = U.^2.*(K.^2 + phi^2 + (1-phi)^2)./(1+K).^2 ./ (2*U - U.^2);

mesh(U, K, ratio);
xlabel('U (mutation rates)')
ylabel('K (transition/transversions)')
zlabel('Convergent Mutations/SNPS')