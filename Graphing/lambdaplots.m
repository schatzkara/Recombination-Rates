mu = 0.0001:0.0001:0.1;
L = 100:100:100000; 
L = L';
lambdas = (L*mu)';

qmat = arrayfun(@(x) poissinv(0.99999999,x),lambdas);

%contour plot
figure;
[C,h] = contourf(L,mu,qmat,25);
clabel(C,h)
colormap(cool)
colorbar
title("Poisson critical values given mu and L for 99.999999% Significance")
xlabel("L")
ylabel("Mu")