syms x y

f = (x^5)*exp(-x^2-y^2);

% function representation
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
saveas(gcf, 'Figure_1.fig')


% Χρησιμοποιούμε την fminsearch για να βρούμε το μέγιστο της συνάρτησης

