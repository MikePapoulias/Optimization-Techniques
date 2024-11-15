syms x y
format long


f = (1/3)*(x^2)+3*(y^2);

fsurf(f, [-15 15 -15 15])
xlabel("x")
ylabel("y")
saveas(gcf, 'Figure_1.fig')




fprintf("Megisti kathodos gia gk = 0.1:")
[X, k] = megisti_kathodos_dosmeno_gk(f, [2 2], 0.001, 0.1);
X
k

fprintf("Megisti kathodos gia gk = 0.3:")
[X, k] = megisti_kathodos_dosmeno_gk(f, [2 2], 0.001, 0.3);
X
k

fprintf("Megisti kathodos gia gk = 3:")
[X, k] = megisti_kathodos_dosmeno_gk(f, [2 2], 0.001, 3);

fprintf("Megisti kathodos gia gk = 5:")
[X, k] = megisti_kathodos_dosmeno_gk(f, [2 2], 0.001, 5);







function [X, k] = megisti_kathodos_dosmeno_gk(f, X, e, gk)
syms x y
k = 1;
while(true)
    n = vpa(subs(jacobian(f), {x,y}, {X}));
    if(norm(n) < e)
        break
    else
        dk = -n;
        X = X + gk*dk;
        k = k+1;
    end
    if(k >= 600)
        fprintf("\nO algorithmos den sigklinei gia to sigkekrimeno gk.\n\n")
        break
    end
end
end