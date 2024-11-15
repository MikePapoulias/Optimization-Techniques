syms x y
format long


f = (1/3)*(x^2)+3*(y^2);


X_pinakas = 0;




% Thema_2:
fprintf("\nSteepest Descent Projection for [5 -5], gk = 0.5, sk = 5 and e = 0.01:")
[X k, X_pinakas] = megisti_kathodos_proboli(f, [5 -5], 0.01, 5, 0.5, 600);
X
k
% Arxiko simeio: [5 -5] - Steepest Descent Projection
figure("Name", sprintf("Megisti Kathodos Proboli_[5 -5]"))
fsurf(f, [-15 15 -15 15])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
% Ektypwsi tou synolou periorismou:
hold on
x = [-10 -10 5 5 -10]; 
y = [12 -8 -8 12 12];
plot(x,y,'r', 'LineWidth',2)
saveas(gcf, 'Figure_Steepest Descent Proboli[5 -5]_siglisi.fig')




% Thema_3:
fprintf("\nSteepest Descent Projection for [-5 10], gk = 0.1, sk = 15 and e = 0.01:")
[X k, X_pinakas] = megisti_kathodos_proboli(f, [-5 10], 0.01, 15, 0.1, 600);
X
k
% Arxiko simeio: [-5 10] - Steepest Descent Projection
figure("Name", sprintf("Megisti Kathodos Proboli_[-5 10]"))
fsurf(f, [-15 15 -15 15])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
% Ektypwsi tou synolou periorismou:
hold on
x = [-10 -10 5 5 -10]; 
y = [12 -8 -8 12 12];
plot(x,y,'r', 'LineWidth',2)
saveas(gcf, 'Figure_Steepest Descent Proboli[-5 10]_siglisi.fig')





% Thema_4:
fprintf("\nSteepest Descent Projection for [8 -10], gk = 0.2, sk = 0.1 and e = 0.01:")
[X k, X_pinakas] = megisti_kathodos_proboli(f, [8 -10], 0.01, 0.1, 0.2, 600);
X
k
% Arxiko simeio: [8 -10] - Steepest Descent Projection
figure("Name", sprintf("Megisti Kathodos Proboli_[8 -10]"))
fsurf(f, [-15 15 -15 15])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
% Ektypwsi tou synolou periorismou:
hold on
x = [-10 -10 5 5 -10]; 
y = [12 -8 -8 12 12];
plot(x,y,'r', 'LineWidth',2)
saveas(gcf, 'Figure_Steepest Descent Proboli[8 -10]_siglisi.fig')











function [X, k, X_pinakas] = megisti_kathodos_proboli(f, X, e, sk, gk, limit)
syms x y
k = 1;
% Ypologismos probolis 1hs sunistwsas symeiou ekkinisis
if(X(1) < -10)
    X(1) = -10;
elseif(X(1) > 5)
    X(1) = 5;
end
% Ypologismos probolis 1hs sunistwsas symeiou ekkinisis
if(X(2) < -8)
    X(2) = -8;
elseif(X(2) > 12)
    X(2) = 12;
end
X_pinakas(k,:) = X;
while(true)
    n = vpa(subs(jacobian(f), {x,y}, {X}));
    if(norm(n) < e)
        break
    elseif(k >= limit)
        fprintf("\nO algorithmos den syglinei gia o sygkekrimeno sk.")
        break
    else
        dk = -n;
        X_arg = vpa(X + sk*dk);
        % Ypologismos probolis 1hs sunistwsas tou xk
        if(X_arg(1) < -10)
            X_arg(1) = -10;
        elseif(X_arg(1) > 5)
            X_arg(1) = 5;
        end
        X_bar(1) = X_arg(1);
        % Ypologismos probolis 2hs sunistwsas tou xk
        if(X_arg(2) < -8)
            X_arg(2) = -8;
        elseif(X_arg(2) > 12)
            X_arg(2) = 12;
        end
        X_bar(2) = X_arg(2);
        X = vpa(X + gk*(X_bar - X));
        k = k+1;
        X_pinakas(k,:) = X;
    end
end
end

