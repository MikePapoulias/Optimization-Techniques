syms x y
format long

f = (x^5)*exp(-x^2-y^2);

X_pinakas = 0;




% Arxiko simeio: [-1 1] - Newton dosmeno gk
fprintf("Newton dosmeno gk_[-1 1]:")
[X, k, X_pinakas] = Niouton_dosmeno_gk(f, [-1,1], 0.0001);
X
k
figure("Name", sprintf("Niouton Dosmeno_gk_[-1 1]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton dosmeno gk[-1 1]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Dosmeno_gk f-k_[-1 1]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas(i,1),X_pinakas(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([1 k])
saveas(gcf, 'Figure_Newton dosmeno gk[-1 1]_f-k.fig')

% Arxiko simeio: [-1 1] - Newton Elaxistopoiisi f
fprintf("Newton Elaxistopoiisi f_[-1 1]:")
[X, k, X_pinakas] = Niouton_elaxistopoiisi_f(f, [-1,1], 0.0001);
X
k
figure("Name", sprintf("Niouton Elaxistopoiisi_f_[-1 1]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton Elaxistopoiisi_f [-1 1]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Elaxistopoiisi_f f-k_[-1 1]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas(i,1),X_pinakas(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([1 k])
saveas(gcf, 'Figure_Newton Elaxistopoiisi_f[-1 1]_f-k.fig')

% Arxiko simeio: [-1 1] - Newton Armijo
fprintf("Newton Armijo_[-1 1]:")
[X, k, X_pinakas] = Niouton_armijo(f, [-1,1], 0.0001);
X
k
figure("Name", sprintf("Niouton Armijo_[-1 1]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton Armijo [-1 1]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Armijo f-k_[-1 1]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas(i,1),X_pinakas(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([0 k+1])
saveas(gcf, 'Figure_Newton Armijo[-1 1]_f-k.fig')






% Arxiko simeio: [0 0] - Newton dosmeno gk
fprintf("Newton dosmeno gk_[0 0]:")
[X, k, X_pinakas] = Niouton_dosmeno_gk(f, [0,0], 0.0001);
X
k
figure("Name", sprintf("Niouton Dosmeno_gk_[0 0]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton dosmeno gk[0 0]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Dosmeno_gk f-k_[0 0]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas(i,1),X_pinakas(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([0 k+1])
saveas(gcf, 'Figure_Newton dosmeno gk[0 0]_f-k.fig')

% Arxiko simeio: [0 0] - Newton Elaxistopoiisi f
fprintf("Newton Elaxistopoiisi f_[0 0]:")
[X, k, X_pinakas] = Niouton_elaxistopoiisi_f(f, [0,0], 0.0001);
X
k
figure("Name", sprintf("Niouton Elaxistopoiisi_f_[0 0]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton Elaxistopoiisi_f [0 0]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Elaxistopoiisi_f f-k_[0 0]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas(i,1),X_pinakas(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([0 k+1])
saveas(gcf, 'Figure_Newton Elaxistopoiisi_f[0 0]_f-k.fig')

% Arxiko simeio: [0 0] - Newton Armijo
fprintf("Steepest Descent Armijo_[0 0]:")
[X, k, X_pinakas] = Niouton_armijo(f, [0,0], 0.0001);
X
k
figure("Name", sprintf("Niouton Armijo_[0 0]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas(i,1),X_pinakas(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton Armijo [0 0]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Armijo f-k_[0 0]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas(i,1),X_pinakas(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([0 k+1])
saveas(gcf, 'Figure_Newton Armijo[0 0]_f-k.fig')







% Arxiko simeio: [1 -1] - Newton dosmeno gk
fprintf("Newton dosmeno gk_[1 -1]:")
[X, k, X_pinakas_1] = Niouton_dosmeno_gk(f, [1,-1], 0.0001);
X
k
figure("Name", sprintf("Niouton Dosmeno_gk_[1 -1]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas_1(i,1),X_pinakas_1(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton dosmeno gk[1 -1]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Niouton Dosmeno_gk f-k_[1 -1]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas_1(i,1),X_pinakas_1(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([1 k])
saveas(gcf, 'Figure_Newton dosmeno gk[1 -1]_f-k.fig')

% Arxiko simeio: [1 -1] - Newton Elaxistopoiisi f
fprintf("Newton Elaxistopoiisi f_[1 -1]:")
[X, k, X_pinakas_1] = Niouton_elaxistopoiisi_f(f, [1,-1], 0.0001);
X
k
figure("Name", sprintf("Newton Elaxistopoiisi_f_[1 -1]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas_1(i,1),X_pinakas_1(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton Elaxistopoiisi_f [1 -1]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Newton Elaxistopoiisi_f f-k_[1 -1]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas_1(i,1),X_pinakas_1(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([1 k])
saveas(gcf, 'Figure_Newton Elaxistopoiisi_f[1 -1]_f-k.fig')

% Arxiko simeio: [1 -1] - Newton Armijo
fprintf("Newton Armijo_[1 -1]:")
[X, k, X_pinakas_1] = Niouton_armijo(f, [1,-1], 0.0001);
X
k
figure("Name", sprintf("Niouton Armijo_[1 -1]"))
fsurf(f, [-5 5 -5 5])
xlabel("x")
ylabel("y")
hold on
for(i=1:k)
    plot(X_pinakas_1(i,1),X_pinakas_1(i,2), 'ro')
end
saveas(gcf, 'Figure_Newton Armijo [1 -1]_siglisi.fig')

K = 0;
F_sub = 0;
figure("Name", sprintf("Newton Armijo f-k_[1 -1]"))
for(i=1:k)
    K(i) = i;
    F_sub(i) = subs(f, {x,y}, {X_pinakas_1(i,1),X_pinakas_1(i,2)});
end
plot(K,F_sub,"b-o")
xlabel("k")
ylabel("f(x,y)")
xlim([1 k])
saveas(gcf, 'Figure_Newton Armijo[1 -1]_f-k.fig')









function [X, k, X_pinakas] = Niouton_dosmeno_gk(f, X, e)
syms x y
apoklisi = 0;
k = 1;
X_pinakas(k,:) = X;
while(true)
    n = vpa(subs(jacobian(f), {x,y}, {X}));
    if(norm(n) < e)
        break
    else
        o = vpa(subs(jacobian(jacobian(f)),{x,y}, {X}));
        dk = -inv(o)*n';
        gk = 0.5;
        X_current = X;
        X = X + gk*dk';
        if(vpa(subs(f,{x,y},{X})) > vpa(subs(f,{x,y},{X_current})))
            apoklisi = apoklisi+1;
        elseif(vpa(subs(f,{x,y},{X})) == vpa(subs(f,{x,y},{X_current})))
            apoklisi = apoklisi+1;
            break
        end
        k = k+1;
        X_pinakas(k,:) = X;
    end
end
if(apoklisi > 0)
    fprintf(" (O algorithmos apoklinei)")
end
end





function [X, k, X_pinakas] = Niouton_elaxistopoiisi_f(f, X, e)
syms x y z
k = 1;
X_pinakas(k,:) = X;
apoklisi = 0;
while(true)
    n = vpa(subs(jacobian(f), {x,y}, {X}));
    if(norm(n) < e)
        break
    else
        o = vpa(subs(jacobian(jacobian(f)),{x,y}, {X}));
        dk = -inv(o)*n';
        gk = solve(vpa(subs(jacobian(f), {x,y}, {X+z*dk'})*(dk)) == 0);
        temp = 1/0;
        temp_gk = 0;
        for i=1:length(gk)
            if(gk(i) > 0)
                if(subs(f, {x,y}, {X+gk(i)*dk'}) < temp)
                    temp_gk = gk(i);
                    temp = subs(f, {x,y}, {X+temp_gk*dk'});
                end
            end
        end
        gk = temp_gk;
        X_current = X;
        X = X + gk*dk';
        if(vpa(subs(f,{x,y},{X})) > vpa(subs(f,{x,y},{X_current})))
            apoklisi = apoklisi+1;
        elseif(vpa(subs(f,{x,y},{X})) == vpa(subs(f,{x,y},{X_current})))
            apoklisi = apoklisi+1;
            break
        end
        k = k+1;
        X_pinakas(k,:) = X;
    end
end
if(apoklisi > 0)
    fprintf(" (O algorithmos apoklinei)")
end
end





function [X, k, X_pinakas] = Niouton_armijo(f, X, e)
syms x y
apoklisi = 0;
k = 1;
X_pinakas(k,:) = X;
a = 0.05;
b = 0.3;
s = 0.6;
while(true)
    n = vpa(subs(jacobian(f), {x,y}, {X}));
    if(norm(n) < e)
        break
    else
        o = vpa(subs(jacobian(jacobian(f)),{x,y}, {X}));
        dk = -inv(o)*n';
        mk = 0;
        gk = s*b^mk;
        while(true)
            if(vpa(subs(f, {x,y}, {X}) - subs(f, {x,y}, {X+gk*dk'})) >= vpa(-a*b^mk*s*(dk'*n')))
                break
            end
            mk = mk + 1;
            gk = s*b^mk;
        end
        X_current = X;
        X = X + gk*dk';
        if(vpa(subs(f,{x,y},{X})) > vpa(subs(f,{x,y},{X_current})))
            apoklisi = apoklisi+1;
        elseif(vpa(subs(f,{x,y},{X})) == vpa(subs(f,{x,y},{X_current})))
            apoklisi = apoklisi+1;
            break
        end
        k = k+1;
        X_pinakas(k,:) = X;
    end
end
if(apoklisi > 0)
    fprintf(" (O algorithmos apoklinei)")
end
end

