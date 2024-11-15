format long
syms x;

f1 = (x-2)^2 + x*log(x+3);
f2 = 5^x + (2-cos(x))^2;
f3 = exp(x)*(x^3-1)+(x-1)*sin(x);

% Midenismoi twn pinakwn
X_1_1 = 0; X_1_2 = 0; X_1_3 = 0; X_2_1 = 0; X_2_2 = 0; X_2_3 = 0; X_3_1 = 0; X_3_2 = 0; X_3_3 = 0;
Y_1_1 = 0; Y_1_2 = 0; Y_1_3 = 0; Y_2_1 = 0; Y_2_2 = 0; Y_2_3 = 0; Y_3_1 = 0; Y_3_2 = 0; Y_3_3 = 0;
A_1_1 = 0; A_1_2 = 0; A_1_3 = 0; A_2_1 = 0; A_2_2 = 0; A_2_3 = 0; A_3_1 = 0; A_3_2 = 0; A_3_3 = 0;


fprintf("Bisector Method for f1:")
[a, b, k, num_of_calcul] = dixotomos(f1, -1, 3, 0.001, 0.01)

fprintf("\nBisector Method for f2:")
[a, b, k, num_of_calcul] = dixotomos(f2, -1, 3, 0.001, 0.01)

fprintf("\nBisector Method for f3:")
[a, b, k, num_of_calcul] = dixotomos(f3, -1, 3, 0.001, 0.01)



% First figure
figure("Name",sprintf("Bisector Method Change parameters"))

% Metaboli ypologismwn f1 gia l=0.01
for i=1:48
    X_1_1(1,i) = 0.0001*i;
end
for i=1:48
    [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f1, -1, 3, X_1_1(1,i), 0.01);
    Y_1_1(1,i) = num_of_calcul;
end
subplot(2,3,1);
plot(X_1_1,Y_1_1)
xlim([0.0001 0.0048])
title("Bisector Method f1 (l=0.01)")
xlabel("distance ε")
ylabel("subs(f1(x))")

% Metaboli ypologismwn f1 gia e=0.001
for i=4:50
    X_1_2(1,i) = 0.001*i;
end
for i=4:50
    [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f1, -1, 3, 0.001, X_1_2(1,i));
    Y_1_2(1,i) = num_of_calcul;
end
subplot(2,3,4);
plot(X_1_2,Y_1_2)
xlim([0.004 0.05])
title("Bisector Method f1 (ε=0.001)")
xlabel("accuracy l")
ylabel("subs(f1(x))")



% Metaboli ypologismwn f2 gia l=0.01
for i=1:48
    X_2_1(1,i) = 0.0001*i;
end
for i=1:48
    [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f2, -1, 3, X_2_1(1,i), 0.01);
    Y_2_1(1,i) = num_of_calcul;
end
subplot(2,3,2);
plot(X_2_1,Y_2_1)
xlim([0.0001 0.0048])
title("Bisector Method f2 (l=0.01)")
xlabel("distance ε")
ylabel("subs(f2(x))")

% Metaboli ypologismwn f2 gia e=0.001
for i=4:50
    X_2_2(1,i) = 0.001*i;
end
for i=4:50
    [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f2, -1, 3, 0.001, X_2_2(1,i));
    Y_2_2(1,i) = num_of_calcul;
end
subplot(2,3,5);
plot(X_2_2,Y_2_2)
xlim([0.004 0.05])
title("Bisector Method f2 (ε=0.001)")
xlabel("accuracy l")
ylabel("subs(f2(x))")



% Metaboli ypologismwn f3 gia l=0.01
for i=1:48
    X_3_1(1,i) = 0.0001*i;
end
for i=1:48
    [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f3, -1, 3, X_3_1(1,i), 0.01);
    Y_3_1(1,i) = num_of_calcul;
end
subplot(2,3,3);
plot(X_3_1,Y_3_1)
xlim([0.0001 0.0048])
title("Bisector Method f3 (l=0.01)")
xlabel("distance ε")
ylabel("subs(f3(x))")

% Metaboli ypologismwn f3 gia e=0.001
for i=4:50
    X_3_2(1,i) = 0.001*i;
end
for i=4:50
    [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f3, -1, 3, 0.001, X_3_2(1,i));
    Y_3_2(1,i) = num_of_calcul;
end
subplot(2,3,6);
plot(X_3_2,Y_3_2)
xlim([0.004 0.05])
title("Bisector Method f3 (ε=0.001)")
xlabel("accuracy l")
ylabel("subs(f3(x))")





% Second figure
figure("Name",sprintf("Bisector Method space limits"))

% Metaboli diastimatos [a, b] tis f1 gia l=0.005
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f1, -1, 3, 0.001, 0.005);
for i=1:k
    A_1_1(1,i) = i;
end
subplot(3,3,1);
plot(A_1_1,Ka,"*r")
hold on
plot(A_1_1,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f1 (l=0.005)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f1 gia l=0.05
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f1, -1, 3, 0.001, 0.05);
for i=1:k
    A_1_2(1,i) = i;
end
subplot(3,3,4);
plot(A_1_2,Ka,"*r")
hold on
plot(A_1_2,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f1 (l=0.05)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f1 gia l=0.5
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f1, -1, 3, 0.001, 0.5);
for i=1:k
    A_1_3(1,i) = i;
end
subplot(3,3,7);
plot(A_1_3,Ka,"*r")
hold on
plot(A_1_3,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f1 (l=0.5)")
xlabel("k")
ylabel("[a,b]")



% Metaboli diastimatos [a, b] tis f2 gia l=0.005
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f2, -1, 3, 0.001, 0.005);
for i=1:k
    A_2_1(1,i) = i;
end
subplot(3,3,2);
plot(A_2_1,Ka,"*r")
hold on
plot(A_2_1,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f2 (l=0.005)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f2 gia l=0.05
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f2, -1, 3, 0.001, 0.05);
for i=1:k
    A_2_2(1,i) = i;
end
subplot(3,3,5);
plot(A_2_2,Ka,"*r")
hold on
plot(A_2_2,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f2 (l=0.05)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f2 gia l=0.5
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f2, -1, 3, 0.001, 0.5);
for i=1:k
    A_2_3(1,i) = i;
end
subplot(3,3,8);
plot(A_2_3,Ka,"*r")
hold on
plot(A_2_3,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f2 (l=0.5)")
xlabel("k")
ylabel("[a,b]")



% Metaboli diastimatos [a, b] tis f3 gia l=0.005
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f3, -1, 3, 0.001, 0.005);
for i=1:k
    A_3_1(1,i) = i;
end
subplot(3,3,3);
plot(A_3_1,Ka,"*r")
hold on
plot(A_3_1,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f3 (l=0.005)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f3 gia l=0.05
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f3, -1, 3, 0.001, 0.05);
for i=1:k
    A_3_2(1,i) = i;
end
subplot(3,3,6);
plot(A_3_2,Ka,"*r")
hold on
plot(A_3_2,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f3 (l=0.05)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f3 gia l=0.5
[a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f3, -1, 3, 0.001, 0.5);
for i=1:k
    A_3_3(1,i) = i;
end
subplot(3,3,9);
plot(A_3_3,Ka,"*r")
hold on
plot(A_3_3,Kb,"ob")
xlim([1 15])
ylim([-2 4])
title("Spaces for f3 (l=0.5)")
xlabel("k")
ylabel("[a,b]")




function [a, b, k, num_of_calcul, Ka, Kb] = dixotomos(f, a, b, e, l)
k = 1;
Ka(1,1) = a;
ia = 1;
Kb(1,1) = b;
ib = 1;
num_of_calcul = 0;       %number of calculations of f
while(true)
    ia = ia+1;
    ib = ib+1;
    if(b-a<l)
        break
    else
        x1 = (a+b)/2 - e;
        x2 = (a+b)/2 + e;
        if(subs(f,{x1}) < subs(f,{x2}))
            b = x2;
            Kb(1,ib) = b;
            Ka(1,ib) = Ka(1,ib-1);
            num_of_calcul = num_of_calcul + 2;    %calcul += 2
        else
            a = x1;
            Ka(1,ia) = a;
            Kb(1,ia) = Kb(1,ia-1);
            num_of_calcul = num_of_calcul + 2;    %calcul += 2
        end
    end
    k = k + 1;
end
end