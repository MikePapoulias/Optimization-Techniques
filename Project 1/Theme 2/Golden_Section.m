format long
syms x;

f1 = (x-2)^2 + x*log(x+3);
f2 = 5^x + (2-cos(x))^2;
f3 = exp(x)*(x^3-1)+(x-1)*sin(x);

% Midenismoi twn pinakwn
X_1 = 0; X_2 = 0; X_3 = 0; 
Y_1 = 0; Y_2 = 0; Y_3 = 0;
A_1_1 = 0; A_1_2 = 0; A_1_3 = 0; A_2_1 = 0; A_2_2 = 0; A_2_3 = 0; A_3_1 = 0; A_3_2 = 0; A_3_3 = 0;


fprintf("\nGolden Section for f1:")
[a, b, k, num_of_calcul] = xrisos_tomeas(f1, -1, 3, 0.01)

fprintf("\nGolden Section for f2:")
[a, b, k, num_of_calcul] = xrisos_tomeas(f2, -1, 3, 0.01)

fprintf("\nGolden Section for f3:")
[a, b, k, num_of_calcul] = xrisos_tomeas(f3, -1, 3, 0.01)



% First figure
figure("Name",sprintf("Golden Section Change parameter"))

% Metaboli ypologismwn f1 gia metablito l
for i=4:50
    X_1(1,i) = 0.001*i;
end
for i=4:50
    [a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f1, -1, 3, X_1(1,i));
    Y_1(1,i) = num_of_calcul;
end
subplot(1,3,1);
plot(X_1,Y_1)
xlim([0.004 0.05])
title("Golden Section f1")
xlabel("accuracy l")
ylabel("subs(f1(x))")


% Metaboli ypologismwn f2 gia metablito l
for i=4:50
    X_2(1,i) = 0.001*i;
end
for i=4:50
    [a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f2, -1, 3, X_2(1,i));
    Y_2(1,i) = num_of_calcul;
end
subplot(1,3,2);
plot(X_2,Y_2)
xlim([0.004 0.05])
title("Golden Section f2")
xlabel("accuracy l")
ylabel("subs(f2(x))")



% Metaboli ypologismwn f3 gia metablito l
for i=4:50
    X_3(1,i) = 0.001*i;
end
for i=4:50
    [a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f3, -1, 3, X_3(1,i));
    Y_3(1,i) = num_of_calcul;
end
subplot(1,3,3);
plot(X_3,Y_3)
xlim([0.004 0.05])
title("Golden Section f3")
xlabel("accuracy l")
ylabel("subs(f3(x))")




% Second figure
figure("Name",sprintf("Golden Section space limits"))

% Metaboli diastimatos [a, b] tis f1 gia l=0.005
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f1, -1, 3, 0.005);
for i=1:k
    A_1_1(1,i) = i;
end
subplot(3,3,1);
plot(A_1_1,Ka,"*r")
hold on
plot(A_1_1,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f1 (l=0.005)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f1 gia l=0.05
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f1, -1, 3, 0.05);
for i=1:k
    A_1_2(1,i) = i;
end
subplot(3,3,4);
plot(A_1_2,Ka,"*r")
hold on
plot(A_1_2,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f1 (l=0.05)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f1 gia l=0.5
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f1, -1, 3, 0.5);
for i=1:k
    A_1_3(1,i) = i;
end
subplot(3,3,7);
plot(A_1_3,Ka,"*r")
hold on
plot(A_1_3,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f1 (l=0.5)")
xlabel("k")
ylabel("[a,b]")



% Metaboli diastimatos [a, b] tis f2 gia l=0.005
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f2, -1, 3, 0.005);
for i=1:k
    A_2_1(1,i) = i;
end
subplot(3,3,2);
plot(A_2_1,Ka,"*r")
hold on
plot(A_2_1,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f2 (l=0.005)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f2 gia l=0.05
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f2, -1, 3, 0.05);
for i=1:k
    A_2_2(1,i) = i;
end
subplot(3,3,5);
plot(A_2_2,Ka,"*r")
hold on
plot(A_2_2,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f2 (l=0.05)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f2 gia l=0.5
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f2, -1, 3, 0.5);
for i=1:k
    A_2_3(1,i) = i;
end
subplot(3,3,8);
plot(A_2_3,Ka,"*r")
hold on
plot(A_2_3,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f2 (l=0.5)")
xlabel("k")
ylabel("[a,b]")



% Metaboli diastimatos [a, b] tis f3 gia l=0.005
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f3, -1, 3, 0.005);
for i=1:k
    A_3_1(1,i) = i;
end
subplot(3,3,3);
plot(A_3_1,Ka,"*r")
hold on
plot(A_3_1,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f3 (l=0.005)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f3 gia l=0.05
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f3, -1, 3, 0.05);
for i=1:k
    A_3_2(1,i) = i;
end
subplot(3,3,6);
plot(A_3_2,Ka,"*r")
hold on
plot(A_3_2,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f3 (l=0.05)")
xlabel("k")
ylabel("[a,b]")

% Metaboli diastimatos [a, b] tis f3 gia l=0.5
[a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f3, -1, 3, 0.5);
for i=1:k
    A_3_3(1,i) = i;
end
subplot(3,3,9);
plot(A_3_3,Ka,"*r")
hold on
plot(A_3_3,Kb,"ob")
xlim([1 20])
ylim([-2 4])
title("Spaces for f3 (l=0.5)")
xlabel("k")
ylabel("[a,b]")




function [a, b, k, num_of_calcul, Ka, Kb] = xrisos_tomeas(f, a, b, l)
x1 = a + (1-0.618)*(b-a);     %initialization
x2 = a + 0.618*(b-a);      %initialization
num_of_calcul = 0;       %number of calculations of f
y1 = subs(f,{x1}); 
num_of_calcul = num_of_calcul + 1;    %calcul += 1
y2 = subs(f,{x2}); 
num_of_calcul = num_of_calcul + 1;    %calcul += 1
Ka(1,1) = a;
ia = 1;
Kb(1,1) = b;
ib = 1;
k = 1;
while(true)
    ia = ia+1;
    ib = ib+1;
    if(b-a<l)
        break
    else
        if(y1 > y2)
            a = x1;
            Ka(1,ia) = a;
            Kb(1,ia) = Kb(1,ia-1);
            x1 = x2;
            y1 = y2;
            x2 = a + 0.618*(b-a);
            y2 = subs(f,{x2});  
            num_of_calcul = num_of_calcul + 1;   %calcul+=1
            k = k + 1;
        else
            b = x2;
            Kb(1,ib) = b;
            Ka(1,ib) = Ka(1,ib-1);
            x2 = x1;
            y2 = y1;
            x1 = a + (1-0.618)*(b-a);
            y1 = subs(f,{x1});
            num_of_calcul = num_of_calcul + 1;   %calcul+=1
            k = k + 1;
        end
    end
end
end