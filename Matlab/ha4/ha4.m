 close all
 clear all
disp('-----Excercise 1-----')

f = @(x) log(x).*exp(-sqrt(x)); % function to be avaluated
lims = [1,pi];
Inum = integral(f, lims(1),lims(2)) %built in numerical sol

Ns = linspace(1,5*1e5+1, 1001);

for i = [1:length(Ns)]
x= lims(1) + (lims(2)-lims(1))*rand(1,Ns(i));
Imc(i) = ((lims(2)-lims(1))/Ns(i)) * sum(f(x));
if i>1
err(i-1) =(Imc(i)-Imc(i-1)) ;
end
end
figure(1)
hold on
plot(Ns,Imc)
plot(Ns, ones(length(Ns))*Inum,'r')
ylim([0.3,0.35])
title('Integral')
xlabel('N')
legend('Monte-Carlo', 'Correct solution')
hold off

figure(2)
bar(Ns(2:end),err)
ylim([-0.05,0.05])
hold on
plot(Ns, ones(length(Ns))*1e-3,'r')
plot(Ns, -ones(length(Ns))*1e-3,'r')
title('Error of MC integral')
legend('error', '1e-3','-1e-3')
hold off
for i = 1:length(err)
    if err(i:end)<1e-3
        disp('Error completely stabalizes below 1e-3 for N larger than:')
        disp(Ns(i+1))
        break
    end
end
find(abs(err)< 1e-3);

disp(' ')
disp('-----Excercise 2-----')
del =0.3;%delta
lam = 2;%lambda
th=0.7;%theta
sig=0.5;%sigma - reinvestment fraction
beta =0.8;% beta-discount factor
%Initial conditions
k0 =1;
k(1) = (1-del)*k0;
y(1) = k(1)^(th)*lam;
i(1) = sig*y(1);%reinvest
p(1) = y(1)-i(1); %profit
v(1) = beta*log10(p(1));%joy 
months=60;
% base on IC, calculate the rest
for m = 2:months
    k(m) = (1-del)*k(m-1) + i(m-1);
    y(m) = k(m)^(th)*lam;
    i(m) = sig*y(m);
    p(m) = y(m)-i(m);
    v(m) = beta^(m)*log10(p(m));
end
Vp = sum(v); %cumulative joy


figure(3)
plot(1:months, p)
title('Profit p for first 60 months if \sigma = 0.5')
xlabel('months')
ylabel('p')


disp('b)---')
disp('Cumulative joy: ')
disp(Vp)

disp('c)---')
V=[];
sigmas =linspace(0,1,10000);
%find optimal sigma using for-loop
for sig =sigmas 
    del =0.3;
    lam = 2;
    th=0.7;
    beta =0.8;
    k0 =1;
    k(1) = (1-del)*k0;
    y(1) = k(1)^th*lam;
    i(1) = sig*y(1);
    p(1) = y(1)-i(1);
    v(1) = beta*log10(p(1)); 
    for m = 2:60
        k(m) = (1-del)*k(m-1)+i(m-1);
        y(m) = k(m)^th*lam;
        i(m) = sig*y(m);
        p(m) = y(m)-i(m);
        v(m) = beta^(m)*log10(p(m));
    end
V = [V, sum(v)];
 
end
figure(4)
plot(sigmas,V)
title('Cumulative joy as afunction of \sigma')
xlabel('\sigma')
ylabel('')
disp('Max cumulative joy:')
disp(max(V))
disp('at sigma: ')
disp(sigmas(V==max(V)))

