% get_stocks.m

clear all; close all;
format short
q=dlmread('stocks.dat');
td=q(:,1);   % trading days
S1=q(:,2);   % Apple, 
S2=q(:,3);   % SP500
S3=q(:,4);   % Cola
hold on

%a)
figure(1)
set(gca,'FontSize',15) 
plot(td,S1/S1(1))
plot(td,S2/S2(1))
plot(td,S3/S3(3))
title('The relative evolvement of the stocks')
legend('Apple','SP500', 'Coca Cola' )
ylabel('value relative to index')
xlabel('day')


%b)
rj1 = (S1(2:end) - S1(1:end-1))./ S1(1:end-1);
rj2 = (S2(2:end) - S2(1:end-1))./ S2(1:end-1);
rj3 = (S3(2:end) - S3(1:end-1))./S3(1:end-1);
annual = [(1+mean(rj1))^365-1, (1+mean(rj2))^365-1,(1+mean(rj3))^365-1];

figure(2)
set(gca,'FontSize',15) 
plot(td(2:end),[rj1,rj2,rj3])
title('day-to-day returns')
legend('Apple','SP500', 'Coca Cola' )
ylabel('return relative previous day')
xlabel('day')

 disp('   Apple mean day-to-day return: ')
 disp(mean(rj1))
 disp('   SP500 mean day-to-day return: ')
 disp(mean(rj2))
 disp('   Coca-Cola mean day-to-day return: ')
 disp(mean(rj3))

ee=ones(3,1);
vec=[0;0;1];
ra=[mean(rj1); mean(rj2) ; mean(rj3)];  
rp=[rj1-mean(rj1),rj2-mean(rj2),rj3-mean(rj3)];

C=rp'*rp/length(rp);
CC = inv(C);

A=[ra'*CC*ra , ee'*CC*ra ,vec'*CC*ra;
    ra'*CC*ee , ee'*CC*ee,vec'*CC*ee;
   ra'*CC*vec , ee'*CC*vec, vec'*CC*vec ]; %new A vector with updated theory
AA=inv(A);

rhos = [1.05,1.1,1.15].^(1/(365))-1;%day-to-day return 
K=0;
 for rho=rhos;  %...loop over the desired returns
     lambda=AA*[rho;1;0.3];
     w=lambda(1)*CC*ra+lambda(2)*CC*ee + lambda(3)*CC*vec;    % eq. 5.12
     K=K+1;
     rrho(K)=rho;
     sig(K)=sqrt(w'*C*w);    % definition of sigma
    disp(w)
 end
figure(3)
plot(sig,rrho,'*', 'linewidth', 18)
xlabel('Volatility \sigma')
ylabel('Portfolio return \rho')
