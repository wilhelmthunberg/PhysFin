close all
clear all
disp('-------------PhysFin-HA2-----------')
format short
data=dlmread('ex2_3.dat');
s =data(:,1);
y = data(:,2);


A = [s.^3, s.^2 , s, s.^0 ]; %define A matrix

coeff = inv(A'*A)*A'*y; %find coefficents, eqution 7.4


p3 = coeff(1).*s.^3 +coeff(2).*s.^2 +coeff(3).*s +coeff(4); % polynomial

rms = sqrt((1/length(s))* sum((p3-y).^2)); % root mean square to find error bars
lambda = diag(s.^0 .* (1/rms)); %lambda
C= inv((A'*lambda^2*A)); %covariance matrix, equation7.7

disp('3rd degree polyfit coefficients:')
disp(coeff')

disp('rms:')
disp(rms)
coeff_unc =sqrt(diag(C)); %coeff unc are diag entries of C

disp ('Covariance matrix:')
disp(C)
disp('3rd degree polyfit coefficient uncetainties:')
disp(coeff_unc')



figure(1)
plot(s,p3) 
hold on
plot(data(:,1), data(:,2), '*')
legend( '3rd degree fit', 'scatter')

%From coeff and coeff uncertainties it appears that the coefficient for the linear term may be reoved
q=4; 
p=3;
m = length(y)-q; %length of y is number of data points
n=1;
chi_q = sum (((y- p3)/rms).^2)  %equation 7.3

p3_p = coeff(1).*s.^3 +coeff(2).*s.^2  +coeff(4);



A_p = [s.^3, s.^2 , s.^0 ];
coeff_p = inv(A_p'*A_p)*A_p'*y;
p3_p = s.^3.*coeff_p(1)+ s.^2.*coeff_p(2) +coeff_p(3);

chi_p = sum (((y- p3_p)/rms).^2)
F_value = ((chi_p - chi_q)/(q-p))/(chi_q/m) %equation 7.37

xhat=n*F_value/(m+n*F_value);
p=1-betainc(xhat,n/2,m/2)
disp('--------------------------------------------------')