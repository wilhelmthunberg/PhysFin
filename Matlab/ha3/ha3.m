clear all
close all
disp('----PhysFin-HA3------------------')

disp('----------------Ex1--------------')
d=importdata('ex3_1.dat'); 
t=d(:,1); co2=d(:,2); month=d(:,3);
p=polyfit(t,co2,1); 
trend=polyval(p,t);
co2_nt = co2-trend; %co2 notrend
figure(1)
plot(t,co2_nt)
title('CO2 levels without trend')
ylabel('CO2')
xlabel('time')

co2_nt_ns = co2_nt(13:end)-co2_nt(1:end-12); %notrend no seasonality
figure(2)
plot(t(13:end),co2_nt_ns)
title('CO2 levels without trend and seasonality')
ylabel('CO2')
xlabel('time')
ind_shift = [0:length(co2_nt_ns)-1];

for i=ind_shift
gamma(i+1) = 1/length(co2_nt_ns) * sum(co2_nt_ns(1+i:end).*co2_nt_ns(1:end-i));

end
figure(3)
rhos = gamma/gamma(1);
bar(ind_shift,rhos)
hold on
plot(ind_shift,ones(1,length(ind_shift))*(2/sqrt(length(ind_shift))),'r' )
plot(ind_shift,-ones(1,length(ind_shift))*(2/sqrt(length(ind_shift))),'r' )
title('Autocorrelations')
hold off

%PACF
ps=[1:4];
for p=ps
    rho = rhos(1:p+1) ;
    rho_mat =[rho(1:end-1)'] ;
    for i = 1:p-1
        rho_mat=[rho_mat circshift(rho((1:end-1)),i)'];
    end
    ph = inv(rho_mat)*rho(2:end)';
    phi_jj(p) = ph(end);
end
figure(4)
bar(ps,phi_jj)
hold on
plot(ps,ones(1,length(phi_jj))*(2/sqrt(length(ind_shift))),'r' )
plot(ps,-ones(1,length(phi_jj))*(2/sqrt(length(ind_shift))),'r' )
title('Partial Autocorrelation')
inds = [(nonzeros((phi_jj >2/sqrt(length(ind_shift))).*[1:length(phi_jj)]))]';
phi_jj = nonzeros(phi_jj.*(phi_jj >2/sqrt(length(ind_shift))));

A=[];
y=co2_nt_ns;
for i=[1:length(inds)]
A = [A y(inds(end)+1-inds(i):end-inds(i))];
end

disp('Coefficients')
phi = inv(A'*A)*A' * y(inds(end)+1:end)

new_y = A*phi;
figure(5)
plot(t(13:end),co2_nt_ns,'*')
title('fit vs data')
hold on
plot(t(13+max(inds):end),new_y)
legend('data','fit')
hold off


disp('----------------Ex2--------------')
disp(' ')
disp(' ')
clear all
D= @(x) 1/(3*gamma(3)).*exp(-x.^(1/3)); %distribution function
C = @(x) gammainc(x.^(1/3),3); %cumulative
Xi =@(x,n) n*C(x).^(n-1) .* D(x); 
n=[0 1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10];

for i = [1:length(n)]
xn(i) = integral(@(x)(x.*Xi(x,n(i))),0,inf); %use built integrate function
end

figure(6)
semilogx(n,xn.^(1/3))
title('The average maximum expected value')
ylabel('maximum expected value')
xlabel('n')

disp('---finished----------------------')


