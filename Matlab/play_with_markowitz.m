% play_with_markowitz.m

clear all
close('all')

N=1000;  % number of trading days
r1=0.01*ones(N,1)+0.01*randn(N,1);      %..make some test data
r2=0.02*ones(N,1)+0.02*randn(N,1)-0.3*r1;
r3=0.03*ones(N,1)+0.03*randn(N,1)+0.3*r1-0.5*r2;
%mm=1:N; plot(mm,r1,mm,r2)
ee=ones(3,1);   % three stocks
ra=[mean(r1); mean(r2) ; mean(r3)];   % average return

rp=[r1-mean(r1),r2-mean(r2),r3-mean(r3)]; % for the covariance matrix
C=rp'*rp/N; % covariance matrix

%load ra
%load C

CC=inv(C);  % and its inverse

%.....................the 2x2 matrix to determine the lambdas
A=[ra'*CC*ra , ee'*CC*ra ;
    ra'*CC*ee , ee'*CC*ee]
AA=inv(A);

K=0;
for rho=0:0.0002:0.04;  %...loop over the desired returns
    lambda=AA*[rho;1];
    w=lambda(1)*CC*ra+lambda(2)*CC*ee;    % eq. 5.12
    K=K+1;
    rrho(K)=rho;
    sig(K)=sqrt(w'*C*w);    % definition of sigma
end

plot(sig,rrho)
xlabel('Volatility \sigma')
ylabel('Portfolio return \rho')

hold on
rhomin=-AA(1,2)/AA(1,1);
sigmin=sqrt(det(AA)/AA(1,1))     %= sigmin2=sqrt(1/A(2,2))
plot(sigmin,rhomin,'*')
text(sigmin+0.001,rhomin,'\leftarrow Minimum risk')


%...add a zero risk asset
r0=0.008;  % rate of return of the zero risk asset
deltar=ra-r0*ee;
K=0;
for rho=r0:0.0001:5*r0
    ww=(rho-r0)/(deltar'*CC*deltar)*CC*deltar;  % eq.5.22
    K=K+1;
    rrho2(K)=rho;
    sig2(K)=sqrt(ww'*C*ww);
end
plot(sig2,rrho2,'r')

rhot=r0+(deltar'*CC*deltar)/(ee'*CC*deltar)
sigmat=sqrt((deltar'*CC*deltar)/(ee'*CC*deltar)^2)
plot(sigmat,rhot,'r*')
text(sigmat+0.001,rhot-0.001,'\leftarrow Tangent/Market portfolio')

text(0.031,0.027,'\leftarrow Efficient frontier');
text(0.015,0.033,'Capital market line \rightarrow')


plot(sqrt(C(1,1)),ra(1),'k*',sqrt(C(2,2)),ra(2),'k*',sqrt(C(3,3)),ra(3),'k*')


