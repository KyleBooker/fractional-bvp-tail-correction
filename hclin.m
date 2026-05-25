% hclin.m -- post-processing script for output produced by homoclinic.c
% (or any of the *_derivative.c scheme variants).
%
% Loads homoclinic.in to recover the run parameters, loads homoclinic.out
% to get the computed u(y) on the FFT mesh, overlays the analytic homoclinic
% (3a) and ground state (3b) reference solutions, and fits a power law to
% the tail to estimate the decay constant a(p,gamma) of eq. (4).
%
% Reference: Booker & Nec (2019), Math. Model. Nat. Phenom.

clear all

prm=load('homoclinic.in');
mt=prm(1); g=prm(3); p=prm(4); n=prm(5);
N=prm(9); eps=prm(11); tau0=prm(12);
u=load('homoclinic.out');
y=2*pi*n*((0:N-1)-N/2)'/(eps*N);
un=((p+1)/2)^(1/(p-1))*sech((p-1)/2*y).^(2/(p-1));
uv=2./(1+y.^2);

ly=log(y);
i_tail=round(0.75*N):N;
p1=polyfit(ly(i_tail),log(uv(i_tail)),1)
tail=polyval(p1,ly(i_tail));

%plot(y,u,'k',y,uv,'r--')
plot(ly,log(u),'k:',ly(i_tail),tail,'r--')
