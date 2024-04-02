
clc
clear
r_source=[25,30,0]; % source position, units mm
r_det=[55,30,0];

x_local_center=40; % local region where the activation is
y_local_center=30;
z_local_center=30;

xspan=40;
yspan=40;
zspan=20;
r_span=[xspan,yspan,zspan];
r_act=[x_local_center, y_local_center,z_local_center];% the size of the local activation


act=0.3; %stength of activation, rDb-1
lambda=800*10^(-6);% mm

musp=1; %1/mm
mua=0;
Db=10^(-6); %mm^2
vRBC=0;
alpha=1; % 
tau=[0,10.^([-6:0.1:0])];


Grsrdt0=CalculateG12Pt(r_source+[0,0,1/musp],r_det,musp,mua,tau, Db,vRBC,alpha,lambda);

% Rytov approximation

voxelsize=1; % mm^3
%Radius_activation=15;


deltaDb=Db*act; % the activation strength

[Phi]=CalculatePhis2(r_source,r_det,musp,mua,deltaDb,Db,r_act,r_span,voxelsize,lambda,tau);
Phi=Phi./Grsrdt0;

Grsrdt1=Grsrdt0.*exp(Phi); % results from Rytov

figure
semilogx(tau,abs(Grsrdt0)./abs(Grsrdt0(1)),'b');hold on
semilogx(tau,abs(Grsrdt1)./abs(Grsrdt1(1)),'r');

legend('baseline','activation')
ylabel('g_1(\tau)')
xlabel('\tau (s)')

g1=abs(Grsrdt0)./abs(Grsrdt0(1));
g1_act_Rytov=abs(Grsrdt1)./abs(Grsrdt1(1));


% covert to K2
betaval=1;
nn=0;
for ii=2:length(tau)
 nn=nn+1;
    K2base(nn)=2*betaval./tau(ii).*trapz(tau(1:ii),g1(1:ii).^2.*(1-tau(1:ii)./tau(ii))); %baseline
  Texp=tau(ii);  
% fun = @(tauc)beta*tauc/Texp.*(1+tauc./2/Texp.*(exp(-2*Texp./tauc)-1))-K2base(nn);
% x0 = [1*10^(-5)];
% x = fsolve(fun,x0);
% tauc_base(nn)=x;
K2act_Rytov(nn)=2*betaval./tau(ii).*trapz(tau(1:ii),g1_act_Rytov(1:ii).^2.*(1-tau(1:ii)./tau(ii)));% Rytov

end
figure
semilogx(tau(2:end),K2base);hold on
semilogx(tau(2:end),K2act_Rytov);

legend('baseline','activation')
ylabel('K^2')
xlabel('T_{exp} (s)')

