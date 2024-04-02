function Grsrdt=CalculateG12Pt(r_source,r_det,musp,mua,tau, Db,vRBC,alpha,lambda) %% Boas 2016

%rs=r_source+[0,0,1/musp]; % location of the point source

rs=r_source;

rsi=[rs(1),rs(2),-rs(3)]-[0,0,4/(3*musp)]; % - the third component, location of the image source 

c=3*10^11; %speed of light in the medium mm/s
n=1.33; %refractive index of tissue;
v=c/n; % speed of light in the medium. we probably don't need this but I kept it here.
D=v/(3*musp);% photon diffusion coefficient
K0=2*pi/lambda; % in mm^-1;
%K2=3.*mua.*musp+musp.^2.*K0.^2.*alpha.*(6*Db.*tau+vRBC.^2.*tau.^2);
K2=v*(mua+1/3*musp*K0.^2.*alpha.*(6*Db.*tau+vRBC.^2.*tau.^2))/D;
%Grsrd=v/(4*pi*D)*(exp(-sqrt(3*musp*mua)*norm(rs-r_det))/norm(rs-r_det)- exp(-sqrt(3*musp*mua)*norm(rsi-r_det))/norm(rsi-r_det));
%Grsrdt=(3*musp/(4*pi)*(exp(-sqrt(K2).*norm(rs-r_det))./norm(rs-r_det)- exp(-sqrt(K2).*norm(rsi-r_det))./norm(rsi-r_det)));
Grsrdt=abs(3*musp/(4*pi)*(exp(-sqrt(K2).*norm(rs-r_det))./norm(rs-r_det)- exp(-sqrt(K2).*norm(rsi-r_det))./norm(rsi-r_det)));

end




