function [Phi]=CalculatePhis2(r_source,r_det,musp,mua,deltaDb,Db,r_act,r_span,voxelsize,lambda,tau)

c=3*10^11; %speed of light in the medium mm/s
n=1.33; %refractive index of tissue;
v=c/n; % speed of light in the medium. we probably don't need this but I kept it here.

K0=2*pi/lambda; % in mm^-1;
Phi=0;
%xregion=[r_act(1)-r_span(1):voxelsize:r_act(1)+r_span(1)];
%yregion=[r_act(2)-r_span(2):voxelsize:r_act(2)+r_span(2)];
%zregion=[r_act(3)-r_span(3):voxelsize:r_act(3)+r_span(3)];

xregion=[r_act(1)-r_span(1)+voxelsize:voxelsize:r_act(1)+r_span(1)];
yregion=[r_act(2)-r_span(2)+voxelsize:voxelsize:r_act(2)+r_span(2)];
zregion=[r_act(3)-r_span(3)+voxelsize:voxelsize:r_act(3)+r_span(3)];
Dr0=v/(3*musp);
%Phi_pert=0;

for ii=1:length(xregion)
    for jj=1:length(yregion)
        for kk=1:length(zregion)
            rnow=[xregion(ii),yregion(jj),zregion(kk)];
            %if norm(rnow-r_act)<=Radius_activation
                Phi=Phi-voxelsize^3*(deltaDb*2*musp.*tau*K0^2*v/Dr0).*CalculateG12Pt(rnow,r_det,musp,mua,tau, Db,0,1,lambda)...
                    .*CalculateG12Pt(r_source+[0,0,1/musp],rnow,musp,mua,tau,Db,0,1,lambda)/v*Dr0;
            %end
            
        end
    end
end


end