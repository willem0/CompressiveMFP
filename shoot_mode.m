function [u0,u,z]=shoot_mode(freq,kr,h,R_rho,cs,cz)
%
%Setup Profile Parameters
%
%R_rho=1.5;%1.957;%Ratio of Rho_s / Rho_w
%cs=1600;%1.145*c_z(h);%base Speed of Sound in Sediment (m/s)
%
%Shooting Solution
%
w=freq*2*pi;
%Initial Conditions
kn=kr;
mu=sqrt(kn^2-(w/cs)^2);
u_s=1;
dudz=-u_s/R_rho*mu;
options = odeset('AbsTol',1e-6,'RelTol',1e-7); %%%SET ODE SOLVER TOLERENCES
[z,R]=ode113(@(z,R) normal(z,R,kn,w,cz),[h 0],[u_s dudz],options);
u0=R(end,1);
u=R(:,1);
end
%%
function res = normal(z,R,kn,w,cz)
%Unpack Vector
u=R(1);
Q=R(2);
%Coupled ODE's
dudz=Q;
dQdz=(kn^2-(w/(cz*ones(size(z))))^2)*u;
%Repack vector
res=[dudz;dQdz];
end