function [N_modes,kr,Uw,z,idx]=normal_modes(f);
%
%[Num of Modes,Normalized Mode Shape,Depth]=normal_modes(freq(Hz),h max depth);
%
% %%%%%%%%%%%%%%%%%%%%%
% %                   %
% %   Shallow Water   %
% %    Normal Mode    %
% % (Shooting Method) %
% %        by         %
% %  Shaun Anderson   %
% %                   %
% %%%%%%%%%%%%%%%%%%%%%
%
%
global h rho_w rho_s R_rho cs cz h
  N_modes=[];
   kr=[];
   Uw=[];
   z=[];
   idx=[];
%
%Find Mode Shapes
%
options=optimset('Display','off','TolX',1e-6); %%%SET ROOT FINDER TOLERENCES

kw=2*pi*f/cz;
ks=2*pi*f/cs;
if ks < kw
    kn=linspace(ks,kw,150);%[ks:step:kw kw];
%     kn=(1-[0 logspace(-3,0,200)])*(kw-ks)+ks;
else
    kn=linspace(kw,ks,150);%[ks:step:kw kw];
% 	kn=(1-[0 logspace(-3,0,200)])*(ks-kw)+kw;
end
% step=kn(end)-kn(end-1);
kr=[];
for i=1:length(kn)
    try [res(i),fval(i),flag(i),option]=fzero(@(kr)shoot_mode(f,kr,h,R_rho,cs,cz),[kn(i) kn(i+1)],options);
    catch res(i)=NaN;
    end
end
res(isnan(res)) = [];
kr=unique(res*10^7)/10^7;%Possible Area to improve Accuracy%
kr(find(kr>kw | kr<ks))=0;
kr=fliplr(kr);
clear i idx id kn res U
N_modes=length(kr);
x=0:.5:h;
idx2=length(x);
%   
%Normalize
%
for mode=1:N_modes
    [uo,U_temp,z_temp]=shoot_mode(f,kr(mode),h,R_rho,cs,cz); %Calculate Mode Shape In Water
    idx(mode)=length(U_temp);
    U(1:idx(mode),mode)=(U_temp);
    z(1:idx(mode),mode)=z_temp;
    gamma=sqrt(kr(mode)^2-(2*pi*f/cs)^2);
    Us(1:idx2,mode)=U(1,mode)/exp(-gamma*h)*exp(-gamma*(x+h)); %Calculate Mode Shape in Sediment
    normalize=1/(trapz(flipud(z(1:idx(mode),mode)),flipud(U(1:idx(mode),mode).^2))+1/R_rho*trapz(x,Us(1:idx2,mode).^2));
    Uw(1:idx(mode),mode)=U(1:idx(mode),mode)*normalize.^.5;
    Usd(1:idx2,mode)=Us(1:idx2,mode)*normalize.^.5;
if Uw(idx(mode)-1,mode)<0;
        Uw(:,mode)=-Uw(:,mode);Usd(:,mode)=-Usd(:,mode);
    end
end
for mode=1:N_modes
    Uw(:,mode)=interp1(z(1:idx(mode),mode),Uw(1:idx(mode),mode),linspace(h,0,length(z)));
end
z=linspace(0,h,length(z));