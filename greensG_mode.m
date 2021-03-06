function G = greenstG_mode(psi,z,N_modes,modes,rho_w,zr,sources)
%% Inputs
if nargin<7
    zr=2:5:200;                    %(m) Receiver Depth
    zs=10:12:190;                  %(m) Source Depth
    rs=5000+(0:31)*100;            %(m) Source Range
    trs = (0*zs'+1)*rs;
    tzs = zs'*(0*rs+1);
    sources = [trs(:) tzs(:)]';
end;
lzr = length(zr);              %    Corresponding Lengths
lf = size(psi,3);

M = size(sources,2);

%% Green's Function Calculation
G = zeros(lzr,M,lf);
for ii=1:lf
    p_temp=0;
    for kk=1:N_modes(ii)
        psi_zs=interp1(z,psi(:,kk,ii),sources(2,:));
        psi_zr=interp1(z,psi(:,kk,ii),zr);
        p_temp = p_temp + (j./  (rho_w.*4).*psi_zr'*  (psi_zs.*besselh(0,1,modes(kk,ii).*sources(1,:))) );  %Pressure(r,z)
    end
    G(:,:,ii)=(p_temp);
    % GGp = GGp + p_temp'*Gg((0:lzr-1)*lf+ii);
end
