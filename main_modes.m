% %%%%%%%%%%%%%%%%%%%%%
% %                   %
% %   Vector Sensor   %
% %    Normal Mode    %
% %        by         %
% %  Shaun Anderson   %
% %                   %
% %%%%%%%%%%%%%%%%%%%%%

clear; tic;

%%% Simplified to Isovelocity
global h rho_w rho_s R_rho cs cz h
%
% Model Setup
%

h=200;    %m               %(m) Sediment Depth
%
% Environment Variables
%
rho_w=1000;        %kg/m^3 %Density of Water
cz=1520;           %m/s
R_rho=1.5;         %rho_s/rho_w
rho_s=rho_w*R_rho; %Density of Sediment
cs=1600;           %Speed of Sediment


fmin=101;%Hz
fmax=200;%Hz
df=1;%Hz
freq=fmin:df:fmax;

%%%Calcs
psi=zeros(length(0:h),3,length(freq));
wb = waitbar(0,'running simulation...');
for ii=1:length(freq)
    waitbar(ii/length(freq),wb);
    
    [N_modes(ii),kn,Uw,z,idx]=normal_modes(freq(ii));
    modes(1:length(kn),ii)=kn;
    for kk=1:N_modes(ii)
        psi(:,kk,ii)=interp1(z,Uw(:,kk),0:h);
    end
end
clear Uw idx kn

close(wb);
z=0:h;
N_modes=N_modes-[0 diff(N_modes)];%% Skip 1 places for Modes near cutoff Freq

toc, beep
save state psi;
save(['states/state psi']);

disp(['saving states/state_' num2str(fmin) '_' num2str(fmax) '_' num2str(cz) ' psi z N_modes modes rho_w freq']);
eval(['save   states/state_' num2str(fmin) '_' num2str(fmax) '_' num2str(cz) ' psi z N_modes modes rho_w freq']);
