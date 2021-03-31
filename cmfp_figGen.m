% cmfp_figGen.m
clear;
close all;

set(0,'DefaultAxesFontsize',18);

phi_type = 'QR-Gaussian';
% phi_type = 'Gaussian';
% phi_type = 'Bernoulli';
pp = 1-linspace(1e-5,1-1e-5,1000);
shg;


%AMBIGUITY FUNCTIONS
scenariolist = {'single','coherent'};
for k0=1:length(scenariolist)
    scenario = scenariolist{k0};
load(['data/amb_' scenario '_' phi_type]);
cmap = colormap('hot'); cmap = flipud(cmap);
lzs = length(zs);
lrs = length(rs);
b = round(lzs/2)+1;
c = round(lrs/2)+1;

%this quantity hn is really |h_normalized(r)|^2, so we use 10 instead of 20
dbimn = 10*log10(hn /max(hn (:)));     dbimn(dbimn<-20) = -20;
dbim1 = 10*log10(ht1/max(ht1(:)));     dbim1(dbim1<-20) = -20;
dbim2 = 10*log10(ht2/max(ht2(:)));     dbim2(dbim2<-20) = -20;
colormap(cmap);
imagesc(rs,zs,dbimn);  colorbar;  clipboard('copy',['ambn_' scenario '_' phi_type]);  pause;
imagesc(rs,zs,dbim1);  colorbar;  clipboard('copy',['amb1_' scenario '_' phi_type]);  pause;
imagesc(rs,zs,dbim2);  colorbar;  clipboard('copy',['amb2_' scenario '_' phi_type]);  pause;

plot(zs,dbimn(:,b),zs,dbim1(:,b),zs,dbim2(:,b),...
    target(2)+rad(1)*[-1 -1 1 1],[0 -20 -20 0],'--', target(2)+rad(2)*[-1 -1 1 1],[0 -20 -20 0],'--');
legend({'nMFP',...
       ['cMFP: M=' num2str(Mlist(1))],...
       ['cMFP: M=' num2str(Mlist(2))],...
       'Target Ellipse',...
       'Mainlobe Ellipse'});
clipboard('copy',['depth_' scenario '_' phi_type]);
pause;

plot(rs,dbimn(c,:),rs,dbim1(c,:),rs,dbim2(c,:),...
    target(1)+rad(3)*[-1 -1 1 1],[0 -20 -20 0],'--', target(1)+rad(4)*[-1 -1 1 1],[0 -20 -20 0],'--');
legend({'nMFP',...
       ['cMFP: M=' num2str(Mlist(1))],...
       ['cMFP: M=' num2str(Mlist(2))],...
       'Target Ellipse',...
       'Mainlobe Ellipse'});
clipboard('copy',['range_' scenario '_' phi_type]);
pause;
end;
% coherent values, need all radii parameters in the data file



%TAILS
scenariolist = {'single','coherent','incoherent'};
for k0=1:length(scenariolist)
    scenario = scenariolist{k0};
    
    load(['data/tail_' scenario '_' phi_type]);
    ppm = 1-linspace(1e-5,1-1e-5,length(sm));
    ppc = 1-linspace(1e-5,1-1e-5,length(scm));
%     semilogy(snm,pp,'--',sm,pp,'--',scm,1-linspace(0,1,size(scm,1)));
    semilogy(sm,ppm,'--',snm,ppm,'--',scm,ppc);
    axis([0 1 .01 1]);
    clipboard('copy',['tail_' scenario '_' phi_type]);
    
%     Rewrite this to depend on Mlist
    if strcmp(scenario,'single'),      legend({'MFP','nMFP','cMFP: M=3','cMFP: M=6','cMFP: M=20'});  end;
    if strcmp(scenario,'incoherent'),  legend({'MFP','nMFP','cMFP: M=2','cMFP: M=5','cMFP: M=10'});  end;
    if strcmp(scenario,'coherent'),    legend({'MFP','nMFP','cMFP: M=1','cMFP: M=2','cMFP: M=5 '});  end;

    pause;
end;


%NOISE Error
scenariolist = {'single','coherent','incoherent'};
for k0=1:length(scenariolist)
scenario = scenariolist{k0};
load(['data/noise_' scenario '_' phi_type]);
for k=1:5
    scm20(:,k) = scml(k).scm(:,20);
    scmM(:,k) = mean(scml(k).scm>1)';
end;
semilogy(scm20,pp);
axis([0 1 .01 1]);
legend({' 16 dB','12 dB','8 dB','4 dB','0 dB'});
clipboard('copy',['noise_' scenario '_' phi_type]);
pause;

semilogy(scmM);
axis([min(Mlist) max(Mlist) .01 1]);
legend({' 16 dB','12 dB','8 dB','4 dB','0 dB'});
clipboard('copy',['noise_vsM_' scenario '_' phi_type]);
pause;
end;


% LOBE?
scenario = 'single';
phi_list = {'QR-Bernoulli','Gaussian','Bernoulli','QR-Gaussian'};
for k0=1:length(phi_list)
    phi_type = phi_list{k0};
    load(['data/lobe_' scenario '_' phi_type]);
    mcm(:,k0) = median(cmfpr)';
%     mcnm(:,k0) = median(nmfpr)';
end;
clipboard('copy',['lobe_' scenario '_' phi_type]);
plot(Mlist,mfpr(1)*ones(size(Mlist)),'--',Mlist,mcm);
legend({'MFP','cMFP: QR-Bernoulli','cMFP: Gaussian','cMFP: Bernoulli','cMFP: QR-Gaussian'});
pause;


scenario = 'coherent';
for k0=1:length(phi_list)
    phi_type = phi_list{k0};
    load(['data/lobe_' scenario '_' phi_type]);
    mcm(:,k0) = median(cmfpr)';
%     mcnm(:,k0) = median(nmfpr)';
end;
clipboard('copy',['lobe_' scenario '_' phi_type]);
plot(Mlist,mfpr(1)*ones(size(Mlist)),'--',Mlist,mcm);
legend({'MFP','cMFP: QR-Bernoulli','cMFP: Gaussian','cMFP: Bernoulli','cMFP: QR-Gaussian'});
pause;




%MODERR
scenario = 'coherent';
load(['data/moderr_' scenario '_' phi_type]);
for k=1:6
    mfp(k)  = mean(scml(k).snm);
    cmfp(k) = mean(scml(k).scm);
end;
plot(1520:2:1530,mfp,'-+',1520:2:1530,cmfp,'-o');
legend({'MFP','cMFP'});
clipboard('copy',['moderr_' scenario '_' phi_type]);
pause;
clear mfp cmfp;



%PATH
scenario = 'coherent';
load(['data/path_' scenario '_' phi_type]);
plot(targets(1,:),targets(2,:),'--',path_nmfp(1,:),path_nmfp(2,:),'x',path_cmfp(1,:),path_cmfp(2,:),'+');
axis ij;
axis([min(targets(1,:)) max(targets(1,:)) 0 200]);
legend({'Trajectory','cMFP','MFP'});
clipboard('copy',['path1_' scenario '_' phi_type]);
% max(normByCols(path_cmfp-targets))
% max(normByCols(path_mfp-targets))
pause;

load(['data/path2_' scenario '_' phi_type]);
plot(targets(1,:),targets(2,:),'--',path_nmfp(1,:),path_nmfp(2,:),'x',path_cmfp(1,:),path_cmfp(2,:),'+');
axis ij;
axis([min(targets(1,:)) max(targets(1,:)) 0 200]);
legend({'Trajectory','cMFP','MFP'});
clipboard('copy',['path2_' scenario '_' phi_type]);
% max(normByCols(path_cmfp-targets))
% max(normByCols(path_mfp-targets))
