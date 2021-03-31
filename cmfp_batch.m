clear;
t0=tic;

% phi_type = 'QR-Bernoulli';
% phi_type = 'QR-Gaussian';
% phi_type = 'Bernoulli';
phi_type = 'Gaussian';

%defaults:
nlevel = 10^(1/5-1);
Mlist = [3 6 20];
numSources = 100;
scenariolist = {'single','coherent','incoherent'};
phi_list = {'QR-Bernoulli','Gaussian','Bernoulli','QR-Gaussian'};

if 1
    load states/state_141_160_1520.mat
else;
    main_modes;
end;



%TAILS
for k0=1:length(scenariolist)
    scenario = scenariolist{k0};
    
    %sources uniformly distributed in the unit square
    sources = rand(2,numSources) - 1/2;
    
    if strcmp(scenario,'single'),      Mlist = [3 6 20];       end;
    if strcmp(scenario,'incoherent'),  Mlist = [2 5 10];       end;
    if strcmp(scenario,'coherent'),    Mlist = [1 2 5];        end;

    if strcmp(scenario,'single'),      erad2 = [36; 3];             end;
    if strcmp(scenario,'incoherent'),  erad2 = [36; 3];             end;
    if strcmp(scenario,'coherent'),    erad2 = [12; 3];             end;

    cmfp_simulate;
    save(['data/tail_' scenario '_' phi_type],'Mlist','sm','scm','snm');
    
    pp = 1-linspace(1e-5,1-1e-5,length(sm));
    semilogy(snm,pp,'--',sm,pp,'--',scm,1-linspace(0,1,size(scm,1)));
    axis([0 1 .01 1]);
end;


textme('starting noise');

% if 0  %***************CUT************

% % % % % % % % Time Hog
%NOISE Error
for k0=1:length(scenariolist)
    scenario = scenariolist{k0};
    % scenario = 'single';
    erad2 = [36; 3];
    Mlist = 1:20;
    for ji=1:5
        sources = rand(2,numSources) - 1/2;
        nlevel = 10^(ji/5-1);
        cmfp_simulate;
    end;
    save(['data/noise_' scenario '_' phi_type],'Mlist','scml');
    clear scml;
end;

textme('finished noise');


% LOBE
nlevel = 0;
sources = zeros(2,numSources);
Mlist = 1:37;
erad2 = [180; 16]; %target ellipse (20 freq, single)
scenario = 'single';
phi_list = {'QR-Bernoulli','Gaussian','Bernoulli','QR-Gaussian'};
for k0=1:length(phi_list)
    phi_type = phi_list{k0};
    cmfp_simulate;
    save(['data/lobe_' scenario '_' phi_type],'Mlist','mfpr','cmfpr','nmfpr');
end;


% % % % % % % % Time Hog
erad2 = [72; 16]; %target ellipse (20 freq, coherent)
scenario = 'coherent';
for k0=1:length(phi_list)
    phi_type = phi_list{k0};
    cmfp_simulate;
    save(['data/lobe_' scenario '_' phi_type],'Mlist','mfpr','cmfpr','nmfpr');
end;
% end; %***************CUT************



%CROSS
scenario = 'single';
nlevel = 0;
sources = [0; 0];
numSources = 1;
numsims = 1;
Mlist = [10];  cmfp_simulate;  ht1 = ht;
Mlist = [30];  cmfp_simulate;  ht2 = ht;
Mlist = [10 30];
rad = [3 16 36 180];
save(['data/amb_' scenario '_' phi_type],'hn','ht1','ht2','rs','zs','target','Mlist','rad');
scenario = 'coherent';
Mlist = [2];   cmfp_simulate;  ht1 = ht;
Mlist = [20];  cmfp_simulate;  ht2 = ht;
Mlist = [2 20];
rad = [3 16 12 72];
save(['data/amb_' scenario '_' phi_type],'hn','ht1','ht2','rs','zs','target','Mlist','rad');


%MODERR
scenario = 'coherent';
nlevel = 0;
erad2 = [1; 1];
Mlist = [2];
numSources = 100;
numsims = 1;
sources = zeros(2,numSources);
for ji=1:6
    load(['states/state_141_160_' num2str(1518+2*ji) '.mat']);
    cmfp_simulate;
end;
save(['data/moderr_' scenario '_' phi_type],'Mlist','scml');
clear scml;


%PATH
scenario = 'coherent';
numsims = 1;
Mlist = [2];
numSources = 100;
nx = (1:numSources)/numSources-1/2;
sources = [nx; 3*(nx).^2-1/3];
nlevel = 10^(1/5-1);   cmfp_simulate;
save(['data/path_' scenario '_' phi_type],'Mlist','targets','path_nmfp','path_mfp','path_cmfp');
nlevel = 10^(3/5-1);   cmfp_simulate;
save(['data/path2_' scenario '_' phi_type],'Mlist','targets','path_nmfp','path_mfp','path_cmfp');
clear numsims numSources sources;


disp(['Final elapsed time is ' num2str(toc(t0)/3600) ' hours.']);
datestr(now)
return;

%"final" run:
% x noise plot for incoherent, coherent
% normlized info for lobe, x moderr, x path
% x extra values for Mlist for coherent, incoherent
% x more noise for path
