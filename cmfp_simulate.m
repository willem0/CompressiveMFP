% clear;
close all;
tic;
clear cmnx;
%grid points: range/depth
%dls:         range/depth
%sig:         range/depth

% Defaults
if ~exist('ji');          ji=1;                             end;
if ~exist('Mlist');       Mlist = [2];                      end;
if ~exist('scenario');    scenario='single';                end;
if ~exist('nlevel');      nlevel = 10^(1/5-1);              end;
if ~exist('numsims');     numsims = 10;                     end;
if ~exist('numSources');  numSources = 100;                 end;
if ~exist('sources');     sources = rand(2,numSources);     end;
if ~exist('erad2');       erad2 = [1; 1];                   end;

numP = 90;                      %Number of Points
cMFPratio = [];
nMFPratio = [];
 MFPratio = [];


if strcmp(scenario,'single'),      rs = 5000+(1:(numP))*9;      end;
if strcmp(scenario,'incoherent'),  rs = 5000+(1:(numP))*9;      end;
if strcmp(scenario,'coherent'),    rs = 5000+(1:(numP))*3;      end;
zs = 10+(1:(numP*2));             %Source Depth (m)
zr=10:5:190;                    %Receiver Depth (m)
dls = [mean(diff(rs)); mean(diff(zs))];

lrs = length(rs);
lzs = length(zs);
lzr = length(zr);
lf = length(freq);

trs = (0*zs'+1)*rs;
tzs = zs'*(0*rs+1);
grid_pts = [trs(:) tzs(:)]';

if strcmp(scenario,'single'),      flist = lf/2;                end;
if strcmp(scenario,'incoherent'),  flist = 1:lf;                end;
if strcmp(scenario,'coherent'),    flist = 1:lf;                end;
shmat = diag(erad2.^(-1));

%constrain to the centered unit square
sources(sources<-1/2) = -1/2;
sources(sources> 1/2) =  1/2;
targets = diag([max(rs)-min(rs); max(zs)-min(zs)])*sources +...
               [mean(rs); mean(zs)]*ones(1,numSources);

scm = zeros(numSources,numsims,length(Mlist));
 sm = zeros(numSources,1);
snm = zeros(numSources,1);

G = greensG_mode(psi,z,N_modes,modes,rho_w,zr,grid_pts);
G = G / max(abs(G(:)));

%G is only computed once, there is opportunity to reload the state
%information at every stage here. It is possible to then run the modeling
%error simulation here.




for kj=1:numSources


%reload for the modeling error sims, unnecessary in general
load states/state_141_160_1520.mat; 

target = targets(:,kj);

a=sum((shmat*(grid_pts-target*ones(1,length(grid_pts)))).^2)<1;
% main lobe for the lobe simulation

Gg = greens_mode(psi,z,N_modes,modes,rho_w,zr,target);
Gg = Gg / max(abs(Gg(:)));
nn = (randn(size(Gg))+j*randn(size(Gg)))/sqrt(2*length(Gg));
Ggn = Gg + nn * nlevel*norm(Gg(:)); %*C*

rGgn = reshape(conj(Ggn),lf,lzr,1)'; %conjugate necessary

num  = 0;
num2 = 0;
den2 = 0;
h = 0;
hn = 0;
% Traditional MFP
for kf=1:length(flist);
    k = flist(kf);
    hamb(:,k) = G(:,:,k)'*rGgn(:,k); %can almost delete this
    
    if strcmp(scenario,'coherent')
        num = num + G(:,:,k)'*rGgn(:,k);
        den2 = den2 + sum(abs(G(:,:,k)).^2)';
    else
        num2 = abs(G(:,:,k)'*rGgn(:,k)).^2;
        h  =  h  + num2;
        hn =  hn + num2 ./ sum(abs(G(:,:,k)).^2)';
    end;
end;
if strcmp(scenario,'coherent')
    h  = abs(num).^2;
    hn = h./den2;
end;
[trash,mnx]  = max(h);
[trash,nmnx] = max(hn);
 MFPratio(kj) = sqrt(max(h )/max(h .*~a'));
nMFPratio(kj) = sqrt(max(hn)/max(hn.*~a'));
h = reshape(h,lzs,lrs);
hn = reshape(hn,lzs,lrs);

% imagesc(reshape(abs(sum(hamb,2)),lzs,lrs));
imagesc(rs,zs,h);
hold on;
tt = linspace(0,2*pi,1000);
plot(target(1)+cos(tt)*erad2(1),target(2)+sin(tt)*erad2(2),'k');
pause(.01);
hold off;

% Compressed MFP
for ii=1:length(Mlist) %*C*
    M = Mlist(ii);  %number of test-measurements *C*
    loclist = zeros(numsims,1);
    for jj=1:numsims
        %beware of conjugate transpose issues
        num = 0;
        den2 = 0;
        ht = 0;
        for kf=1:length(flist);
            k = flist(kf);
            
            %Generate Phik
            if findstr(phi_type,'Bernoulli'),     Phik = sign(randn(lzr));
            elseif findstr(phi_type,'Gaussian'),  Phik = randn(lzr);
            end;
            
            if findstr(phi_type,'QR'),
                [Phik,R]=qr(Phik');
            end;
            Phik=Phik(:,1:M)';
            
            PG = Phik * G(:,:,k);
            
            hwt = PG' * (Phik*rGgn(:,k));
            num = num + hwt;
            if strcmp(scenario,'coherent')
                den2 = den2 + sum(abs(PG).^2);
            else
                ht = ht + abs( hwt ).^2./sum(abs(PG).^2)';
            end;
        end;
        if strcmp(scenario,'coherent')
            ht = abs(num).^2./den2';
        end;
        [trash,cmnx(jj)] = max(ht);
        anum = abs(num);
        cMFPratio(kj,jj,ii) = max(anum)/max(anum.*~a');
    end;
    
    cMFPerr = grid_pts(:,cmnx)-target*ones(1,numsims);
     MFPerr = grid_pts(:, mnx)-target;
    nMFPerr = grid_pts(:,nmnx)-target;
    path_cmfp(:,kj) = grid_pts(:,mode(cmnx));
    path_mfp(:,kj)  = grid_pts(:,mode( mnx));
    path_nmfp(:,kj) = grid_pts(:,mode(nmnx));
    
    scm(kj,:,ii) = sort(sqrt(  sum((shmat*cMFPerr).^2)  ));
     sm(kj)      =             norm(shmat* MFPerr);
    snm(kj)      =             norm(shmat*nMFPerr);
    disp(['With M=' num2str(M) ' tests, ' num2str(toc/60) ' minutes have elapsed.']);
end;

disp(['*** Source Number: ' num2str(kj)]);
end;
toc
ht = reshape(ht,lzs,lrs);

ji

scm = reshape(scm,size(scm,1)*size(scm,2),length(Mlist));
sm = sort(sm);
scm = sort(scm);
snm = sort(snm);

% pp = 1-linspace(1e-5,1-1e-5,length(sm));
% semilogy(sm,pp,'-+',snm,pp,'--',scm,1-linspace(0,1,size(scm,1)));
% axis([0 1 .01 1]);

L = size(scm,1);

scml(ji).sm  = sm;
scml(ji).scm = scm;
scml(ji).snm = snm;

 mfpr = 20*log10( MFPratio );
nmfpr = 20*log10(nMFPratio );
cmfpr = 20*log10( reshape(cMFPratio,L,length(Mlist)) );

save msT Mlist sm scm snm scml L



% if exist('target_path') %for the case of the path sim
%     target = targets(:,kj);
% else  %in general
%     target = [mean(rs)+std(rs)*randn/2;...
%               mean(zs)+std(zs)*randn/2];
% end;

% 'mfpr','cmfpr','nmfpr');
