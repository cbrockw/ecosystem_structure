% Data Preparation 
close all; clear all;

%mdlfile='/Data/cobialab/conserveandclimate/model_data/tte_biomass_genpred_ap_20171017.nc';
mdlfile='/Data/cobialab/conserveandclimate/model_data/tte_biomass_corr2_20170817.nc';

total_mass=ncread(mdlfile,'total_mass');
mass_dens=ncread(mdlfile,'mass_dens');
orgmass=ncread(mdlfile,'org_mass');
tte=ncread(mdlfile,'tte');
tweb=ncread(mdlfile,'flow_matrix');

ncat=length(orgmass);
numscen=4;
niter = length(total_mass);
niter2 = floor(niter/numscen);

lmassdens=log10(mass_dens);
lmassdens(isinf(lmassdens))=NaN;

%initialize computed variables
k=zeros(niter2,numscen);ttem=zeros(niter2,numscen);
ld=zeros(niter2,numscen);conn=zeros(niter2,numscen);
ppmr=zeros(niter2,numscen);

MD=zeros(niter2,numscen);
TD=zeros(ncat,niter2,numscen);

for jj=1:numscen
    MD(:,jj)=nansum(mass_dens(:,((jj-1)*niter2+1):(jj*niter2)));
    TD(:,:,jj)=mass_dens(:,(jj-1)*niter2+1:(jj*niter2));
end

for ii=1:niter2
    for jj=1:numscen
    %combined simulations
        temp=ii+(jj-1)*niter2;
        S(ii,jj)=length(orgmass(mass_dens(:,temp)>0));
        I=find(isfinite(lmassdens(:,temp)));
        [P,S1]=polyfit(log10(orgmass(I)),lmassdens(I,temp),1);
        k(ii,jj)=P(1);
        ttem(ii,jj)=nanmean(tte(I(2:end),temp));
        L=length(find(squeeze(tweb(:,:,temp)>0)))-1;
        if L>0
            conn(ii,jj)=2*L/(S(ii,jj)*(S(ii,jj)-1));
            ld(ii,jj)=conn(ii,jj).*S(ii,jj);
        else
            conn(ii,jj)=0;
            ld(ii,jj)=0;
        end
    
        TL=zeros(size(I));
        TL(1)=1;
        for kk=2:length(I)
            TL(kk)=sum((tweb(I(kk),I,temp)./sum(tweb(I(kk),I,temp)))'.*TL)+1;
        end
        trophiclevel(ii,jj)=nanmax(TL);
        clear TL
    
        PPMR=zeros(size(I));
        PPMR(1)=0;
        for kk=2:length(I)
            PPMR(kk)=orgmass(kk)./sum((tweb(I(kk),I,temp)./sum(tweb(I(kk),I,temp)))'.*orgmass(I));
        end
        ppmr(ii,jj)=geomean(PPMR(2:end));
        
        clear PPMR
    end
    %progress reporting
    if rem(ii,(niter2/20))==0
        disp([num2str(ii/niter2*100),'% done...']);
    end
end

k_t=0.25+log(ttem)./(log(ppmr));

%% Woodson et al 2018 Figure 2 (Size Spectral Plots for Linear, LGP, GSC, and LGP + GSC)

figure('position',[400 400 800 800]);
cc=get(gca,'colororder');
orgfit=[1e-10 1e1];

I=find(tte(2,1:1000)>=.0999 & tte(2,1:1000)<=.1001);

subplot(221)
scatter(orgmass,mass_dens(:,I(2)),80,cc(2,:),'filled'); hold on;
set(gca,'xscale','log','yscale','log','fontsize',16,'fontweight','bold');
box on; grid on; axis([1e-10 1e5 1e-4 1e0]);
[P,S]=polyfit(log10(orgmass),log10(mass_dens(:,I(2))),1);
massfit1=10.^(P(2)).*orgmass.^P(1);
plot(orgmass,massfit1,'color',cc(2,:),'linewidth',1);
ylabel('Biomass $B_j$ [kg m$^{-2}$]','Interpreter','latex','fontsize',24);
%text(5e-10,2e-1,['$B_{tot}$ = ',num2str(sum(mass_dens(:,43)),'%0.2f'),' kg m$^{-2}$'],'Interpreter','latex','fontsize',16);
text(5e-10,5e-1,'a','fontsize',20);

I=find(tte(2,1001:2000)>=.0999 & tte(2,1001:2000)<=.1001)+1000;

subplot(222)
scatter(orgmass,mass_dens(:,I(2)),80,cc(1,:),'filled'); hold on;
set(gca,'xscale','log','yscale','log','fontsize',16,'fontweight','bold');
box on; grid on; axis([1e-10 1e5 1e-4 1e0]);
[P,S]=polyfit(log10(orgmass),log10(mass_dens(:,I(2))),1);
massfit=10.^(P(2)).*orgmass.^P(1);
plot(orgmass,massfit,'color',cc(1,:),'linewidth',1);

[P,S]=polyfit(log10(orgmass(1:end-2)),log10(mass_dens(1:end-2,I(2))),1);
massfit1=10.^(P(2)).*orgfit.^P(1);
plot(orgfit,massfit1,'color',cc(1,:),'linestyle','--');
%text(5e-10,2e-1,['$B_{tot}$ = ',num2str(sum(mass_dens(:,1102)),'%0.2f'),' kg m$^{-2}$'],'Interpreter','latex','fontsize',16);
text(5e-10,5e-1,'b','fontsize',20);

I=find(tte(2,2001:3000)>=.0999 & tte(2,2001:3000)<=.1001)+2000;

subplot(223)     
scatter(orgmass,mass_dens(:,I(2)),80,cc(3,:),'filled'); hold on;
set(gca,'xscale','log','yscale','log','fontsize',16,'fontweight','bold');
box on; grid on; axis([1e-10 1e5 1e-4 1e0]);
[P,S]=polyfit(log10(orgmass),log10(mass_dens(:,I(2))),1);
massfit=10.^(P(2)).*orgmass.^P(1);
plot(orgmass,massfit,'color',cc(3,:),'linewidth',1);

[P,S]=polyfit(log10(orgmass(1:end-1)),log10(mass_dens(1:end-1,I(2))),1);
massfit1=10.^(P(2)).*orgfit.^P(1);
plot(orgfit,massfit1,'color',cc(3,:),'linestyle','--');
xlabel('Individual Mass $M_j$ [kg]','Interpreter','latex','fontsize',24);
ylabel('Biomass $B_j$ [kg m$^{-2}]$','Interpreter','latex','fontsize',24);
%text(5e-10,2e-1,['$B_{tot}$ = ',num2str(sum(mass_dens(:,2384)),'%0.2f'),' kg m$^{-2}$'],'Interpreter','latex','fontsize',16);
text(5e-10,5e-1,'c','fontsize',20);

I=find(tte(2,3001:4000)>=.0999 & tte(2,3001:4000)<=.1001)+3000;

subplot(224)
scatter(orgmass,mass_dens(:,I(2)),80,cc(5,:),'filled'); hold on;
set(gca,'xscale','log','yscale','log','fontsize',16,'fontweight','bold');
box on; grid on; axis([1e-10 1e5 1e-4 1e0]);
[P,S]=polyfit(log10(orgmass),log10(mass_dens(:,I(2))),1);
massfit=10.^(P(2)).*orgmass.^P(1);
plot(orgmass,massfit,'color',cc(5,:),'linewidth',1);

[P,S]=polyfit(log10(orgmass(1:end-2)),log10(mass_dens(1:end-2,I(2))),1);
massfit1=10.^(P(2)).*orgfit.^P(1);
plot(orgfit,massfit1,'color',cc(5,:),'linestyle','--');
xlabel('Individual Mass $M_j$  [kg]','Interpreter','latex','fontsize',24);
%text(5e-10,2e-1,['$B_{tot}$ = ',num2str(sum(mass_dens(:,3186)),'%0.2f'),' kg m$^{-2}$'],'Interpreter','latex','fontsize',16);
text(5e-10,5e-1,'d','fontsize',20);



%% Woodson et al Figure 3 (Biomass and k versus theoretical for size structured ecosystems)

for ii=1:numscen
    [N,edges,idx]=histcounts(ttem(:,ii),40);
    conn_mean(:,ii) = accumarray(idx(:),conn(:,ii),[],@mean);
    mass_mean(:,ii) = accumarray(idx(:),MD(:,ii),[],@mean);
    mass_std(:,ii) = accumarray(idx(:),MD(:,ii),[],@std);
    k_mean(:,ii) = accumarray(idx(:),k(:,ii),[],@mean);
    k_std(:,ii) = accumarray(idx(:),k(:,ii),[],@std);
    te_mean(:,ii) = accumarray(idx(:),ttem(:,ii),[],@mean);
    tl_mean(:,ii) = accumarray(idx(:),trophiclevel(:,ii),[],@mean);
    tl_max(:,ii) = accumarray(idx(:),trophiclevel(:,ii),[],@max);
    sp_mean(:,ii) = accumarray(idx(:),S(:,ii),[],@mean);
    ppmr_mean(:,ii) = accumarray(idx(:),ppmr(:,ii),[],@mean);
    k_theor(:,ii)=accumarray(idx(:),k_t(:,ii),[],@mean);
end

for ii=1:numscen
    te_diff(:,ii)=4.*(abs((te_mean(:,ii)-nanmean(te_mean(:,ii)))./nanstd(te_mean(:,ii)))).^(-1.25);
end

fa=0.2;

figure('position',[400 400 1200 600]);
cc=get(gca,'colororder');
subplot(121);

%scatter(te_mean(2:end,1),mass_mean(2:end,1),100,cc(2,:),'filled'); hold on;
%scatter(te_mean(2:end,2),mass_mean(2:end,2),100,cc(1,:),'filled');
patch([te_mean(2:end,1); flipud(te_mean(2:end,1)); te_mean(2,1)],[mass_mean(2:end,1)-mass_std(2:end,1); flipud(mass_std(2:end,1)+mass_mean(2:end,1)); -mass_std(2,1)+mass_mean(2,1)],cc(2,:),'facealpha',fa,'edgecolor','none'); hold on;
patch([te_mean(2:end,2); flipud(te_mean(2:end,2)); te_mean(2,2)],[mass_mean(2:end,2)-mass_std(2:end,2); flipud(mass_std(2:end,2)+mass_mean(2:end,2)); -mass_std(2,2)+mass_mean(2,2)],cc(1,:),'facealpha',fa,'edgecolor','none');
patch([te_mean(2:end,3); flipud(te_mean(2:end,3)); te_mean(2,3)],[mass_mean(2:end,3)-mass_std(2:end,3); flipud(mass_std(2:end,3)+mass_mean(2:end,3)); -mass_std(2,3)+mass_mean(2,3)],cc(3,:),'facealpha',fa,'edgecolor','none');
patch([te_mean(2:end,4); flipud(te_mean(2:end,4)); te_mean(2,4)],[mass_mean(2:end,4)-mass_std(2:end,4); flipud(mass_std(2:end,4)+mass_mean(2:end,4)); -mass_std(2,4)+mass_mean(2,4)],cc(5,:),'facealpha',fa,'edgecolor','none');
%patch([te_mean(2:end,5); flipud(te_mean(2:end,5)); te_mean(2,5)],[mass_mean(2:end,5)-mass_std(2:end,5); flipud(mass_std(2:end,4)+mass_mean(2:end,5)); -mass_std(2,5)+mass_mean(2,5)],cc(4,:),'facealpha',fa,'edgecolor','none');
patch([.05 .15 .15 .05 .05],[.01 .01 3 3 .01],[.6 .6 .6],'facealpha',0.1,'edgecolor','none'); hold on;
patch([.085 .115 .115 .085 .085],[.01 .01 3 3 .01],[.6 .6 .6],'facealpha',0.2,'edgecolor','none'); hold on;

plot(te_mean(2:end,1),mass_mean(2:end,1),'color',cc(2,:),'linewidth',4);
plot(te_mean(2:end,2),mass_mean(2:end,2),'color',cc(1,:),'linewidth',4);
plot(te_mean(2:end,3),mass_mean(2:end,3),'color',cc(3,:),'linewidth',4);
plot(te_mean(2:end,4),mass_mean(2:end,4),'color',cc(5,:),'linewidth',4);
%plot(te_mean(2:end,5),mass_mean(2:end,5),'color',cc(4,:),'linewidth',4);
axis square; box on; axis([0 0.2 0 3]);
set(gca,'fontsize',16,'fontweight','bold','yscale','log');
xlabel('$cTE$','Interpreter','latex','fontsize',20,'fontweight','bold');
ylabel('Biomass [kgC m$^{-2}$]','Interpreter','latex','fontsize',24,'fontweight','bold');grid on;
text(.01,2.2,'a','fontsize',24,'fontweight','bold');
leg=legend('base','LGP','GSC','LGP+GSC','loop');
set(leg,'Position',[0.44 0.042 0.141 0.0675]);

SC=2;

subplot(122);
%patch([-.2 -.04 -.04 -.2 -.2],[-.75 -.75 .25 .25 -.75],[.6 .6 .6],'facealpha',0.1,'edgecolor','none'); hold on;
%patch([.07 .16 .16 .07 .07],[-.75 -.75 .25 .25 -.75],[.6 .6 .6],'facealpha',0.2,'edgecolor','none'); hold on;
%patch([-.12 -.08 -.08 -.12 -.12],[-.75 -.75 .25 .25 -.75],[.6 .6 .6],'facealpha',0.2,'edgecolor','none'); hold on;
plot([-.75 .35],[-.75 .35],'k','linewidth',2); hold on;
plot([-.75 .35],[0 0],'k','linewidth',2);
h(1)=scatter(k_theor(:,1),k_mean(:,1),SC.*te_diff(:,1),cc(2,:),'filled');
h(2)=scatter(k_theor(:,2),k_mean(:,2),SC.*te_diff(:,2),cc(1,:),'filled');
h(3)=scatter(k_theor(:,3),k_mean(:,3),SC.*te_diff(:,3),cc(3,:),'filled');
h(4)=scatter(k_theor(:,4),k_mean(:,4),SC.*te_diff(:,4),cc(5,:),'filled');
%h(4)=scatter(k_theor(:,5),k_mean(:,5),SC.*te_diff(:,5),cc(4,:),'filled');

% h(1)=scatter(k_theor(:,1),k_mean(:,1),100,cc(2,:),'filled'); 
% h(2)=scatter(k_theor(:,2),k_mean(:,2),100,cc(1,:),'filled'); 
% h(3)=scatter(k_theor(:,3),k_mean(:,3),100,cc(3,:),'filled'); 
% h(4)=scatter(k_theor(:,4),k_mean(:,4),100,cc(5,:),'filled'); 
set(gca,'fontsize',16,'fontweight','bold');
xlabel('$k_{t}=0.25+log(TE)/log(PPMR)$','Interpreter','latex','fontsize',24,'fontweight','bold'); ylabel('$k_{f}$','Interpreter','latex','fontsize',24,'fontweight','bold');
axis square; grid on; box on; axis([-.25 .1 -.25 .1]);
text(-.235,.08,'b','fontsize',24,'fontweight','bold');

%export_fig NatComm_FIG3.png -png -r300 -transparent


%% Woodson et al Figure 4 (Bar Plot with biomass and k for Linear, LGP, GSC, and Combined

figure('position',[400 400 1200 600]);
cc=get(gca,'colororder');

subplot(121)
for ii=1:4
    [lo(ii), up(ii)]=confint(MD(:,ii),'mean');
    lo(ii)=nanmean(MD(:,ii))-lo(ii); up(ii)=up(ii)-nanmean(MD(:,ii));
    mbar(ii)=bar(ii,nanmean(MD(:,ii)),'edgecolor','none'); hold on;
    mbar(ii).FaceColor=cc(ii,:);
end
mbar(1).FaceColor=cc(2,:);
mbar(2).FaceColor=cc(1,:);
mbar(4).FaceColor=cc(5,:);

errorbar([1:4],nanmean(MD),lo,up,'k','linewidth',2,'linestyle','none');
set(gca,'fontsize',16,'fontweight','bold');
set(gca,'XTick',[1:4],'XTickLabel',{'BASE','LGP','GSC','LGP+GSC'});
ylabel('Biomass [kgC m$^{-2}$]','Interpreter','latex','fontsize',24,'fontweight','bold'); axis square; box on;

subplot(122)
for ii=1:4
    [lo(ii), up(ii)]=confint(k(:,ii),'mean');
    lo(ii)=nanmean(k(:,ii))-lo(ii); up(ii)=up(ii)-nanmean(k(:,ii));
    kbar(ii)=bar(ii,nanmean(k(:,ii)),'edgecolor','none'); hold on;
    kbar(ii).FaceColor=cc(ii,:);
end
kbar(1).FaceColor=cc(2,:);
kbar(2).FaceColor=cc(1,:);
kbar(4).FaceColor=cc(5,:);

errorbar([1:4],nanmean(k),lo,up,'k','linewidth',2,'linestyle','none');
set(gca,'fontsize',16,'fontweight','bold');
set(gca,'fontsize',16,'fontweight','bold');
set(gca,'XTick',[1:4],'XTickLabel',{'BASE','LGP','GSC','LGP+GSC'});
ylabel('$k$','Interpreter','latex','fontsize',24,'fontweight','bold'); axis square; box on;
axis([0.5 4.5 -.125 .05]);


%% NWoodson et al Figure 5 (Bar Plot Comparison with Palmyra)

scale=50/80;

xx=mass_dens(2:5,1101); xx=xx./xx(1);
yy=scale.*[75; 119-75; 140-119; 244-140]; yy=yy./yy(1); %Bradley et al 2017
zz=scale.*[75; 119-75; 140-119; 372-140]; zz=zz./zz(1); %Sandin et al 

figure()
cc=get(gca,'colororder');
patch([-xx(1)/2 -xx(1)/2 xx(1)/2 xx(1)/2 -xx(1)/2],[0 1 1 0 0],cc(5,:),'edgecolor','none');
patch([-xx(2)/2 -xx(2)/2 xx(2)/2 xx(2)/2 -xx(2)/2],[1 2 2 1 1],cc(2,:),'edgecolor','none');
patch([-xx(3)/2 -xx(3)/2 xx(3)/2 xx(3)/2 -xx(3)/2],[2 3 3 2 2],cc(3,:),'edgecolor','none');
patch([-xx(4)/2 -xx(4)/2 xx(4)/2 xx(4)/2 -xx(4)/2],[3 4 4 3 3],cc(1,:),'facealpha',0.3,'edgecolor','none');

offset=2;
patch([offset-yy(1)/2 offset-yy(1)/2 offset+yy(1)/2 offset+yy(1)/2 offset-yy(1)/2],[0 1 1 0 0],cc(5,:),'edgecolor','none');
patch([offset-yy(2)/2 offset-yy(2)/2 offset+yy(2)/2 offset+yy(2)/2 offset-yy(2)/2],[1 2 2 1 1],cc(2,:),'edgecolor','none');
patch([offset-yy(3)/2 offset-yy(3)/2 offset+yy(3)/2 offset+yy(3)/2 offset-yy(3)/2],[2 3 3 2 2],cc(3,:),'edgecolor','none');
patch([offset-yy(4)/2 offset-yy(4)/2 offset+yy(4)/2 offset+yy(4)/2 offset-yy(4)/2],[3 4 4 3 3],cc(1,:),'facealpha',0.3,'edgecolor','none');

box on;
set(gca,'YTick',[-10 10]);
text([1 1 1 1],[0.5 1.5 2.5 3.5],{'Herbivores','Planktivores','Piscivores','LGPs (Sharks)'},'fontsize',16,'fontweight','bold','horizontalalignment','center');
set(gca,'XTick',[0 2],'XTickLabel',{'Model','Palmyra^{13}'},'fontsize',24,'fontweight','bold')
set(gca,'fontsize',16,'fontweight','bold');

%% Woodson et al 2018 Supplementary Figure S3 (Biomass vs Ind Mass LGP + GSC with AP)
figure('position',[400 400 400 400]);
cc=get(gca,'colororder');
orgfit=[1e-10 1e1];

I=find(tte(2,3001:4000)>=.0999 & tte(2,3001:4000)<=.1001)+3000;

scatter(orgmass,mass_dens(:,I(2)),80,cc(6,:),'filled'); hold on;
set(gca,'xscale','log','yscale','log','fontsize',16,'fontweight','bold');
box on; grid on; axis([1e-10 1e5 1e-4 1e0]);
[P,S]=polyfit(log10(orgmass),log10(mass_dens(:,I(2))),1);
massfit=10.^(P(2)).*orgmass.^P(1);
plot(orgmass,massfit,'color',cc(6,:),'linewidth',1);
axis square;

[P,S]=polyfit(log10(orgmass(1:end-3)),log10(mass_dens(1:end-3,I(2))),1);
massfit1=10.^(P(2)).*orgfit.^P(1);
plot(orgfit,massfit1,'color',cc(6,:),'linestyle','--');
xlabel('Individual Mass $M_j$  [kg]','Interpreter','latex','fontsize',24);
%text(5e-10,2e-1,['$B_{tot}$ = ',num2str(sum(mass_dens(:,3186)),'%0.2f'),' kg m$^{-2}$'],'Interpreter','latex','fontsize',16);
%text(5e-10,5e-1,'d','fontsize',20);
ylabel('Biomass $B_j$ [kg m$^{-2}]$','Interpreter','latex','fontsize',24);
%% Woodson et al 2018 Supplementary Figure S5 (Convergence Plot) 
% Need to use different Data File
clear all; close all;

mdlfile='/Data/cobialab/conserveandclimate/model_data/tte_biomass_corr2_20161122.nc';
%mdlfile='/Data/cobialab/conserveandclimate/model_data/tte_biomass_genpred_20170817.nc';

total_mass=ncread(mdlfile,'total_mass');
mass_dens=ncread(mdlfile,'mass_dens');
orgmass=ncread(mdlfile,'org_mass');
tte=ncread(mdlfile,'tte');
tweb=ncread(mdlfile,'flow_matrix');

ncat=length(orgmass);
numscen=4;
niter = length(total_mass);
niter2 = floor(niter/numscen);

lmassdens=log10(mass_dens);
lmassdens(isinf(lmassdens))=NaN;

%initialize computed variables
k=zeros(niter2,numscen);ttem=zeros(niter2,numscen);
ld=zeros(niter2,numscen);conn=zeros(niter2,numscen);
ppmr=zeros(niter2,numscen);

MD=zeros(niter2,numscen);
TD=zeros(ncat,niter2,numscen);

for jj=1:numscen
    MD(:,jj)=nansum(mass_dens(:,((jj-1)*niter2+1):(jj*niter2)));
    TD(:,:,jj)=mass_dens(:,(jj-1)*niter2+1:(jj*niter2));
end

for ii=1:niter2
    for jj=1:numscen
    %combined simulations
        temp=ii+(jj-1)*niter2;
        S(ii,jj)=length(orgmass(mass_dens(:,temp)>0));
        I=find(isfinite(lmassdens(:,temp)));
        [P,S1]=polyfit(log10(orgmass(I)),lmassdens(I,temp),1);
        k(ii,jj)=P(1);

    %progress reporting
        if rem(ii,(niter2/20))==0
            disp([num2str(ii/niter2*100),'% done...']);
        end
    end
end

temp=1;
clear meantest stdtest;
ndata=[2 3 4 5 6 7 8 9 10 20 40 60 80 100 200 400 600 800 1000 2000 4000 6000 8000 10000];

for jj=1:length(ndata)
    for ii=1:10     
        idx=randi(25000,ndata(jj),1);
        meantest(ii,temp)=nanmean(k(idx,4));
        stdtest(ii,temp)=nanstd(k(idx,4));
    end
    temp=temp+1;
end

figure('position',[400 400 1200 600]);
cc=get(gca,'colororder');

patch([ndata fliplr(ndata) ndata(1)],[min(meantest) fliplr(max(meantest)) min(meantest(:,1))],cc(1,:),'facealpha',0.2,'edgecolor','none'); hold on;
patch([ndata fliplr(ndata) ndata(1)],[min(stdtest)./sqrt(ndata) fliplr(max(stdtest)./sqrt(ndata)) min(stdtest(:,1))./sqrt(ndata(1))],cc(2,:),'facealpha',0.2,'edgecolor','none'); hold on;

plot(ndata,nanmean(meantest),'linewidth',2); hold on;
plot(ndata,nanmean(stdtest)./sqrt(ndata),'linewidth',2);
plot([2000 2000],[-.4 .2],'k','linewidth',2);
%set(gca,'yscale','log','fontsize',16,'fontweight','bold');
set(gca,'xscale','log','fontsize',24,'fontweight','bold');
grid on; box on; axis([2 1e4 -.2 .2]);
xlabel('$n$','Interpreter','latex','fontsize',32); ylabel('$\mu, SE$','Interpreter','latex','fontsize',32);
lg1=legend('$\mu$','$SE$');
set(lg1,'Interpreter','latex');






