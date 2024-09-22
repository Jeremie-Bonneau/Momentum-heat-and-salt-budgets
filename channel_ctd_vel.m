%% Jill current and CTD profiles
%ctd profiles 
ds=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Ellesmere_2015\CTD\Milne2015Jul_ctdd.mat');

%% StC (closest to offshore)
% currents
%directly from Jill's thesis p77
z=[2.5 3.5   4.5   5.5   6.5   7.5   8.5   9.5  10.5];
v=[0   0.37  0.49  0.60  0.40  0.2   0.11  0.05 0];
w=[0   44 64 70    78    80    84    87   90];  %adjusted to take into account the surface topography 

c_c.min=2.5; %ice draft at StC
c_c.z=[c_c.min:0.1:18];
c_c.v=interp1(z,v,c_c.z);
c_c.w=interp1(z,w,c_c.z);

% temp and sal at StC
ctdc=ds(10);
r=2;
c_c.CT=interp1(movmean(-ctdc.z(10:400),r),movmean(ctdc.CT(10:400),r),c_c.z);
c_c.SA=interp1(movmean(-ctdc.z(10:400),r),movmean(ctdc.SA(10:400),r),c_c.z);
c_c.CTf=interp1(movmean(-ctdc.z(10:400),r),movmean(ctdc.CTf(10:400),r),c_c.z);
c_c.rho=interp1(movmean(-ctdc.z(10:400),r),movmean(ctdc.rho(10:400),r),c_c.z);


%% StB, mis mooring
% currents at b
%directly from Jill's thesis p78
z=[7 7.5  8.5  9.8  10.5  11.5  12.5  13.5  14.5  15.5  16.5];
v=[0   0.13 0.33 0.46 0.45  0.42  0.35  0.23  0.1   0.05  0];
w=[0   17.1 20.5 21.4 23.43 24.28 24.28 25.33 26.33 27.33 28.33];
c_b.min=7.0; %ice draft
c_b.z=[c_b.min:0.1:18];
c_b.v=interp1(z,v,c_b.z);
c_b.w=interp1(z,w,c_b.z);

% temp and sal
ctdb=ds(5);
r=2;
c_b.CT=interp1(movmean(-ctdb.z,r),movmean(ctdb.CT,r),c_b.z);
c_b.SA=interp1(movmean(-ctdb.z,r),movmean(ctdb.SA,r),c_b.z);
c_b.CTf=interp1(movmean(-ctdb.z,r),movmean(ctdb.CTf,r),c_b.z);
c_b.rho=interp1(movmean(-ctdb.z,r),movmean(ctdb.rho,r),c_b.z);


%% StA (mel)
ctda=ds(21);
c_a.min=2.2; %ice draft
c_a.z=[c_a.min:0.1:18];
r=2;
[ctda.CT,~]=CTD_deloop(ctda.CT,ctda.press);
[ctda.SA,~]=CTD_deloop(ctda.SA,ctda.press);
[ctda.CTf,~]=CTD_deloop(ctda.CTf,ctda.press);
[ctda.rho,~]=CTD_deloop(ctda.rho,ctda.press);
[ctda.z,~]=CTD_deloop(ctda.z,ctda.press);

c_a.CT=interp1(movmean(-ctda.z,r),movmean(ctda.CT,r),c_a.z);
c_a.SA=interp1(movmean(-ctda.z,r),movmean(ctda.SA,r),c_a.z);
c_a.CTf=interp1(movmean(-ctda.z,r),movmean(ctda.CTf,r),c_a.z);
c_a.rho=interp1(movmean(-ctda.z,r),movmean(ctda.rho,r),c_a.z);


%% offshore
ctdo=ds(12);
c_o.min=2.0; %ice draft
c_o.z=[c_o.min:0.1:18];
r=2;

[ctdo.CT,~]=CTD_deloop(ctdo.CT,ctdo.press);
[ctdo.SA,~]=CTD_deloop(ctdo.SA,ctdo.press);
[ctdo.CTf,~]=CTD_deloop(ctdo.CTf,ctdo.press);
[ctdo.rho,~]=CTD_deloop(ctdo.rho,ctdo.press);
[ctdo.z,~]=CTD_deloop(ctdo.z,ctdo.press);

c_o.CT=interp1(movmean(-ctdo.z,r),movmean(ctdo.CT,r),c_o.z);
c_o.SA=interp1(movmean(-ctdo.z,r),movmean(ctdo.SA,r),c_o.z);
c_o.CTf=interp1(movmean(-ctdo.z,r),movmean(ctdo.CTf,r),c_o.z);
c_o.rho=interp1(movmean(-ctdo.z,r),movmean(ctdo.rho,r),c_o.z);




%% Figure
cmap=parula(7);
figure
subplot(3,4,1)
plot(c_a.CT,c_a.z,'color','c','linewidth',1.5); hold on
plot(c_a.CTf,c_a.z,'--','color','c','linewidth',1); hold on
plot([-2 5],[c_a.min c_a.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([-1.8 3.5])
ylabel('Depth [m]','interpreter','latex')
xlabel('$\Theta\; [^\circ \textrm{C}] $','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,2)
plot(c_b.CT,c_b.z,'color',cmap(3,:),'linewidth',1.5); hold on
plot(c_b.CTf,c_b.z,'--','color',cmap(3,:),'linewidth',1); hold on
plot([-2 5],[c_b.min c_b.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([-1.8 3.5])
yticklabels('')
xlabel('$\Theta\; [^\circ \textrm{C}] $','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,3)
plot(c_c.CT,c_c.z,'color',cmap(1,:),'linewidth',1.5); hold on
plot(c_c.CTf,c_c.z,'--','color',cmap(1,:),'linewidth',1); hold on
plot([-2 5],[c_c.min c_c.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([-1.8 3.5])
yticklabels('')
xlabel('$\Theta\; [^\circ \textrm{C}] $','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,4)
plot(ctdo.CT,-ctdo.z,'-k','linewidth',1.5); hold on
plot(ctdo.CTf,-ctdo.z,'--k','linewidth',1); hold on
plot([-2 50],[c_o.min c_o.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([-1.8 3.5])
yticklabels('')
xlabel('$\Theta\; [^\circ \textrm{C}] $','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij


subplot(3,4,5)
plot(c_a.SA,c_a.z,'color','c','linewidth',1.5); hold on
plot([-2 50],[c_a.min c_a.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([0 35])
ylabel('Depth [m]','interpreter','latex')
xlabel('$S_A$ [g kg\textsuperscript{-1}]','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,6)
plot(c_b.SA,c_b.z,'color',cmap(3,:),'linewidth',1.5); hold on
plot([-2 50],[c_b.min c_b.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([0 35])
yticklabels('')
xlabel('$S_A$ [g kg\textsuperscript{-1}]','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,7)
plot(c_c.SA,c_c.z,'color',cmap(1,:),'linewidth',1.5); hold on
plot([-2 50],[c_c.min c_c.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([0 35])
yticklabels('')
xlabel('$S_A$ [g kg\textsuperscript{-1}]','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,8)
plot(ctdo.SA,-ctdo.z,'-k','linewidth',1.5); hold on
plot([-2 50],[c_o.min c_o.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([0 35])
yticklabels('')
xlabel('$S_A$ [g kg\textsuperscript{-1}]','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij



subplot(3,4,10)
plot(c_b.v,c_b.z,'color',cmap(3,:),'linewidth',1.5); hold on
plot([-2 50],[c_b.min c_b.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([0 0.6])
ylabel('Depth [m]','interpreter','latex')
xlabel('$U$ [m s\textsuperscript{-1}]','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij

subplot(3,4,11)
plot(c_c.v,c_c.z,'color',cmap(1,:),'linewidth',1.5); hold on
plot([-2 50],[c_c.min c_c.min],'color',[0.5 0.5 0.5],'linewidth',2)
ylim([0 18])
xlim([0 0.6])
yticklabels('')
xlabel('$U$ [m s\textsuperscript{-1}]','Interpreter','latex')
set(gca,'linewidth',1.2,'fontsize',10)
x=get(gca,'position');
x(3)=0.18;
set(gca,'position',x);
axis ij



%% computing KE and PE
% KE
KEb=sum((c_b.v.^2./2).*c_b.w,'omitnan').*0.1; % at StB
KEc=sum((c_c.v.^2./2).*c_c.w,'omitnan').*0.1; % at StC


%PE
%fill top with nan
c_o.rho=[nan(1,20) c_o.rho];
c_a.rho=[nan(1,22) c_a.rho];
c_b.rho=[nan(1,70) c_b.rho];
c_c.rho=[nan(1,25) c_c.rho];


c_b.w=[nan(1,70) c_b.w];
c_c.w=[nan(1,25) c_c.w];
c_a.w=c_b.w;
c_a.w(1:72)=25;

%compute
PEb=(c_o.rho-c_b.rho).*c_b.w.*[0:0.1:18];
PEb=sum(PEb,'omitnan').*0.1.*9.8./1020;
%
PEc=(c_o.rho-c_c.rho).*c_c.w.*[0:0.1:18];
PEc=sum(PEc,'omitnan').*0.1.*9.8./1020;
% 

% PEa=(c_o.rho-c_a.rho).*c_a.w.*[0:0.1:18]; not reliable
% PEa=sum(PEa,'omitnan').*0.1.*9.8./1020; not reliable



%% bounds on PE
mis=importdata('mis_hourly_avg.mat');
adcp=importdata('adcp_hourly_avg.mat');
mel=importdata('mel_hourly_avg.mat');


%PEa
c_b.z=[nan(1,70) c_b.z];
rho=gsw_rho(mis.salinity,mis.temp,mis.depth);
for i=1:8401
    temp=interp1(mis.depth,rho(:,i),c_b.z);
    temp=temp.*c_b.w.*c_b.z;
    temp=sum(temp(66:165).*0.1.*9.8./1020,'omitnan');
    PE_ds(i)=temp;
end

figure; plot(diff(PE_ds)) %of diff(PE total) because this way we don't nned to compare to a (~constant) reference (offshore)


