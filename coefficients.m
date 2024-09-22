% coefficient

mis=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\mis_hourly_avg.mat'); %path to mis_hourly_avg
mel=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\mel_hourly_avg.mat'); %path to mis_hourly_avg
adcp=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\adcp_hourly_avg.mat'); %path to adcp_hourly_avg
m=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\m.mat'); %path to m
Hmis=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\Hmis.mat'); %path to Hmis
rho_w=1015;
c_w=4000;

% adcp.speed=adcp.speed./2; to investigate the impact of changing the channel width

%% Using a GammaT/GammaS ratio of 30
% using a ratio enable us to compute Tb and Sb using 3eq formulation
% we start by using an Cd^(1/2)GammaT estimation then computing m from the
% 3eq param, then integrating over the channel and comparing the total
% melting with m obtained from heat budget.

c=[0.0003 0.0005 0.0007 0.0009 0.0011];
for k=1:5
    CdGammaT=c(k);
    CdGammaS=CdGammaT/30;
    Cd=0.0077; % does not matter
    G_T=CdGammaT./(Cd^0.5);
    G_S=CdGammaS./(Cd^0.5);
    
    F_eq_1=zeros(1,length(mel.time));
    F_eq_2=zeros(1,length(mel.time));
    a_mis=nan(42,8401);
    a_mel=nan(42,8401);
    
    for i=1:length(adcp.time)
        % top at mis (x=L)
        U_w=(adcp.speed(2,i));
        T_w=mis.temp(2,i);
        S_w=mis.salinity(2,i);
        [aa1(i),~,~,~] = melt_3eq(Cd,G_T,G_S,U_w,T_w,S_w,8,2);
        % top at mel (x=0)
        U_w=adcp.speed(2,i);
        T_w=mel.temp(2,i);
        S_w=mel.salinity(2,i);
        [aa2(i),T_b(i),S_b(i),a_error(i)] = melt_3eq(Cd,G_T,G_S,U_w,T_w,S_w,8,2);
        %top together
        F_eq_1(i)=(aa1(i)+aa2(i))./3600./24./2.*18;
    
        for j=1:adcp.H(i)
            % sides at mis (x=L)
            U_w=adcp.speed(j,i);
            T_w=mis.temp(j,i);
            S_w=mis.salinity(j,i);
            [a_mis(j,i),~,~,~] = melt_3eq(Cd,G_T,G_S,U_w,T_w,S_w,adcp.depth(j),2);
            % sides at mel (x=0)
            U_w=adcp.speed(j,i);
            T_w=mel.temp(j,i);
            S_w=mel.salinity(j,i);
            [a_mel(j,i),~,~,~] = melt_3eq(Cd,G_T,G_S,U_w,T_w,S_w,adcp.depth(j),2);
            % together
            F_eq_2(i)=F_eq_2(i)+((a_mel(j,i)+a_mis(j,i))./3600./24./2.*0.62);

        end
    end

    m_coef(k,:)=(F_eq_1+F_eq_2).*16000.*3600.*24; %in m^3/day
end

%% frazil
frazil_L=8000; %length over which there is frazil. 8000 is half the channel
C=7e-5; %frazil concentration
r=0.0005; %frazil crystals radius

f_frazil=nan(42,length(adcp.time));
T_b=nan(42,length(adcp.time));
for i=1:length(adcp.time)
    i;
    for j=1:adcp.H(i)
        T_w=mis.temp(j,i);
        S_w=mis.salinity(j,i);
        [f_frazil(j,i),T_b(j,i),~,~] = melt_3eq_frazil(C,r,T_w,S_w,adcp.depth(j),3);
    end
end

f_frazil(f_frazil>0)=0;
frazil=sum(f_frazil.*adcp.w,'omitnan').*0.2; %in m^3 per day
frazil=frazil.*frazil_L;

%%
for i=1:length(adcp.time)
    H(i)=adcp.depth(adcp.H(i)); %depth of outflow layer (for plots)
end

figure

subplot(4,1,2)
avg=48;
temp=movmean(m,avg);
frazil=movmean(frazil,avg);
plot(adcp.time(avg:end-avg),temp(avg:end-avg),'-k','linewidth',1.5); hold on
cmap=cmocean('ice',9);
lw=[1.1 1.5 1.1 1.1 1.1];
for i=2
    plot(mel.time(avg:end-avg),movmean(m_coef(i,avg:end-avg),avg),'color','r','linewidth',1.5);
    plot(mel.time(avg:end-avg),movmean(m_coef(i,avg:end-avg),avg)+frazil(avg:end-avg),'--','color','r','linewidth',lw(i));
end
plot(mel.time([avg end-avg]),[0 0],'--k')
xlim([mel.time(avg) mel.time(end-avg)])
legend('$\dot{m}$ from budget',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0005',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0005+frazil','interpreter','latex');
datetick('x','keeplimits')
set(gca,'linewidth',1.3,'fontsize',10)
xh=get(gca,'position');
xh(3)=0.70;
set(gca,'position',xh);
ylabel('Melt [m\textsuperscript{3} d\textsuperscript{-1}]','Interpreter','latex')
xticklabels('');


subplot(4,1,1)
avg=48;
temp=movmean(m,avg);
frazil=movmean(frazil,avg);
plot(adcp.time(avg:end-avg),temp(avg:end-avg),'-k','linewidth',1.5); hold on
i=1;
plot(mel.time(avg:end-avg),movmean(m_coef(i,avg:end-avg),avg),'color',cmap(i+1,:),'linewidth',1);
i=2;
plot(mel.time(avg:end-avg),movmean(m_coef(i,avg:end-avg),avg),'color','r','linewidth',1.3);
for i=[3 4 5]
    plot(mel.time(avg:end-avg),movmean(m_coef(i,avg:end-avg),avg),'color',cmap(i+1,:),'linewidth',1);
end
plot(adcp.time(avg:end-avg),temp(avg:end-avg),'-k','linewidth',1.5); hold on

xlim([mel.time(avg) mel.time(avg)+52])
datetick('x','keeplimits')
legend('$\dot{m}$ from budget',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0003',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0005',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0007',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0009',...
    '$\dot{m}$ from $C_d^{0.5}\Gamma_{\Theta}$=0.0011','','interpreter','latex','numcolumns',1);
set(gca,'linewidth',1.3,'fontsize',10)
xh=get(gca,'position');
xh(3)=0.45;
set(gca,'position',xh);
ylabel('Melt [m\textsuperscript{3} d\textsuperscript{-1}]','Interpreter','latex')


subplot(4,1,3)
pcolor(mis.time,mis.depth,a_mel.*365); hold on
shading flat
caxis([-100 100])
cmocean('balance')
plot(adcp.time,H,'k','linewidth',0.8); hold on
axis ij
datetick('x','keeplimits')
x=get(gca,'position');
x(3)=0.7;
x(4)=0.18;
set(gca,'position',x);
h=colorbar;
set(h,'linewidth',1.3)
x=get(h,'position');
x(3)=0.02;
x(1)=0.84;
set(h,'position',x)
xticklabels('');
ylabel('Depth [m]','interpreter','latex')
ylabel(h,'$a_B$ [m a\textsuperscript{-1}]','Interpreter','latex','fontsize',10)
set(gca,'linewidth',1.3,'fontsize',10)
set(gca,'Layer','top')


subplot(4,1,4)
pcolor(mis.time,mis.depth,a_mis.*365); hold on
shading flat
caxis([-5 5])
cmocean('balance')
plot(adcp.time,H,'k','linewidth',0.8); hold on
axis ij
datetick('x','keeplimits')
x=get(gca,'position');
x(3)=0.7;
x(4)=0.18;
set(gca,'position',x);
h=colorbar;
set(h,'linewidth',1.3)
x=get(h,'position');
x(3)=0.02;
x(1)=0.84;
set(h,'position',x)
ylabel('Depth [m]','interpreter','latex')
ylabel(h,'$a_B$ [m a\textsuperscript{-1}]','Interpreter','latex','fontsize',10)
set(gca,'linewidth',1.3,'fontsize',10)
set(gca,'Layer','top')
