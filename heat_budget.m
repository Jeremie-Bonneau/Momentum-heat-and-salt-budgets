% heat budget
mis=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\mis_hourly_avg.mat'); %path to mis_hourly_avg
mel=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\mel_hourly_avg.mat'); %path to mis_hourly_avg
adcp=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\adcp_hourly_avg.mat'); %path to adcp_hourly_avg
K=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\K.mat'); %path to K

%%
rho_w=1015;
rho_i=920;
c_w=4000;
L=330000;
c_i=2000;

% Hout - Hin

% %check influence of cooling on sides of MEL
% mel.CTf=gsw_CT_freezing(mel.salinity,mel.depth);
% mel.temp=mel.temp-0.2;
% mel.temp(mel.temp<mel.CTf)=mel.CTf(mel.temp<mel.CTf);

h=mis.temp-mel.temp;
h=h.*adcp.along.*0.2;
h=h.*adcp.w;
h = sum(h).*rho_w.*c_w; 

% Hmix
dT_dz_mel=diff(mel.temp)./0.2; %0.2 is the vertical bin size
dT_dz_mel=movmean(dT_dz_mel,10); %smoothing gradients because otherwise very chunky
dT_dz_mis=diff(mis.temp)./0.2; %0.2 is the vertical bin size
dT_dz_mis=movmean(dT_dz_mis,10); %smoothing gradients because otherwise very chunky
% bottom flux
for i=1:length(adcp.time)
    grad_mel(i)=(mel.temp(adcp.H(i),i)-mel.temp(adcp.H(i)-1,i))./0.2;
    grad_mis(i)=(mis.temp(adcp.H(i),i)-mis.temp(adcp.H(i)-1,i))./0.2;
    grad(i)=(grad_mel(i)+grad_mis(i))./2;
    W=adcp.w(adcp.H(i)); % width of the channel at that depth
    area_grad=grad(i).*W.*16000;
    h_bot(i)=rho_w.*c_w.*K(i).*area_grad;
end

% total
Hmis=(abs(h)-abs(h_bot));
m=(abs(h)-abs(h_bot))./(L+c_i.*(1.5+10)); %kg/s
m=m.*3600.*24./920; %m3 per day

% visual
figure("Position",[50 50 800 800])
subplot(5,1,1)
plot(mel.time,movmean(K,1),'Color',[0.45 0.45 0.45 0.2],'Linewidth',1)
hold on
plot(mel.time,movmean(K,7*24),'color',[0.45 0.45 0.45 1.0],'Linewidth',1.7)
ylabel('$K$ [m\textsuperscript{2} s\textsuperscript{-1}]','Interpreter','latex')
xlim([mel.time(1) mel.time(end)])
datetick('x','keeplimits')
xticklabels('');
x=get(gca,'position');
x(3)=0.8;
x(4)=0.15;
set(gca,'position',x);
set(gca, 'LineWidth',1.3,'FontSize',10)


subplot(5,1,2)
plot(mel.time,-movmean(h_bot./1e9,1),'color',[0.2457    0.2960    0.5885 0.2],'Linewidth',1); hold on
plot(mel.time,-movmean(h,1)./1e9,'color',[0.3822    0.6575    0.7820 0.2],'linewidth',1); hold on

plot(mel.time,-movmean(h./1e9,7*24),'color',[0.3822    0.6575    0.7820 1.0],'linewidth',1.5)
plot(mel.time,-movmean(h_bot./1e9,7*24),'color',[0.2457    0.2960    0.5885 1.0],'Linewidth',1.5)
xticklabels('');
x=get(gca,'position');
x(3)=0.8;
x(4)=0.15;
set(gca,'position',x);
ylabel('$H$ [GW]','Interpreter','latex')
legend('','','$H_{in}-H_{out}$','$H_{mix}$','interpreter','latex')
xlim([mel.time(1) mel.time(end)])
datetick('x','keeplimits')
xticklabels('');
set(gca, 'LineWidth',1.3,'FontSize',10)


subplot(5,1,3)
plot(mel.time,movmean(m,1),'Color',[0.45 0.45 0.45 0.2],'Linewidth',1); hold on
plot(mel.time,movmean(m,7*24),'-','color',[0.45 0.45 0.45 1.0],'Linewidth',1.7);
plot(mel.time([1 end]),[0 0],'--k')
x=get(gca,'position');
x(3)=0.8;
x(4)=0.15;
set(gca,'position',x);
ylabel('Melt [m\textsuperscript{3} d\textsuperscript{-1}]','Interpreter','latex')
xlim([mel.time(1) mel.time(end)])
datetick('x','keeplimits')
set(gca, 'LineWidth',1.3,'FontSize',10)
