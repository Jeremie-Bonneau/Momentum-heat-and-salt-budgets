% Salt budget 
mis=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\mis_hourly_avg.mat'); %path to mis_hourly_avg
mel=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\mel_hourly_avg.mat'); %path to mis_hourly_avg
adcp=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\adcp_hourly_avg.mat'); %path to adcp_hourly_avg
m=importdata('C:\Users\Jérémie\OneDrive\Documents\UBC\Directed Study\m.mat'); %path to m, to check influence of melting on budget


%% Sout - Sin
s=mis.salinity-mel.salinity; %g/kg
s=s.*1020; %g/m^3    1020 is density
s=s.*adcp.along;  % g/s/m^2
s=s.*adcp.w.*0.2; % g/s
s=sum(s);

%% K
grad_mel=nan(1,length(mis.time)); %salinity gradient at the bottom of outflow layer at mis mooring
grad_mis=nan(1,length(mis.time)); %salinity gradient at the bottom of outflow layer at mel mooring
grad=nan(1,length(mis.time)); % average salinity gradient at the bottom of outflow layer

for i=1:length(mis.time)
    grad_mel(i)=(mel.salinity(adcp.H(i),i)-mel.salinity(adcp.H(i)-1,i))./0.2; %g/kg/m
    grad_mis(i)=(mis.salinity(adcp.H(i),i)-mis.salinity(adcp.H(i)-1,i))./0.2; %g/kg/m
    grad(i)=(grad_mel(i)+grad_mis(i))/2; %average (assuming linear)
    W=adcp.w(adcp.H(i)); % width of the channel at that depth 
    A(i)=W*16000; %area over which salt is entrained inside outflowing layer 
    K(i)=s(i)./(grad(i)*A(i)); %kg/m/s
    K(i)=K(i)./1020; %m^2/s
end
% 16000 is the length of the channel

%% m
%estimating s flux from submarine melting
m_f=m./3600./24; %m^3/s
for i=1:length(adcp.time)
    s_m(i)=mean(mis.salinity(1:adcp.H(i),i)+mel.salinity(1:adcp.H(i),i))./2; %mean outflow salinity
end

m_f=m_f.*s_m; %salt flux in g m^3/kg/s
m_f=m_f.*1020; %in g/s


m_surf=0.4*2*16*1000*1000/75/24/3600; %average melt
m_surf=s_m.*m_surf;
m_surf=m_surf*1020;
m_surf(adcp.time>datenum(2018,09,01))=0;


%% plot budget
cmap=cmocean('ice',5);
figure
plot(downsample(adcp.time,48),downsample(movmean(s,48),48),'color',cmap(2,:),'linewidth',1.2); hold on
plot(adcp.time,m_surf,'color',cmap(3,:),'linewidth',1.2); hold on
plot(downsample(adcp.time,48),downsample(movmean(m_f,48),48),'color',cmap(4,:),'linewidth',1.2)
set(gca,'linewidth',1.2,'fontsize',12)
xlim([mel.time(48) mel.time(end-48)])
datetick('x','keeplimits')
legend('$S_{in}-S_{out}$','$S_{\dot{m}}$','$S_{surf}$','interpreter','latex')
ylabel('[g s\textsuperscript{-1}]','interpreter','latex')
