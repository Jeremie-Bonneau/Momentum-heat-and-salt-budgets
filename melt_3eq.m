function [a,T_b,S_b,a_error] = melt_3eq(Cd,G_T,G_S,U_w,T_w,S_w,P,precision)
%Function to calculate the melting using Jenkins 2010 3 equations model.
%Although here we use the full TEOS-10 equation for freezing point instead 
% of a simplified version.

% V1.2 Jérémie Bonneau   May 10 2023
% jeremie_bonneau@outlook.com

% INPUT: 
% Cd   : Drag coefficient
% G_T  : Transfer coefficient for heat
% G_S  : Transfer coefficient for salt 
% U_W  : Far field water velocity magnitude
% T_W  : Far field water temperature
% S_W  : Far field water salinity
% P    : Pressure at which the ice-ocean interaction happens (depth)
% Precision : choice of 1, 2 or 3. This dictates the size of the array of
% possible boundary layer salinity values. 1=37, 2=370, 3=3700. Recommend
% using 3 unless big calculations. Takes about 1 min to compute 1M values
% on one core at precision 3...

% OUTPUT: 
% a   : melt rate [m/day]
% T_S : Boundary layer temperature [degC]
% S_b : Boundary layer salinity [g/kg]
% a_error: Max error due to discrete S_b values [m/day]

% Equation for freezing point
% 1) T_f=f(S,P)    freezing point, we use TEOS-10
% 
% Equation for heat transfer
% 2) rho_i*(a*L + c_i*a*(T_b-T_i)) = rho_w*c_w*Cd^0.5*U_w*G_T*(T_b-T_b)
% 
% Equation for salt transfer
% 3) rho_i*a*(S_b-S_i) = rho_w*Cd^0.5*U_w*G_S*(S_b-S_b)
% We use S_i=0
%
% Variables
% T   : Temperature [degC]
% S   : Salinity [g/kg]
% P   : Pressure [dbar]
% a   : melting rate [m/s in eq, m/d as output]
%
% i index is for ice, b is for boundary layer, w is for far field water
% (outside the boundary layer), f is for freezing. 
% 
% Coefficients
% L      : Latent heat of fusion [3.3e5 J/kg]
% c_i    : Heat capacity of ice [2000 J/kg/degC]
% T_i    : Ice temperature [-10 degC]  Could be modified for colder ice
% rho_i  : Heat capacity of seawater [920 kg/m^3]
% rho_w  : density of seawater [1030 kg/m^3]
% Cd     : Drag coefficient : user input
% G_T    : Transfer coefficient for heat: user input 
% G_S    : Transfer coefficient for salt: user input 
% 
% Three unknowns: a, T_b and S_b, three equations so ok
% Technically, T_b is at freezing temperature, so T_b in(2) can be replaced
% by T_f. So what we are doing here is to compute the freezing temperature
% for the possible salinity in the boundary layer (Sarray) and calculate
% the resulting melt rate of (2) and (3). The salinity for which both melt
% rates (from (2) and (3)) are equal is the right boundary layer salinity.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example (uncomment and paste in command window run or in script)
% Tm=[-1.8:0.1:4]; % array of possible Temp
% Sm=[25:0.1:35];  % array of possible Salinity
% Cd=0.0025;  %value from Cowton 2015
% G_T=0.022;  %value from Cowton 2015
% G_S=0.00062;%value from Cowton 2015
% P=40; %at ~40m depth
% precision=3;
% 
% %Different velocities to try
% U(1)=0.002;
% U(2)=0.02;
% U(3)=0.2;
% 
% [Tm,Sm]=meshgrid(Tm,Sm);
% 
% % calculations
% am=nan(length(Tm(:,1)),length(Tm(1,:)),3);
% for k=1:3
%     for i=1:length(Tm(:,1))
%         for j=1:length(Tm(1,:))
%             [a,~,~,~]=melt_3eq(Cd,G_T,G_S,U(k),Tm(i,j),Sm(i,j),P,precision);
%             am(i,j,k)=a;
%         end
%     end
% end
% 
% figure; 
% subplot(1,3,1)
% contourf(Sm,Tm,squeeze(am(:,:,1)));
% title('Melting [m/d] for U_w=0.002') ; colorbar
% subplot(1,3,2)
% contourf(Sm,Tm,squeeze(am(:,:,2)));
% title('Melting [m/d] for U_w=0.02') ; colorbar
% subplot(1,3,3)
% contourf(Sm,Tm,squeeze(am(:,:,3)));
% title('Melting [m/d] for U_w=0.2') ; colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Calculations

% Constants
c_w=4000;
L=330000;
c_i=2000;
T_i=-10;
rho_w=1030; %seawater
rho_i=920; %ice

% have to modify if S is possibly>37 g/kg
if precision==1
    Sarray=[0:37];
elseif precision==2
    Sarray=[0:0.1:37];
elseif precision==3
    Sarray=[0:0.01:37];
elseif precision==4
    Sarray=[0:0.001:37];
end

T_f=gsw_CT_freezing_poly(Sarray,P); %BDL temp according to TEOS-10
a1=((c_w.*Cd^0.5.*U_w.*G_T.*(T_w-T_f)))./((L+c_i.*(T_f-T_i)));  % melting according to heat (1)
a2=((Cd^0.5.*U_w.*G_S.*(S_w-Sarray)))./(Sarray);  % melting according to salt (2) (S_i=0)
[~,idx]=min((a1-a2).^2); % finding at what value (index) of Sarray the system is closed 
a=(a1(idx)+a2(idx))/2.*3600.*24; % in m per day
a_error=abs(a1(idx)-a2(idx)).*3600.*24;
T_b=T_f(idx);
S_b=Sarray(idx);


end