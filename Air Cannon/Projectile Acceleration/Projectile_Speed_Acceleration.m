clear all
clc
clearvars
close all


L=80*2.54/100; %Barrel Length [m]
R=6.1/2*2.54/100; %Chamber Internal Diameter [m]
W=14*2.54/100; %Cylinder Length [m]
Piston=11.5*2.54/100; %Piston Length [m]
D_rod=(1+3/8)*2.54/100; %Piston Rod Diameter [m]
D_FrontPlate=4*2.54/100; %front Plate Diameter [m]
t_FrontPlate=1.5*2.54/100; %Front Plate Tickness [m]
t_Piston=1.5*2.54/100; %Piston Tickness [m]
W_Back=W-Piston; %Back Chamber Length
W_Front=Piston-t_Piston; %Front Chamber Length
V_back=pi*R^2*W_Back; % Back Chamber Volume
V_Front=pi*R^2*W_Front-pi*(D_rod/2)^2*(Piston-t_Piston-t_FrontPlate)-pi*(D_FrontPlate/2)^2*t_FrontPlate; %Front Chamber Volume %V_0
V_Tot=pi*R^2*W; %Chamber Total Volume
d_Barrel=2*2.54/100; %Barrel Internal Diamter [m] %r
A=pi*(0.04267/2)^2; %Golf Ball Cross Section [m]
%A=pi*(d_Barrel/2)^2; %Barrel Cross Section [m]
P_0=294400; %220*6894.76; %Initial Pressure in the Chamber [Pa]
P_atm=100000; %Atmospheric Pressure [Pa]
m=0.5; %Projectile Mass [kg]

f=0.2; %Friction Factor between PTFE and Steel

x=0.1:0.1:10;
v_out_is=sqrt((2/m)*(P_0*V_Front*log(1+A*L/V_Front)-A*L*P_atm));
%v_is=sqrt((2/m)*(P_0*V_Front*log(1+A*x/V_Front)-A*x*P_atm)); %outlet velocity with isoentropic expansion
V_Vec=0:0.0005:V_Tot;
V_Vec(1)=0.0001;
v_out_is_Vol=sqrt((2/m)*(P_0*V_Vec.*log(1+A*L./V_Vec)-A*L*P_atm));

% figure
% plot(V_Vec,v_out_is_Vol,[V_Front,V_Front],[10,70],'k',[V_Tot,V_Tot],[10,70],'r')


p=180000:10000:2500000;
v_is=sqrt((2/m)*(p*V_Front*log(1+A*L/V_Front)-A*L*P_atm)); %outlet velocity with isoentropic expansion
v_is_frict=sqrt((2/m)*(p*V_Front*log(1+A*L/V_Front)-A*L*P_atm-10*f*m*9.81)); %outlet velocity with isoentropic expansion with friction
gamma=7/5; 
v_ad=sqrt((2/m)*((p*V_Front)/(gamma-1)*(1-(V_Front/(A*L+V_Front))^(gamma-1))-A*L*P_atm)); %outlet velocity with adiabatic expansion
v_ad_frict=sqrt((2/m)*((p*V_Front)/(gamma-1)*(1-(V_Front/(A*L+V_Front))^(gamma-1))-A*L*P_atm-10*f*m*9.81)); %outlet velocity with adiabatic expansion with friction

% figure
% plot(p*10^-5,v_is,'k',p*10^-5,v_ad, 'b',p*10^-5,v_is_frict,'--k',p*10^-5,v_ad_frict, '--b')
% %eventually add experimental points
% hold on
% p_exp=[300000,350000,400000,500000];
% dist=5*2.54/100;
% v_exp=[1.025*dist/6.177e-3,1.025*dist/5.070e-3,1.025*dist/4.54e-3,1.025*dist/3.803e-3];
% plot(p_exp*10^-5, v_exp, '*r')
% hold off

%Plot of the Velocity and Accelaration as function of the position 
%NO FRICTION

x=0:0.01:2;
P_Front=[2,2.5,3,5,10,15].*1e5;
time=zeros(length(P_Front),length(x));

for i_press=1:length(P_Front)
    
    Vel_Vector_ad(i_press,:)=sqrt((2/m)*((P_Front(i_press)*V_Front)/(gamma-1)*(1-(V_Front./(A*x+V_Front)).^(gamma-1))-A*x*P_atm));
    Acc_Vector_ad(i_press,:)=(1/m)*(A*(P_Front(i_press)*(V_Front.^(gamma))./((V_Front+A*x).^(gamma))-P_atm));
    
    Inverse_V_fun= @(xx) 1./(sqrt((2/m)*((P_Front(i_press)*V_Front)/(gamma-1)*(1-(V_Front./(A*xx+V_Front)).^(gamma-1))-A*xx*P_atm)));
    
    for i_temp=2:length(x)
        time(i_press,i_temp)=integral(Inverse_V_fun,0,x(i_temp));
    end
    
end

figure
plot(x,Vel_Vector_ad(1,:),x,Vel_Vector_ad(2,:),x,Vel_Vector_ad(3,:),x,Vel_Vector_ad(4,:),x,Vel_Vector_ad(5,:),x,Vel_Vector_ad(6,:));
title('Projectile Velocity vs Position')
xlabel('Position [m]')
ylabel('Velocity [m/s]')
legend('2bar','2.5bar','3bar','5bar','10bar','15bar')
saveas(gcf,'Velocity_vs_Position.png')

figure
plot(x,Acc_Vector_ad(1,:),x,Acc_Vector_ad(2,:),x,Acc_Vector_ad(3,:),x,Acc_Vector_ad(4,:),x,Acc_Vector_ad(5,:),x,Acc_Vector_ad(6,:));
title('Projectile Acceleration vs Position')
xlabel('Position [m]')
ylabel('Acceleration [m/s^2]')
legend('2bar','2.5bar','3bar','5bar','10bar','15bar')
saveas(gcf,'Acceleration_vs_Position.png')

figure
plot(x,Vel_Vector_ad(1,:),x,Vel_Vector_ad(2,:),x,Vel_Vector_ad(3,:));
title('Projectile Velocity vs Position')
xlabel('Position [m]')
ylabel('Velocity [m/s]')
legend('2bar','2.5bar','3bar')
saveas(gcf,'Velocity_vs_Position_Zoom.png')

figure
plot(x,Acc_Vector_ad(1,:),x,Acc_Vector_ad(2,:),x,Acc_Vector_ad(3,:));
title('Projectile Acceleration vs Position')
xlabel('Position [m]')
ylabel('Acceleration [m/s^2]')
legend('2bar','2.5bar','3bar')
saveas(gcf,'Acceleration_vs_Position_Zoom.png')


%% Plotting as function of time

% P_Front=3e5;
% 
% for 
% 
% Inverse_V_fun()= @(xx) 1./(sqrt((2/m)*((P_Front*V_Front)/(gamma-1)*(1-(V_Front./(A*xx+V_Front)).^(gamma-1))-A*xx*P_atm)));
% 
% for i_temp=1:length(x)
%     time(i_temp)=integral(Inverse_V_fun,0,x(i_temp));
% end
% 
% time(:,1)=0;


figure
plot(time(1,:),Vel_Vector_ad(1,:),time(2,:),Vel_Vector_ad(2,:),time(3,:),Vel_Vector_ad(3,:),time(4,:),Vel_Vector_ad(4,:),time(5,:),Vel_Vector_ad(5,:),time(6,:),Vel_Vector_ad(6,:));
title('Projectile Velocity vs Time')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('2bar','2.5bar','3bar','5bar','10bar','15bar')
saveas(gcf,'Velocity_vs_Time.png')

figure
plot(time(1,:),Acc_Vector_ad(1,:),time(2,:),Acc_Vector_ad(2,:),time(3,:),Acc_Vector_ad(3,:),time(4,:),Acc_Vector_ad(4,:),time(5,:),Acc_Vector_ad(5,:),time(6,:),Acc_Vector_ad(6,:));
title('Projectile Acceleration vs Time')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
legend('2bar','2.5bar','3bar','5bar','10bar','15bar')
saveas(gcf,'Acceleration_vs_Time.png')

figure
plot(time(1,:),Vel_Vector_ad(1,:),time(2,:),Vel_Vector_ad(2,:),time(3,:),Vel_Vector_ad(3,:));
title('Projectile Velocity vs Time')
xlabel('Time [s]')
ylabel('Velocity [m/s]')
legend('2bar','2.5bar','3bar')
saveas(gcf,'Velocity_vs_Time_Zoom.png')

figure
plot(time(1,:),Acc_Vector_ad(1,:),time(2,:),Acc_Vector_ad(2,:),time(3,:),Acc_Vector_ad(3,:));
title('Projectile Acceleration vs Time')
xlabel('Time [s]')
ylabel('Acceleration [m/s^2]')
legend('2bar','2.5bar','3bar')
saveas(gcf,'Acceleration_vs_Time_Zoom.png')
