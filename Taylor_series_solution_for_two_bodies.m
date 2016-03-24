%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg



% this program is used for showing a satalite about the earth by tailor serise
close all; clear all; clc;
%% constants
G=6.67428E-11; % gravitational constant
open=1;  dt_old=0;
%% Intial condition
M=[5.972*10^24,7.3477e22];       % [M1,M2]
mu=G*(sum(M));
R1=[0;0;0];                 % position of M(1)
R2=[3.844e8;0;0];                 % position of M(2)
I=[0,0,0];              % location of initial axis
V1=[0;0;0];                 % velocity of M(1)
V2=[0.5;0.5;0]*1.02306e3;                 % velocity of M(2)
RT=R2';
VT=V2';
%% orbit parameter
tf=365*24*360*0.3;   % final time of soution
dt=360;            % time step
Dt=[360,36000];  % [min dt, max dt] for slider
%% RK4 parameter
X0=[R1;R2;V1;V2];
B=[0;0;0;0;0;0;0;0;0;0;0;0];
order=12;
fig=figure();
menu = uicontrol('Parent',fig,'Style','popupmenu','String',{'dt';'Exit'},'Units','centimeters' ,'Position',[17.5,0.25,3,0.5]);
slider = uicontrol('Parent',fig,'Style','slider','Units','centimeters' ,'Position',[0,0,10,0.5],'value',Dt(1),'SliderStep', [0.01,0.1] , 'min',Dt(1), 'max',Dt(2));
set(gcf,'Color','w');
while open==1
    figure(fig);
    S = get(menu,'value');
    if S == 1
        dt=get(slider,'value');
    elseif S == 2
        open=0;
    end
    if dt~=dt_old
        clear sol RT VT R1_x R1_y R1_z R2_x R2_y R2_z V1_x V1_y V1_z V2_x V2_y V2_z error;
        cla;
        sol(1:12,1)=X0;
        RT=R2';
        VT=V2';
        %% solution by RK4
        for n=1:length(0:dt:tf)
            b=G*M(2)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
            c=-G*M(1)/(norm(sol(1:3,n)-sol(4:6,n)))^3;
            A=[0,0,0,0,0,0,1,0,0,0,0,0; ...
                0,0,0,0,0,0,0,1,0,0,0,0; ...
                0,0,0,0,0,0,0,0,1,0,0,0; ...
                0,0,0,0,0,0,0,0,0,1,0,0; ...
                0,0,0,0,0,0,0,0,0,0,1,0; ...
                0,0,0,0,0,0,0,0,0,0,0,1;...
                -b,0,0,b,0,0,0,0,0,0,0,0; ...
                0,-b,0,0,b,0,0,0,0,0,0,0; ...
                0,0,-b,0,0,b,0,0,0,0,0,0; ...
                -c,0,0,c,0,0,0,0,0,0,0,0; ...
                0,-c,0,0,c,0,0,0,0,0,0,0; ...
                0,0,-c,0,0,c,0,0,0,0,0,0 ];
            [ XX ] = RK4( A,B,sol(1:12,n),dt,n*dt,(n+1)*dt,order );
            sol(1:12,n+1)=XX(1:12,2);
        end
        R1_x=sol(1,:);
        R1_y=sol(2,:);
        R1_z=sol(3,:);
        R2_x=sol(4,:);
        R2_y=sol(5,:);
        R2_z=sol(6,:);
        V1_x=sol(7,:);
        V1_y=sol(8,:);
        V1_z=sol(9,:);
        V2_x=sol(10,:);
        V2_y=sol(11,:);
        V2_z=sol(12,:);
        %% solution by four terms of taylor series
        tl4=dt;
        for l4=2:length(0:dt:tf)
            zero=RT(l4-1,1:3);
            first=VT(l4-1,1:3);
            second=-mu*RT(l4-1,1:3)/(norm(RT(l4-1,1:3)))^3;
            third=0*mu/(norm(RT(l4-1,1:3)))^4*(-norm(RT(l4-1,1:3))*VT(l4-1,1:3)+3*RT(l4-1,1:3)*norm(VT(l4-1,1:3)));
            fourth=-mu*((norm(RT(l4-1,1:3)))^3.*second-3*(norm(RT(l4-1,1:3)))^2*norm(VT(l4-1,1:3))*VT(l4-1,1:3))/(norm(RT(l4-1,1:3)))^6+3*mu*((norm(RT(l4-1,1:3)))^4*(norm(second)*RT(l4-1,1:3)+norm(VT(l4-1,1:3))*VT(l4-1,1:3))-4*(norm(RT(l4-1,1:3)))^3*(norm(VT(l4-1,1:3)))^2*RT(l4-1,1:3))/(norm(RT(l4-1,1:3)))^8;
            RT(l4,1:3)=zero+first*tl4+second*tl4^2/2+third*tl4^3/6+fourth*tl4^4/24;
            VT(l4,1:3)=first+second*tl4+third*tl4^2/2+fourth*tl4^3/6;
        end
    end
        %% plotting
        % --------------------------------------------------------------------------------------------------------------------------------------------------------
        subplot(5,1,1:3);
        cla;
        % axes at M1
        view(3);
        hold all;
        % M1 position
        plot3(0,0,0,'o','color','cyan','LineWidth',5);
        % M2 position relative to M1
        plot3(R2_x-R1_x,R2_y-R1_y,R2_z-R1_z,'color','g','LineWidth',2);
        % Taylor solution
        plot3(RT(:,1),RT(:,2),RT(:,3),'color','b','LineWidth',1)
        grid on;
        xlabel('X','Fontsize',18);
        ylabel('Y','Fontsize',18);
        zlabel('Z','Fontsize',18);
        title(['Solution at M1 Center of Mass'  '  dt = ' num2str(dt)],'Fontsize',18);
        legend('M1 position','RK45 solution','Taylor solution');
        dt_old=dt;
%--------------------------------------------------------------------------------------------------------------------------------------------------------
subplot(5,1,5);
%% error plotting
for ttt=1:length(R1_x)-1
    error(ttt)=abs(norm([R2_x(ttt),R2_y(ttt),R2_z(ttt)])-norm(RT(ttt,:)))/norm([R2_x(ttt),R2_y(ttt),R2_z(ttt)]);
end
% figure(2);
% set(gcf,'Color','w');
plot(0:1:length(R1_x)-2,error*100);
grid on;
xlabel('R index','Fontsize',18);
ylabel('|error %|','Fontsize',18);
legend(['Max.(|Error %|) = ' num2str(max(error*100)) ' Mean(|Error %|) = ' num2str(mean(error*100))  ' Std(|Error %|) = ' num2str(std(error*100)) ],'location','northwest');
end
close;