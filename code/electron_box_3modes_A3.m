function electron_box_3modes_A3(mode)
%mode=3;
%% Parameters

%General Program Parameters
global G;
G.Mode=6;
G.Mode2=mode;
% %Mode2=1: Part 1
% %Mode2=2: Part 2
% %Mode2=3: Part 2

%num_elec=10000;
num_elec=10000;
num_elec_plot=10;%64;
pause_time=0.001;
multiplier=5;
old_plot_elec=zeros(num_elec,2);

%Assignment Parameters
m0=9.10938356e-31; % [kg]
mn=0.26*m0; % [kg]
kB=1.38064852e-23; % [J*K^-1]
Temp=300; % [K]
Tau_mn=0.2e-12; % [s]
q_electron=1.60217662e-19; % [C]
elec_concentration=1e19; % [m^-2]
global P;
P.box_size_x=200e-9; % [m]
P.box_size_y=100e-9; % [m]
P.thermal_speed=sqrt(kB*Temp/mn); %FIXME: is for 1d -- may be right
P.time_step=(P.box_size_y/(multiplier*100))/P.thermal_speed;
P.distancemax=6*(P.box_size_y/(multiplier*100));
P.Pscat=1-exp(-P.time_step/Tau_mn);
P.MFP=sqrt(2)*P.thermal_speed*Tau_mn;
P.V_part1=0.1; % [V]
P.E_part1=P.V_part1/P.box_size_x; % []
P.F_part1=q_electron*P.E_part1; % []
P.a_part1=P.F_part1/mn;
P.x_bins=200;
P.y_bins=100;
P.size_div=P.box_size_y/P.y_bins;
P.electron_per_particle=elec_concentration.*P.box_size_x.*P.box_size_y./num_elec;
P.particle_charge=P.electron_per_particle.*q_electron;

% Matrix Names
global S;
S.px=1; S.py=2; S.vx=3; S.vy=4; S.num_dt=5; %S.ax=6; S.ay=7;


% Matrix to Track path lengths
scat_info=zeros(1,2);

%% Pre-Print
if(G.Mode2==2)
    % 1.a)
    fprintf('The electric field is: %4.3s V/m.\n',P.E_part1);
    % 1.b)
    fprintf('The force on each electron is: %4.3s N.\n',P.F_part1);
    % 1.c)
    fprintf('The accelation of each electron is: %4.3s m/s^2.\n',P.a_part1);
end

%% Get Electric Field and Acceleration
if(G.Mode2==2)
    [Ex_scale,Ey_scale]=box_potential_1a_A3(0.1);
elseif(G.Mode2==3)
    [Ex_scale,Ey_scale]=box_potential_2a_A3(0.8);
end
P.ax=q_electron.*Ex_scale./mn;
P.ay=q_electron.*Ey_scale./mn;


%% Generate Boxes
boxes=zeros(1,5);
boxes(1,:)=P.box_size_y.*[-0.5 1 3 1 0]; %top boundary
boxes(2,:)=P.box_size_y.*[-0.5 -1 3 1 0]; %bottom boundary

if(G.Mode2==3)
    boxes(3,:)=P.box_size_y.*[0.8 -0.1 0.4 0.5 1/P.box_size_y];% x0 y0 width height is_diffuse --- (x0 y0) is bottom left
    boxes(4,:)=P.box_size_y.*[0.8 0.6 0.4 0.5 1/P.box_size_y];% x0 y0 width height is_diffuse --- (x0 y0) is bottom left
    %boxes(4,:)=P.box_size_y.*[0.8 0.6 0.4 0.5 0];% x0 y0 width height is_diffuse --- (x0 y0) is bottom left
end


%% Create matrix of electrons
electrons=zeros(num_elec,5);
for i1=1:num_elec
    valid=false;
    while(~valid)
        electrons(i1,:)=create_electon;
        valid=true;%assume valid after creating electron
        for ii=1:size(boxes,1)
            e_y=electrons(i1,S.py);
            e_x=electrons(i1,S.px);
            x0=boxes(ii,1);
            y0=boxes(ii,2);
            x1=boxes(ii,1)+boxes(ii,3);
            y1=boxes(ii,2)+boxes(ii,4);
            inside_box=(e_y>y0)&(e_y<y1)&(e_x>x0)&(e_x<x1);
            valid=~inside_box&valid;
        end
    end
    
end


%% Run for over time and plot
if(G.Mode2==1)
%     figure(1)
%     title('Figure 1: Particle Trajectories for Simple Model')
%     xlabel('x position (m)')
%     ylabel('y position (m)')
elseif(G.Mode2==2)
    figure(1)
    title('Figure 1: Particle Trajectories')
    xlabel('x position (m)')
    ylabel('y position (m)')
elseif(G.Mode2==3)
    figure(7)
    title('Figure 7: Particle Trajectories with Bottle Neck')
    xlabel('x position (m)')
    ylabel('y position (m)')
end
hold on;  


%plot boxes
if(G.Mode2==3)
    for ii=1:size(boxes,1)
        x0=boxes(ii,1);
        y0=boxes(ii,2);
        x1=boxes(ii,1)+boxes(ii,3);
        y1=boxes(ii,2)+boxes(ii,4);
        plot([x0 x1 x1 x0 x0],[y0 y0 y1 y1 y0],'Color','black');
    end
end

%initialize matrix that keeps track of limited number of electrons
for eii=1:num_elec_plot
    old_plot_elec(eii,1)=electrons(eii,S.px);
    old_plot_elec(eii,2)=electrons(eii,S.py);
end
color_map=hsv(num_elec_plot);

%loop over time
for t=1:1000*multiplier 
%for t=1:10*multiplier %FIXME: set back time loop to 1000*multiplier
   time=t*P.time_step;
   [electrons,scat_info]=move_electrons(electrons,boxes,scat_info);

   %Calculate Temperature
   v_speed=sqrt(electrons(:,S.vx).^2+electrons(:,S.vy).^2);
   vx_MS=sum(electrons(:,S.vx).^2)/size(electrons,1);
   vy_MS=sum(electrons(:,S.vy).^2)/size(electrons,1);   
   vt_MS=vx_MS+vy_MS;
   Temp_Calc=mn*vt_MS/(2*kB);
   
   Time_Temp(t)=Temp_Calc;
   
   %Calculate Current
   %Current_Temp(t)=sum(electrons(:,S.vx).*P.particle_charge)/P.box_size_x;
   Current_Temp(t)=mean(electrons(:,S.vx)).*elec_concentration.*q_electron.*P.box_size_y;
   
   %Plot '2-D plot of particle trajectories'
       if(mod(t,multiplier)==0)
           for eii=1:num_elec_plot
                loop_around=abs(old_plot_elec(eii,1)-electrons(eii,S.px))>(P.box_size_x/2);
                old_plot_elec(eii,1)=~loop_around.*old_plot_elec(eii,1)+loop_around.*electrons(eii,S.px);
                plot([old_plot_elec(eii,1) electrons(eii,S.px)],...
                     [old_plot_elec(eii,2) electrons(eii,S.py)],...
                     'LineStyle','-','Color',color_map(eii,:))

                old_plot_elec(eii,1)=electrons(eii,S.px);
                old_plot_elec(eii,2)=electrons(eii,S.py);
           end
       end


       xlim([0 P.box_size_x])
       ylim([0 P.box_size_y])
   

   pause(pause_time)

   
end
hold off; 

%% Perform Final Evaluation
if(G.Mode2==2)

    
    % 1.d)
    figure(2)
    plot(P.time_step:P.time_step:t*P.time_step,Current_Temp)
    title('Figure 2: Current as a Function of Time')
    xlabel('Time (s)')
    ylabel('Current (A)')
    grid on;
    
    % 1.e)
    % Density
    figure(3)
    hist3([electrons(:,S.px),electrons(:,S.py)],[100 50],'CdataMode','auto')
    view(0,90)
    title('Figure 3: Electron Density Map')
    xlabel('x position (m)')
    ylabel('y position (m)')
    cb3=colorbar;
    cb3.Label.String = 'Number of Electrons';
    xlim([0 P.box_size_x])
    ylim([0 P.box_size_y])
    
    % Temperature
    size_div=P.box_size_y/50;
    numb_x=round(P.box_size_x/size_div)+1; numb_y=round(P.box_size_y/size_div)+1;
    temp_map=zeros(numb_x,numb_y,2);
    for ee=1:size(electrons,1)
        x_round=ceil(electrons(ee,S.px)./size_div)+1;
        y_round=ceil(electrons(ee,S.py)./size_div)+1;
        temp_map(x_round,y_round,1)=temp_map(x_round,y_round,1)+(electrons(ee,S.vx).^2)+(electrons(ee,S.vy).^2);
        temp_map(x_round,y_round,2)=temp_map(x_round,y_round,2)+1;
    end
    temp_map_val=(temp_map(:,:,1)./(temp_map(:,:,2)+0.01))*mn/(2*kB);
    
    
    figure(4)
    surf(temp_map_val');
    cb4=colorbar;
    cb4.Label.String = 'Temperature (K)';
    view(0,90)
    title('Figure 4: Temperature Map')
    xlabel('x position bin')
    ylabel('y position bin')
    xlim([2 numb_x])
    ylim([2 numb_y])

elseif(G.Mode2==3)
    figure(8)
    hist3([electrons(:,S.px),electrons(:,S.py)],[100 50],'CdataMode','auto')
    view(0,90)
    title('Figure 8: Electron Density Map with Bottle Neck')
    xlabel('x position (m)')
    ylabel('y position (m)')
    cb3=colorbar;
    cb3.Label.String = 'Number of Electrons';
    xlim([0 P.box_size_x])
    ylim([0 P.box_size_y])
    
end
% if(G.Mode2==1)
%     %fprintf('The measured average speed is: %4.3s m/s.\n',mean(v_speed));
%     fprintf('The theoretical thermal velocity is: %4.3f m/s.\n',sqrt(2)*P.thermal_speed);
%     fprintf('The theoretical mean free path is: %4.3s m.\n',P.MFP);
%     
%     figure(2)
%     plot(P.time_step:P.time_step:t*P.time_step,Time_Temp)
%     title('Figure 2: Temperature vs. Time for Fixed Velocities')
%     xlabel('Time (s)')
%     ylabel('Temperature (K)')
%     grid on;
%     
%     
% elseif(G.Mode2==2)
%     
%     % 2.a)
%     %fprintf('The measured average speed is: %4.3f m/s.\n',mean(v_speed));
%     figure(3)
%     histogram(v_speed,100)
%     xlim([0 6*P.thermal_speed])
%     ylim([0 350])
%     title('Figure 3: Speed Distribution')
%     xlabel('Speed (m/s)')
%     ylabel('Number of Electrons')
%     grid on;
%     
%     % 2.c)
%     figure(5)
%     plot(P.time_step:P.time_step:t*P.time_step,Time_Temp)
%     title('Figure 5: Temperature vs. Time for Scattering')
%     xlabel('Time (s)')
%     ylabel('Temperature (K)')
%     grid on;
%     
%     % 2.d)
%     tmp_scat_info=scat_info(2:end,:);
%     Calc_Tau_mn=mean(P.time_step.*tmp_scat_info(:,1));
%     Calc_MFP=mean(P.time_step.*tmp_scat_info(:,1).*tmp_scat_info(:,2));
%     
%     fprintf('The measured mean time between collisions is: %4.3s s.\n',Calc_Tau_mn);
%     fprintf('The measured mean free path is: %4.3s m.\n',Calc_MFP);
%     
%     
% elseif(G.Mode2==3)
%     
%     % 2.c)
%     figure(7)
%     hist3([electrons(:,S.px),electrons(:,S.py)],[100 50],'CdataMode','auto')
%     view(0,90)
%     title('Figure 7: Electron Density Map')
%     xlabel('x position (m)')
%     ylabel('y position (m)')
%     cb7=colorbar;
%     cb7.Label.String = 'Number of Electrons';
%     xlim([0 P.box_size_x])
%     ylim([0 P.box_size_y])
%     
%     % 2.d)
%     size_div=P.box_size_y/50;
%     numb_x=round(P.box_size_x/size_div)+1; numb_y=round(P.box_size_y/size_div)+1;
%     temp_map=zeros(numb_x,numb_y,2);
%     for ee=1:size(electrons,1)
%         x_round=ceil(electrons(ee,S.px)./size_div)+1;
%         y_round=ceil(electrons(ee,S.py)./size_div)+1;
%         temp_map(x_round,y_round,1)=temp_map(x_round,y_round,1)+(electrons(ee,S.vx).^2)+(electrons(ee,S.vy).^2);
%         temp_map(x_round,y_round,2)=temp_map(x_round,y_round,2)+1;
%     end
%     temp_map_val=(temp_map(:,:,1)./(temp_map(:,:,2)+0.01))*mn/(2*kB);
%     
%     
%     figure(8)
%     surf(temp_map_val');
%     cb8=colorbar;
%     cb8.Label.String = 'Temperature (K)';
%     view(0,90)
%     title('Figure 8: Temperature Map')
%     xlabel('x position bin')
%     ylabel('y position bin')
%     xlim([2 numb_x])
%     ylim([2 numb_y])
%     
%     
% end


end %FIXME: make function

%% Function create_electon
function [electron]=create_electon
    global G;
    global S;
    global P;
    

    electron(1,S.px)=rand()*P.box_size_x;
    electron(1,S.py)=rand()*P.box_size_y;

    if(G.Mode2==1)
        angle=2*pi*rand();
        speed=sqrt(2)*P.thermal_speed;
        electron(1,S.vx)=speed*cos(angle);
        electron(1,S.vy)=speed*sin(angle);
    elseif((G.Mode2==2)||(G.Mode2==3))
        electron(1,S.vx)=P.thermal_speed*randn();
        electron(1,S.vy)=P.thermal_speed*randn();
    end
    electron(1,S.num_dt)=0;

end

%% Function move_electrons
function [electrons,scat_info]=move_electrons(electrons,boxes,scat_info)
    global G;
    global S;
    global P;
    
     x_bin=ceil(electrons(:,S.px)./P.size_div);
     y_bin=ceil(electrons(:,S.py)./P.size_div);
     x_bin=x_bin.*(x_bin>=1).*(x_bin<=P.x_bins)+ones(length(x_bin),1).*(x_bin<1)...
           +P.x_bins.*ones(length(x_bin),1).*(x_bin>P.x_bins);
     y_bin=y_bin.*(y_bin>=1).*(y_bin<=P.y_bins)+ones(length(y_bin),1).*(y_bin<1)...
           +P.y_bins.*ones(length(y_bin),1).*(y_bin>P.y_bins);
     
     
     index_x=sub2ind([P.y_bins P.x_bins],[y_bin.'],[x_bin.']);
     tmp_acc_x=P.ax(index_x);
     tmp_acc_y=P.ay(index_x);
     acc_x=tmp_acc_x.';
     acc_y=tmp_acc_y.';
    
     new_py=electrons(:,S.py)+electrons(:,S.vy)*P.time_step;
     new_px=electrons(:,S.px)+electrons(:,S.vx)*P.time_step;
     
     one=1;
     distancemax=P.distancemax;
     inside_box=zeros(size(electrons,1),1);
     bottom=inside_box; top=inside_box; left=inside_box; right=inside_box; diffusive=inside_box;
     for ii=1:size(boxes,1)
        x0=boxes(ii,1);
        y0=boxes(ii,2);
        x1=boxes(ii,1)+boxes(ii,3);
        y1=boxes(ii,2)+boxes(ii,4);
        inside_box=ii.*((new_py>y0)&(new_py<y1)&(new_px>x0)&(new_px<x1)&(inside_box==0))+one.*inside_box;
            
        bottom=(inside_box==ii)&(abs(new_py-y0)<distancemax)|bottom;
        top=(inside_box==ii)&(abs(new_py-y1)<distancemax)|top;
        left=(inside_box==ii)&(abs(new_px-x0)<distancemax)|left;
        right=(inside_box==ii)&(abs(new_px-x1)<distancemax)|right;
     end
     
     diffusive=(inside_box>=1).*boxes(inside_box+double((inside_box==0)),5);
     
     
    if (G.Mode2==3) 
        b_rand_abs=abs(randn(size(electrons,1),1));
        b_rand_nor=(randn(size(electrons,1),1));
        
        bounce_vertical=(inside_box>=1)&(top|bottom);
        bounce_lateral=(inside_box>=1)&(left|right);       
        
        %x dim
        electrons(:,S.vx)=~diffusive.*(~bounce_lateral.*electrons(:,S.vx)-bounce_lateral.*electrons(:,S.vx))...
                           +diffusive.*((inside_box==0).*electrons(:,S.vx)...
                                       -bounce_lateral.*b_rand_abs.*P.thermal_speed.*(sign(electrons(:,S.vx)))...
                                       +bounce_vertical.*b_rand_nor.*P.thermal_speed)...
                           +acc_x.*P.time_step;
        electrons(:,S.px)=mod(electrons(:,S.px)+electrons(:,S.vx)*P.time_step, P.box_size_x);


        %y dim
        electrons(:,S.vy)=~diffusive.*(~bounce_vertical.*electrons(:,S.vy)-bounce_vertical.*electrons(:,S.vy))...
                          +diffusive.*((inside_box==0).*electrons(:,S.vy)...
                                       -bounce_vertical.*b_rand_abs.*P.thermal_speed.*(sign(electrons(:,S.vy)))...
                                       +bounce_lateral.*b_rand_nor.*P.thermal_speed)...
                           +acc_y.*P.time_step;
        electrons(:,S.py)=electrons(:,S.py)+electrons(:,S.vy)*P.time_step;
        
    else
        %x dim
        bounce_lateral=(inside_box>=1)&(left|right);
        electrons(:,S.vx)=~bounce_lateral.*electrons(:,S.vx)-bounce_lateral.*electrons(:,S.vx)...
                          +acc_x.*P.time_step;
        electrons(:,S.px)=mod(electrons(:,S.px)+electrons(:,S.vx)*P.time_step, P.box_size_x)...
                          +acc_y.*P.time_step;


        %y dim
        bounce_vertical=(inside_box>=1)&(top|bottom);
        electrons(:,S.vy)=~bounce_vertical.*electrons(:,S.vy)-bounce_vertical.*electrons(:,S.vy);
        electrons(:,S.py)=electrons(:,S.py)+electrons(:,S.vy)*P.time_step;
    end

    %% Scattering Functionality
    if((G.Mode2==2)||(G.Mode2==3))
        scattering=(P.Pscat>rand(size(electrons,1),1));
        rand_x=randn(size(electrons,1),1);
        rand_y=randn(size(electrons,1),1);
        
        if (G.Mode2==2)
            for ss=1:length(scattering)
                if(scattering(ss)==1)
                    num_dt=electrons(ss,S.num_dt);
                    speed_1=sqrt(electrons(ss,S.vx).^2+electrons(ss,S.vy).^2);
                    scat_info(end+1,:)=[num_dt speed_1];
                end
            end
        end
        
        electrons(:,S.num_dt)=~scattering.*(electrons(:,S.num_dt)+1)+scattering.*0;
        electrons(:,S.vx)=scattering.*P.thermal_speed.*rand_x+~scattering.*electrons(:,S.vx);
        electrons(:,S.vy)=scattering.*P.thermal_speed.*rand_y+~scattering.*electrons(:,S.vy);

    end

end
