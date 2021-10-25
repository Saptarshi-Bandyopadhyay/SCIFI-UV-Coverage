%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 by California Institute of Technology.  ALL RIGHTS RESERVED. %
% United  States  Government  sponsorship  acknowledged.   Any commercial use %
% must   be  negotiated  with  the  Office  of  Technology  Transfer  at  the %
% California Institute of Technology.                                         %
%                                                                             %
% This software may be subject to  U.S. export control laws  and regulations. %
% By accepting this document,  the user agrees to comply  with all applicable %
% U.S. export laws and regulations.  User  has the responsibility  to  obtain %
% export  licenses,  or  other  export  authority  as may be required  before %
% exporting  such  information  to  foreign  countries or providing access to %
% foreign persons.                                                            %
%                                                                             %
% This  software  is a copy  and  may not be current.  The latest  version is %
% maintained by and may be obtained from the Mobility  and  Robotics  Sytstem %
% Section (347) at the Jet  Propulsion  Laboratory.   Suggestions and patches %
% are welcome and should be sent to the software's maintainer.                %
%                                                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plotting code

color_array = ['c' 'r' 'b' 'g' 'm' 'y' 'k'];
marker_array = ['o' 's' 'd' '^' 'v' '>' '<' 'p' 'h' '+'];
standard_font_size = 25;
title_font_size = 40;

plot_handle = figure(1);
clf
set(plot_handle,'Color',[1 1 1]);
set(plot_handle,'units','normalized','Position',[0 0 1 1])
set(plot_handle,'PaperPositionMode','auto');

subplot(2,2,1)
hold on
plot3(0,0,0,'ok','MarkerSize',10,'MarkerFaceColor','b')
plot3(Sun_Pos_Velocity(1,k), Sun_Pos_Velocity(2,k), Sun_Pos_Velocity(3,k), 'hk','MarkerSize',20,'MarkerFaceColor','y')
plot3(Moon_Pos_Velocity(1,k), Moon_Pos_Velocity(2,k), Moon_Pos_Velocity(3,k), 'ok','MarkerSize',10,'MarkerFaceColor','w')

plot3(1.2*AU*Target_x, 1.2*AU*Target_y, 1.2*AU*Target_z, 'hk','MarkerSize',20,'MarkerFaceColor','r')

plot3(Sun_Pos_Velocity(1,1:k),Sun_Pos_Velocity(2,1:k),Sun_Pos_Velocity(3,1:k),'-y','LineWidth',2)

legend('Earth','Sun','Moon','Target Direction')
axis equal
view(0,90)
xlabel('X axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
ylabel('Y axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
zlabel('Z axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
title('Trajectory of Sun, Moon around Earth (ECI Frame)','fontsize',title_font_size,'FontName','Times New Roman')
set(gca, 'fontsize',standard_font_size,'FontName','Times New Roman')
hold off

k_start = max(k - ceil(2*86400/delta_t), 1);
subplot(2,2,2)
legendCell = [];
hold on
for ns = 1:1:num_SC
    plot3(SC_orbit_data{ns}.pos_earth_array(k,1)-SC_orbit_data{1}.pos_earth_array(k,1), SC_orbit_data{ns}.pos_earth_array(k,2)-SC_orbit_data{1}.pos_earth_array(k,2), SC_orbit_data{ns}.pos_earth_array(k,3)-SC_orbit_data{1}.pos_earth_array(k,3), 'ok','MarkerSize',10,'MarkerFaceColor',color_array(mod(ns,length(color_array))+1) )
    legendCell{ns} = ['SC ',num2str(ns)];
end
for ns = 1:1:num_SC
    plot3(SC_orbit_data{ns}.pos_earth_array(k_start:k,1)-SC_orbit_data{1}.pos_earth_array(k_start:k,1), SC_orbit_data{ns}.pos_earth_array(k_start:k,2)-SC_orbit_data{1}.pos_earth_array(k_start:k,2), SC_orbit_data{ns}.pos_earth_array(k_start:k,3)-SC_orbit_data{1}.pos_earth_array(k_start:k,3), '-','LineWidth',2,'Color', color_array(mod(ns,length(color_array))+1) )
end
legend(legendCell)
axis equal
xlabel('X axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
ylabel('Y axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
zlabel('Z axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
title('Trajectory of SC around chief (Chief Centered Inertial Frame)','fontsize',title_font_size,'FontName','Times New Roman')
set(gca, 'fontsize',standard_font_size,'FontName','Times New Roman')
hold off

subplot(2,2,3)
hold on
time_earth_example = time_array(k); % [sec]
earth_example

for ns = 1:1:num_SC
    plot3(SC_orbit_data{ns}.pos_earth_array(k,1), SC_orbit_data{ns}.pos_earth_array(k,2), SC_orbit_data{ns}.pos_earth_array(k,3), 'ok','MarkerSize',10,'MarkerFaceColor',color_array(mod(ns,length(color_array))+1) )
    legendCell{ns} = ['SC ',num2str(ns)];
end
for ns = 1:1:num_SC
    plot3(SC_orbit_data{ns}.pos_earth_array(k_start:k,1), SC_orbit_data{ns}.pos_earth_array(k_start:k,2), SC_orbit_data{ns}.pos_earth_array(k_start:k,3), '-','LineWidth',2,'Color', color_array(mod(ns,length(color_array))+1) )
end

lla = ecef2lla(SC_orbit_data{1}.pos_earth_array(k,1:3));
% view([lla(2)+90,30])
% view([10,30])

legend(legendCell)
axis equal
title('Trajectory of SC around Earth','fontsize',title_font_size,'FontName','Times New Roman')
xlabel('X axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
ylabel('Y axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
zlabel('Z axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
set(gca, 'fontsize',standard_font_size,'FontName','Times New Roman')
hold off

subplot(2,2,4)
plot(reshape(UV_pos(1:k,:,1),[],1), reshape(UV_pos(1:k,:,2),[],1),'.k')
axis equal
title('UV coverage of Target','fontsize',title_font_size,'FontName','Times New Roman')
xlabel('X axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
ylabel('Y axis [km]','fontsize',standard_font_size,'FontName','Times New Roman')
set(gca, 'fontsize',standard_font_size,'FontName','Times New Roman')
hold off


this_date_str_t0 = date_str_t0;
this_date_str_t0.Second = this_date_str_t0.Second + time_array(k);
this_date_str = datestr(this_date_str_t0);
sgtitle(['Simulation Time = ',num2str(round(time_array(k)/86400, 2, 'significant')),' days, Date = ',this_date_str(1:11)],'fontsize',title_font_size,'FontName','Times New Roman')


