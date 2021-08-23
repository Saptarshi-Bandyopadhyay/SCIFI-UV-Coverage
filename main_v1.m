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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                  UV Coverage for SCIFI JNEXT project                    %
%                                                                         %
% Main file to compute the UV coverage of SC in orbit.                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all, clc, clear all

%% Initialize Variables

% Time
t0_date = '21-AUG-2021 00:00:00'; % [sec] Format = [DD-MMM(words)-YYYY HH-MM-SS], If format is incorrect, cspice_str2et wont work
t0 = 0; % [sec]
delta_t = 10*60; % [sec]
total_t = 30*24*60*60; % [sec]

% Target
Target_RA = 100; % [deg] Note: RA is usually given in hr
Target_Dec = 10; % [deg]

% Spacecraft Orbits
num_SC = 5;
sma_chief = 42000; % [km]
initial_interSC_distance = 1; % [km]

% Obstruction Variables
keep_out_angle_Sun = 30; % [deg]
keep_out_angle_Moon = 30; % [deg]

%% Initialize SPICE and Position
addpath(genpath('../mice'))
cspice_furnsh('../de421.bsp')
cspice_furnsh('../naif0012.tls')

% Time Array
time_array = [t0:delta_t:total_t]';
time_date_array = cspice_str2et(t0_date) + time_array;
date_str_t0 = datetime(t0_date);

% Moon and Sun's Position
Moon_Pos_Velocity = cspice_spkezr('MOON',time_date_array','J2000','NONE','EARTH');
Sun_Pos_Velocity = cspice_spkezr('SUN',time_date_array','J2000','NONE','EARTH');

% Astronomical Target Position
% [X, Y, Z] to [RA, Dec] -> See get_satellite_ra_dec_delta in https://github.com/Bill-Gray/sat_code/blob/master/observe.cpp
Target_z = sind(Target_Dec);
Target_x = sqrt( (1 - Target_z^2) / (1 + (tand(Target_RA))^2 ) );
if (Target_RA > 90) && (Target_RA < 270)
    Target_x = -Target_x;
end
Target_y = Target_x * tand(Target_RA);

% Constants
AU = 1.5e8; % [km]
muEarth = 398600.4415; % [km^3 sec^-2]
radiusEarth = 6400; % [km]

%% Compute Orbits

SC_orbit_data = [];

for ns = 1:1:num_SC
    
    if ns == 1
        r0 = [sma_chief 0 0]';
        v_mag = sqrt(muEarth/sma_chief);
        v0 = [0 v_mag 0]';
        
        omega_z = sqrt(muEarth/sma_chief^3); 
        
    else
        r0 = [sma_chief 0 0]' + initial_interSC_distance*2*(rand(3,1) - 0.5);
        
        % LVLH frame
        x0 = r0(1) - sma_chief;
        y0 = r0(2); 
        z0 = r0(3);
        
        x0_dot = 0;
        y0_dot = -2*omega_z*x0;
        z0_dot = 0;
        
        v_mag = sqrt(muEarth/sma_chief);
        v0 = [0 v_mag 0]' + [x0_dot y0_dot z0_dot]';
        
    end
    
    options = odeset('AbsTol', 1e-8, 'RelTol',1e-8);
    [T, XV] = ode45(@func_two_body_simple, time_array, [r0; v0], options);
    
    SC_orbit_data{ns}.pos_earth_array = XV;
    SC_orbit_data{ns}.units_pos_vel = '[km km/sec]';
end

%% Compute UV Coverage

Target_pos = 100*AU*[Target_x Target_y Target_z];

% First Check for Obstructions
obstruction_array = ones(length(time_array),1);

for k=1:1:length(time_array)
    
    for ns = 1:1:num_SC
        
        this_SC_pos = SC_orbit_data{ns}.pos_earth_array(k,1:3);
        
        % Sun's obstruction
        Sun_SC_vector = Sun_Pos_Velocity(1:3,k)' - this_SC_pos;
        Target_SC_vector = Target_pos - this_SC_pos;
        this_angle = real( acosd( dot( Sun_SC_vector/norm(Sun_SC_vector), Target_SC_vector/norm(Target_SC_vector) ) ) );
        if this_angle <= keep_out_angle_Sun
            obstruction_array(k,1) = 2*obstruction_array(k,1);
        end
        
        % Moon's obstruction
        Moon_SC_vector = Moon_Pos_Velocity(1:3,k)' - this_SC_pos;
        this_angle = real( acosd( dot( Moon_SC_vector/norm(Moon_SC_vector), Target_SC_vector/norm(Target_SC_vector) ) ) );
        if this_angle <= keep_out_angle_Moon
            obstruction_array(k,1) = 3*obstruction_array(k,1);
        end
        
        % Eclipsed by Earth
        % Use https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        x0 = [0 0 0];
        x1 = Target_pos;
        x2 = this_SC_pos;
        
        distance = norm( cross(x0-x1, x0-x2) )/ norm(x2-x1);
        
        if distance > radiusEarth
            % No eclipse
            
        else
            % Check t
            t = - dot( x1-x0, x2-x1 )/(norm(x2-x1))^2;
            if (t>=0) && (t<=1)
                % Eclipsed by Earth
                obstruction_array(k,1) = 5*obstruction_array(k,1);
            end
        end
        
    end
end

% Compute UV
UV_pos = zeros(length(time_array),num_SC*(num_SC-1),2);
for k=1:1:length(time_array)
    
    % Check for Obstruction
    if obstruction_array(k,1) == 1
        uv_counter = 0;
        
        % SC vs SC Loop
        for ns_i = 1:1:num_SC
            for ns_j = 1:1:num_SC
                
                % Only find uv-coverage for different SC
                if ns_i ~= ns_j
                    this_SCi_pos = SC_orbit_data{ns_i}.pos_earth_array(k,1:3);
                    this_SCj_pos = SC_orbit_data{ns_j}.pos_earth_array(k,1:3);
                    SCj_SCi_vector = this_SCj_pos - this_SCi_pos;
                    
                    Target_SCi_vector = Target_pos - this_SCi_pos;
                    Target_SCi_vector_normalized = Target_SCi_vector/norm(Target_SCi_vector);
                    
                    SCj_SCi_vector_component_towards_Target = dot(SCj_SCi_vector, Target_SCi_vector_normalized)*Target_SCi_vector_normalized;
                    SCj_SCi_vector_image_plane = SCj_SCi_vector - SCj_SCi_vector_component_towards_Target;
                    
                    Z_vector_component_towards_Target = dot([0 0 1], Target_SCi_vector_normalized)*Target_SCi_vector_normalized;
                    Z_vector_image_plane = [0 0 1] - Z_vector_component_towards_Target;
                    V_vector_image_plane_normalized = Z_vector_image_plane/norm(Z_vector_image_plane);
                    
                    SCj_SCi_vector_component_towards_V_vector = dot(SCj_SCi_vector_image_plane, V_vector_image_plane_normalized)*V_vector_image_plane_normalized;
                    SCj_SCi_vector_component_towards_U_vector = SCj_SCi_vector_image_plane - SCj_SCi_vector_component_towards_V_vector;
                    
                    U_vector_image_plane_normalized = cross(V_vector_image_plane_normalized, Target_SCi_vector_normalized);
                    U_vector_image_plane_normalized = U_vector_image_plane_normalized/norm(U_vector_image_plane_normalized);
                    
                    uv_counter = uv_counter + 1;
                    if real(acosd( dot(SCj_SCi_vector_component_towards_U_vector/norm(SCj_SCi_vector_component_towards_U_vector), U_vector_image_plane_normalized) )) < 90
                        UV_pos(k,uv_counter,1) = norm(SCj_SCi_vector_component_towards_U_vector);
                    else
                        UV_pos(k,uv_counter,1) = -norm(SCj_SCi_vector_component_towards_U_vector);
                    end
                    
                    if real(acosd( dot(SCj_SCi_vector_component_towards_V_vector/norm(SCj_SCi_vector_component_towards_V_vector), V_vector_image_plane_normalized) )) < 90
                        UV_pos(k,uv_counter,2) = norm(SCj_SCi_vector_component_towards_V_vector);
                    else
                        UV_pos(k,uv_counter,2) = -norm(SCj_SCi_vector_component_towards_V_vector);
                    end
                end
                
            end
        end
    end
end

%% Plotting
close all

k = length(time_array);
plot_code;

%% Video
flag_show_video = 1;
flag_store_video = 0;

if flag_show_video == 1
    close all
    disp('Plot Video')
    if (flag_store_video == 1)
        myVideo = VideoWriter('UV_coverage_video.mp4', 'MPEG-4');
        myVideo.FrameRate = 30;  % Default 30
        myVideo.Quality = 100;    % Default 75
        open(myVideo);
    end
    
    for k=1:10:length(time_array)
        plot_code;
        drawnow limitrate
        if (flag_store_video == 1)
            F = getframe(plot_handle);
            writeVideo(myVideo, F);
        end
    end
    
    if (flag_store_video == 1)
        close(myVideo);
    end
end

disp('Fin!')
