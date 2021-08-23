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

%% Textured 3D Earth example
% Code adapted from: https://www.mathworks.com/matlabcentral/fileexchange/13823-3d-earth-example

%% Options

space_color = 'w';
npanels = 180;   % Number of globe panels around the equator deg/panel = 360/npanels
alpha   = 1; % globe transparency level, 1 = opaque, through 0 = invisible
%GMST0 = []; % Don't set up rotatable globe (ECEF)
GMST0 = 4.89496121282306; % Set up a rotatable globe at J2000.0

% Earth texture image
% Anything imread() will handle, but needs to be a 2:1 unprojected globe
% image.

image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';

% Mean spherical earth

erad    = 1e-3*6371008.7714; % [km] equatorial radius
prad    = 1e-3*6371008.7714; % [km] polar radius 
erot    = 7.2921158553e-5; % earth rotation rate (radians/sec)

%% Create figure

% h1=figure()
% set(h1,'Color', space_color, 'units','normalized','outerposition',[0 0 1 1]);

hold on;

% Turn off the normal axes

% set(gca, 'NextPlot','add', 'Visible','off');

axis equal;
axis auto;

% Set initial view

view(0,30);

axis vis3d;

%% Create wireframe globe

% Create a 3D meshgrid of the sphere points using the ellipsoid function

[x, y, z] = ellipsoid(0, 0, 0, erad, erad, prad, npanels);

globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);

if ~isempty(GMST0)
    hgx = hgtransform;
    set(hgx,'Matrix', makehgtform('zrotate',GMST0));
    set(globe,'Parent',hgx);
end

%% Texturemap the globe

% Load Earth image for texture map

cdata = imread(image_file);

% Set image as color data (cdata) property, and set face color to indicate
% a texturemap, which Matlab expects to be in cdata. Turn off the mesh edges.

% set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', alpha, 'EdgeColor', 'none');

cdata_time = cdata;
idx_time_earth_example = mod( floor(size(cdata,2)*time_earth_example/86400), size(cdata,2)) + 1;
if idx_time_earth_example == 1
    cdata_time = cdata;
else
    %     cdata_time(:,1:size(cdata,2)-idx_time_earth_example,:) = cdata(:,idx_time_earth_example+1:size(cdata,2),:);
    %     cdata_time(:,size(cdata,2)-idx_time_earth_example+1:size(cdata,2),:) = cdata(:,1:idx_time_earth_example,:);
    
    cdata_time(:,1:idx_time_earth_example,:) = cdata(:,size(cdata,2)-idx_time_earth_example+1:size(cdata,2),:);
    cdata_time(:,idx_time_earth_example+1:size(cdata,2),:) = cdata(:,1:size(cdata,2)-idx_time_earth_example,:);
end

set(globe, 'FaceColor', 'texturemap', 'CData', cdata_time, 'FaceAlpha', alpha, 'EdgeColor', 'none');
