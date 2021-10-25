close all
clear all

%file name to save gif
filename = 'V_JPL_power_separation';

%%% V
phi_vec=[-90 -90 -90 -90 -90 -90 90 90 90 90 90 90]; %phase of relative e vector
adi_mag=[50 100 200 400 800 1600  50 100 200 400 800 1600]; %magnitude of i vector
ade_mag= [50 100 200 400 800 1600 50 100 200 400 800 1600]; %magnitude of relative e vector
nu=[-90 -90 -90 -90 -90 -90 -90 -90 -90 -90 -90 -90]; %relative i vector phase
alambda=[0 0 0 0 0 0 0 0 0 0 0 0];
ada = [0 0 0 0 0 0 0 0 0 0 0 0 ];

%%%skew eivec
% phi_vec=[-96 -84 -72 -108 -60 -120 96 84 72 108 60 120]; %phase of relative e vector
% adi_mag=[50 100 200 400 800 1600  50 100 200 400 800 1600]; %magnitude of i vector
% ade_mag= [50 100 200 400 800 1600 50 100 200 400 800 1600]; %magnitude of relative e vector
% nu=[90 90 90 90 90 90 -90 -90 -90 -90 -90 -90]; %relative i vector phase
% alambda=[0 0 0 0 0 0 0 0 0 0 0 0];
% ada = [0 0 0 0 0 0 0 0 0 0 0 0 ];

%%% power in-plane separation
%  ade_mag = [200 200 200 200 400 400 400 400 800 800 800 800];
%  phi_vec = [0 90 180 270 45 135 225 315 0 90 180 270];
% adi_mag=zeros(1,12); %magnitude of i vector
% nu=zeros(1,12); %relative i vector phase
% alambda=zeros(1,12);
% ada = zeros(1,12);









if length(phi_vec) ~= length(adi_mag) ||...
    length(adi_mag) ~= length(ade_mag) ||...
    length(ade_mag) ~= length(nu) ||...
    length(nu) ~= length(alambda)
    fprintf('WARNING: VECTORS NOT EQUAL \n')
end
    

for i=1:length(phi_vec)
    count = 0;
    for k=0:0.001:360
        count = count + 1;
        
        ade_mag_curr=ade_mag(i);
        nu_curr=nu(i);
        phi_curr=phi_vec(i);
        adi_mag_curr=adi_mag(i);
        alambda_curr=alambda(i);
        
        r(i,count)=-ade_mag_curr*cosd(k-phi_curr)+ada(i);
        t(i,count)= alambda_curr+2*ade_mag_curr*sind(k-phi_curr);
        n(i,count)=adi_mag_curr*sind(k-nu_curr);
    end
end

minsep=1000;
count = 0;
for k=0:0.001:360
    count = count + 1;
    
    for index1=1:length(phi_vec);
        r1vec=[r(index1,count), t(index1,count), n(index1,count)];
        for index2=(index1+1):length(phi_vec)
            r2vec=[r(index2,count), t(index2,count), n(index2,count)];
            sep_curr=norm(r1vec-r2vec);
            if sep_curr<minsep
                minsep=sep_curr;
                indexmin=[index1 index2 count];
            end
        end
    end
end
fprintf('Min Separation: %f \n', minsep)



close all
tr=figure();
hold on
xlabel('T (m)')
ylabel('R (m)')
nr=figure();
hold on
xlabel('N (m)')
ylabel('R (m)')
p3D=figure();
hold on
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')
set(gca,'FontSize',20)
plot3(0,0,0,'o','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
grid on




for i=1:length(phi_vec)
    figure(tr)
    plot(t(i,:),r(i,:),'k')
    figure(nr)
    plot(n(i,:),r(i,:),'k')
    figure(p3D)
    plot3(r(i,:),t(i,:),n(i,:),'k')
end

all_view = figure();
for i=1:length(phi_vec)
    figure(all_view)
    hold on
    subplot(1,3,1)
    axis equal
    hold on
    plot(t(i,:),r(i,:),'k')
    xlabel('T (m)')
    ylabel('R (m)')
    subplot(1,3,2)
    axis equal
    hold on
    xlabel('N (m)')
    ylabel('R (m)')
    plot(n(i,:),r(i,:),'k')
    subplot(1,3,3)
    axis equal
    hold on
    xlabel('R (m)')
    ylabel('T (m)')
    zlabel('N (m)')
    plot3(r(i,:),t(i,:),n(i,:),'k')
    view(3)
end
% axis equal

figure(nr)
% axis equal
%3D plotting
figure()
axis equal on
grid on
hold on
for i=1:length(phi_vec)
    plot3(r(i,:),t(i,:),n(i,:), 'k');
end
view(3)
xlabel('R (m)')
ylabel('T (m)')
zlabel('N (m)')


figure()
axis equal on
grid on
for i=1:length(phi_vec)
    hold on
    plot(ade_mag(i)*cosd(phi_vec(i)),ade_mag(i)*sind(phi_vec(i)), 'o','MarkerFaceColor','k');
end
xlabel('a\delta e_x (m)')
ylabel('a\delta e_y (m)')
ylim([-900 900])


figure()
axis equal on
grid on
hold on
for i=1:length(phi_vec)
    plot(adi_mag(i)*cosd(nu(i)),adi_mag(i)*sind(nu(i)), 'o','MarkerFaceColor','k');
end
xlabel('a\delta i_x (m)')
ylabel('a\delta i_y (m)')
ylim([-900 900])



%3D plotting
figure(p3D)
axis equal on

for i=1:length(phi_vec)
    p{i}=plot3(r(i,1),t(i,1),n(i,1), 'o','MarkerFaceColor','k');
end
view(3)

movie_fig = 8;
figure(movie_fig)
axis equal
hold on
set(gca,'FontSize',20)
xlabel('R [m]')
ylabel('T [m]')
zlabel('N [m]')
xmin = min(min(r));
ymin = min(min(t));
zmin = min(min(n));
xmax = max(max(r));
ymax = max(max(t));
zmax = max(max(n));
xlim([xmin-0.1*abs(xmin) xmax+0.1*abs(xmax)])
ylim([ymin-0.1*abs(ymin) ymax+0.1*abs(ymax)])
zlim([zmin-0.1*abs(zmin) zmax+0.1*abs(zmax)])

axis manual
 for i=1:length(phi_vec)
     plot3(r(i,:),t(i,:), n(i,:),'k')
 end

 for i=1:length(phi_vec)
     p{i}=plot3(r(i,1),t(i,1),n(i,1), 'o','MarkerFaceColor','k', 'MarkerEdgeColor','k');
 end
 p{i+1}=plot3(0,0,0, 'o','MarkerFaceColor','k', 'MarkerEdgeColor','k');

view(3)
grid on
F = getframe(gcf);
im = frame2im(F);
[imind, cm] = rgb2ind(im,256);
imwrite(imind, cm, filename, 'gif', 'Loopcount',inf, 'DelayTime',0.00000001)


for index=1:360
    figure(movie_fig)
    for i=1:length(phi_vec)
        p{i}.XData = r(i,index);
        p{i}.YData = t(i,index);
        p{i}.ZData = n(i,index);
    end
    drawnow
    F = getframe(gcf);
im = frame2im(F);
[imind, cm] = rgb2ind(im,256);
imwrite(imind, cm, filename, 'gif', 'WriteMode','append', 'DelayTime',0.00000001)
end
    
i=5;

save('SC_RTN.mat','r','t','n')