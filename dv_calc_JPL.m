%% dv budget
Time = 4*7*24*3600; %one month

s3_const = s3_constants;
J2 = s3_const.J2_EARTH;
GM = s3_const.GM.EARTH;
Re = s3_const.R_EARTH;
a = 7000e3;
n = sqrt(GM/a^3);
eta = 1;
i = deg2rad(60);
ade_max = 1600; %max ade vector length
ade_ave = 800; %average ade vector length

rho =5.63e-13;
Cd = 1.4;
m = 12;
dA = 0.3*0.2*3*(1-cosd(10)); %differental max cross area
dBd = Cd*dA/m;
dadot = rho*a*n*dBd; %drag drift

drift_ada_tot = a*dadot*Time;
dv_da = drift_ada_tot/2*n;


phi_dot = 0.75*J2*Re^2*sqrt(GM)/(a^3.5*eta^4)*(5*cos(i)^2-1); %rotate rate of ade
drift_ade_max = phi_dot*Time*ade_max;
drift_ade_ave = phi_dot*Time*ade_ave;
dv_de_max = drift_ade_max/2*n;
dv_de_ave = drift_ade_ave/2*n;