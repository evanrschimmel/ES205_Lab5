clear all
close all
clc

%% Work with calibration data

hcdata = 9:-1:1;
vcdata1 = [2.42, 2.33, 2.2, 2.1, 1.98, 1.88, 1.77, 1.66, 1.54];
fit1 = polyfit(hcdata,vcdata1, 1);
m1 = fit1(1);
int1 = fit1(2);

vcdata2 = [2.33, 2.23, 2.12, 2.02, 1.91, 1.78, 1.65, 1.55, 1.46];
fit2 = polyfit(hcdata,vcdata2, 1);
m2 = fit2(1);
int2 = fit2(2);

%% Define constants and variables

g = 32.2 * 12; % in/s^2

% Tank 1
Atank1 = 22.94; % cross sectional area of tank, in^2
d1 = 0.335; % diameter of orfice, in
Aout1 = pi * (d1/2)^2; % cross sectional area of exit orfice, in^2
a1 = 3.49; % outlet length, in
b1 = 0.5; % outlet inset, in
Cd1_orig = 0.7; % discharge coefficient, unitless
h1_0_single = 9; % initial water height, in
h1_0_double = 8; % initial water height, in

% Tank 2
Atank2 = 22.94; % cross sectional area of tank, in^2
d2 = 0.240; % diameter of orfice, in
Aout2 = pi * (d2/2)^2; % cross sectional area of exit orfice, in^2
a2 = 3.24; % outlet length, in
b2 = 0.5; % outlet inset, in
Cd2_orig = 0.7; % discharge coefficient, unitless
h2_0_single = 9; % initial water height, in
h2_0_double = 4; % initial water height, in

%% Pull experimental data and sort into vectors, run sim with original CDs

upper_data_single = readmatrix('ES205_Lab_05_Upper_Tank');
lower_data_single = readmatrix('ES205_Lab_05_Lower_Tank');
data_double = readmatrix('ES205_Lab_05_Dual_Tank');

tdata_upper = upper_data_single(:,1);
tdata_lower = lower_data_single(:,1);
tdata_double = data_double(:,1);

h1data_single = (((upper_data_single(:,2)-int1)) ./ m1);
h2data_single = (((lower_data_single(:,3)-int2)) ./ m2);
h1data_double = (((data_double(:,2)-int1)) ./ m1);
h2data_double = (((data_double(:,3)-int2)) ./ m2);

maxstep = 0.01; % max step time, s
tol = 1e-6;

tf_upper = tdata_upper(end);
sim('ES205_lab_5_upper_model_orig')

tf_lower = tdata_lower(end);
sim('ES205_lab_5_lower_model_orig')

tf_double = tdata_double(end);
sim('ES205_lab_5_dual_model_orig')

%% Run loop for best Cd
stepsize = 0.01;

dex = 1;
for Cd1_mod = 0.5:stepsize:1
    sim('ES205_lab_5_upper_model_mod')
values1(dex)=Cd1_mod;
see(dex)= sqrt(sum((h1_single_mod-h1data_single).^2)./(length(h1data_single)-2));
dex = dex+1;
end
[~,si]= min(see,[],'linear');
Cd1_mod=values1(si);
sim('ES205_lab_5_upper_model_mod')

dex = 1;
for Cd2_mod = 0.5:stepsize:1
    sim('ES205_lab_5_lower_model_mod')
values2(dex)=Cd2_mod;
see(dex)= sqrt(sum((h2_single_mod-h2data_single).^2)./(length(h2data_single)-2));
dex = dex+1;
end
[~,si]= min(see,[],'linear');
Cd2_mod=values2(si);
sim('ES205_lab_5_lower_model_mod')

sim('ES205_lab_5_dual_model_mod')

%% Plot

figure
plot(tdata_upper,h1data_single,'k-',t_upper,h1_single_orig,'r--',t_upper,h1_single_mod,'b--')
xlabel('Time [s]')
ylabel('Water height [in]')
legend('Experimental upper tank','Original upper tank model','Refined upper tank model')
set(gcf, 'color', 'w')

figure
plot(tdata_lower,h2data_single,'k-',t_lower,h2_single_orig,'r--',t_lower,h2_single_mod,'b--')
xlabel('Time [s]')
ylabel('Water height [in]')
legend('Experimental lower tank','Original lower tank model','Refined lower tank model')
set(gcf, 'color', 'w')

figure
plot(tdata_double,h1data_double,'k-',tdata_double,h2data_double,'k-',t_double,h1_double_orig,'r--',t_double,h1_double_mod,'b--',t_double,h2_double_orig,'r--',t_double,h2_double_mod,'b--')
xlabel('Time [s]')
ylabel('Water height [in]')
legend('Experimental upper tank','Experimental lower tank','Original upper tank model','Refined upper tank model','Original lower tank model','Refined lower tank model')
set(gcf, 'color', 'w')