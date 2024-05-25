% Create parameters and random users
Users = 20;
W = 20; % Bandwidth Mhz
Pdis = 20; % BS Sector Power Watios
P = Pdis/W; % Spectral power density W/Mhz
Radio = 1000; % Meters


k = 1.3806488e-23;
T = 290; % Temperature in Kelvin
F = 3.6e9; % Hz
c = 3e8; % Meters/second
lamda = c/F; % Meters
azimuth = 2*pi*rand(Users,1); % Radians
radio = Radio*rand(Users,1); % Meters
[x,y] = pol2cart(azimuth,radio);

figure
plot(x,y,'.')
hold on

viscircles([0 0],1000);

% Divide the cell in sectors
AreaUsers = zeros(Users,1);
rhosArea = cell(5,1);
for user = 1:Users
 if (azimuth(user) < 2*pi/5)
 AreaUsers(user) = 1;
 rhosArea{1} = [rhosArea{1} radio(user)];
 elseif ((azimuth(user) < 4*pi/5) && (azimuth(user) > 2*pi/5))
 AreaUsers(user) = 2;
 rhosArea{2} = [rhosArea{2} radio(user)];
 elseif((azimuth(user) < 6*pi/5) && (azimuth(user) > 4*pi/5))
 AreaUsers(user) = 3;
 rhosArea{3} = [rhosArea{3} radio(user)];
 elseif ((azimuth(user) < 8*pi/5) && (azimuth(user) > 6*pi/5))
 AreaUsers(user) = 4;
 rhosArea{4} = [rhosArea{4} radio(user)];
 elseif ((azimuth(user) < 10*pi/5) && (azimuth(user) > 8*pi/5))
 AreaUsers(user) = 5;
 rhosArea{5} = [rhosArea{5} radio(user)];
 end
end

% Sort the users in each sector by distance to the BS
rhosArea{1} = sort(rhosArea{1});
rhosArea{2} = sort(rhosArea{2});
rhosArea{3} = sort(rhosArea{3});
rhosArea{4} = sort(rhosArea{4});
rhosArea{5} = sort(rhosArea{5});

% Count users per sector
Users1 = 0;
Users2 = 0;
Users3 = 0;
Users4 = 0;
Users5 = 0;
for user = 1:Users
 if (AreaUsers(user) == 1)
 Users1 = Users1 + 1;
 elseif (AreaUsers(user) == 2)
 Users2 = Users2 + 1;
 elseif (AreaUsers(user) == 3)
 Users3 = Users3 + 1;
 elseif (AreaUsers(user) == 4)
 Users4 = Users4 + 1;
 elseif (AreaUsers(user) == 5)
 Users5 = Users5 + 1;
 end
end

UsersArea = [Users1 Users2 Users3 Users4 Users5];

% Calculate power coefficients for each area/sector (Max. 10users per sector) and bandwidth per user
% (OFDMA)
powers = cell(5,1);
bandwidths = zeros(5,1);
for area = 1:5
 if (UsersArea(area) == 1)
 power = Pdis/W;
 powers{area} = power;
 bandwidths(area) = W;
 
 elseif (UsersArea(area) == 2)
 syms P1 P2;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20;
 eqs = [eq1 eq2];
 vars = [P1 P2];
 [sol1, sol2] = solve(eqs, vars);
 power = [sol1 sol2];
 powers{area} = power;
 bandwidths(area) = W/2;
 
 elseif (UsersArea(area) == 3)
 syms P1 P2 P3;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20;
 eq3 = P3 == 2*(P2+P1);
 eqs = [eq1 eq2 eq3];
 vars = [P1 P2 P3];
 [sol1, sol2, sol3] = solve(eqs, vars);
 power = [sol1 sol2, sol3];
 powers{area} = power;
 bandwidths(area) = W/3;
 
 elseif (UsersArea(area) == 4)
 syms P1 P2 P3 P4;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eqs = [eq1 eq2 eq3 eq4];
 vars = [P1 P2 P3 P4];
 [sol1, sol2, sol3, sol4] = solve(eqs, vars);
 power = [sol1 sol2 sol3 sol4];
 powers{area} = power;
 bandwidths(area) = W/4;
 
 elseif (UsersArea(area) == 5)
 syms P1 P2 P3 P4 P5;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20 + P5*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eq5 = P5 == 2*(P1+P2+P3+P4);
 eqs = [eq1 eq2 eq3 eq4 eq5];
 vars = [P1 P2 P3 P4 P5];
 [sol1, sol2, sol3, sol4, sol5] = solve(eqs, vars);
 power = [sol1 sol2 sol3 sol4 sol5];
 powers{area} = power;
 bandwidths(area) = W/5;
 
 elseif (UsersArea(area) == 6)    
 syms P1 P2 P3 P4 P5 P6;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20 + P5*20 + P6*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eq5 = P5 == 2*(P1+P2+P3+P4);
 eq6 = P6 == 2*(P1+P2+P3+P4+P5);
 eqs = [eq1 eq2 eq3 eq4 eq5 eq6];
 vars = [P1 P2 P3 P4 P5 P6];
 [sol1, sol2, sol3, sol4, sol5, sol6] = solve(eqs, vars);
 power = [sol1 sol2 sol3 sol4 sol5 sol6];
 powers{area} = power;
 bandwidths(area) = W/6;
 
 elseif (UsersArea(area) == 7)
 syms P1 P2 P3 P4 P5 P6 P7;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20 + P5*20 + P6* 20 + P7*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eq5 = P5 == 2*(P1+P2+P3+P4);
 eq6 = P6 == 2*(P1+P2+P3+P4+P5);
 eq7 = P7 == 2*(P1+P2+P3+P4+P5+P6);
 eqs = [eq1 eq2 eq3 eq4 eq5 eq6 eq7];
 vars = [P1 P2 P3 P4 P5 P6 P7];
 [sol1, sol2, sol3, sol4, sol5, sol6, sol7] = solve(eqs,vars);
 power = [sol1 sol2 sol3 sol4 sol5 sol6 sol7];
 powers{area} = power;
 bandwidths(area) = W/7;
 
 elseif (UsersArea(area) == 8)
 syms P1 P2 P3 P4 P5 P6 P7 P8;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20 + P5*20 + P6*20 + P7*20 + P8*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eq5 = P5 == 2*(P1+P2+P3+P4);
 eq6 = P6 == 2*(P1+P2+P3+P4+P5);
 eq7 = P7 == 2*(P1+P2+P3+P4+P5+P6);
 eq8 = P8 == 2*(P1+P2+P3+P4+P5+P6+P7);
 eqs = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8];
 vars = [P1 P2 P3 P4 P5 P6 P7 P8];
 [sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8] =solve(eqs, vars);
 power = [sol1 sol2 sol3 sol4 sol5 sol6 sol7 sol8];
 powers{area} = power;
 bandwidths(area) = W/8;
 
 elseif (UsersArea(area) == 9)
 syms P1 P2 P3 P4 P5 P6 P7 P8 P9;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20 + P5*20 + P6*20 + P7*20 +P8*20 + P9*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eq5 = P5 == 2*(P1+P2+P3+P4);
 eq6 = P6 == 2*(P1+P2+P3+P4+P5);
 eq7 = P7 == 2*(P1+P2+P3+P4+P5+P6);
 eq8 = P8 == 2*(P1+P2+P3+P4+P5+P6+P7);
 eq9 = P9 == 2*(P1+P2+P3+P4+P5+P6+P7+P8);
 eqs = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9];
 vars = [P1 P2 P3 P4 P5 P6 P7 P8 P9];
 [sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9] =solve(eqs, vars);
 power = [sol1 sol2 sol3 sol4 sol5 sol6 sol7 sol8 sol9];
 powers{area} = power;
 bandwidths(area) = W/9;
 
 elseif (UsersArea(area) == 10)
 syms P1 P2 P3 P4 P5 P6 P7 P8 P9 P10;
 eq1 = P2 == 2*P1;
 eq2 = Pdis == P1*20 + P2*20 + P3*20 + P4*20 + P5*20 + P6*20 + P7*20 +P8*20 + P9*20 + P10*20;
 eq3 = P3 == 2*(P2+P1);
 eq4 = P4 == 2*(P3+P2+P1);
 eq5 = P5 == 2*(P1+P2+P3+P4);
 eq6 = P6 == 2*(P1+P2+P3+P4+P5);
 eq7 = P7 == 2*(P1+P2+P3+P4+P5+P6);
 eq8 = P8 == 2*(P1+P2+P3+P4+P5+P6+P7);
 eq9 = P9 == 2*(P1+P2+P3+P4+P5+P6+P7+P8);
 eq10 = P10 == 2*(P1+P2+P3+P4+P5+P6+P7+P8+P9);
 eqs = [eq1 eq2 eq3 eq4 eq5 eq6 eq7 eq8 eq9 eq10];
 vars = [P1 P2 P3 P4 P5 P6 P7 P8 P9 P10];
 [sol1, sol2, sol3, sol4, sol5, sol6, sol7, sol8, sol9,sol10] = solve(eqs, vars);
 power = [sol1 sol2 sol3 sol4 sol5 sol6 sol7 sol8 sol9
sol10];
 powers{area} = power;
 bandwidths(area) = W/10;
 end
end

% Calculate SNR and capacity per user, then calculate total BS capacity
Cnoma = 0;
Cofdma = 0;
N = k*T;
for area = 1:5
 for user = 1:length(powers{area})
 SNR = P/(N*(4*pi*rhosArea{area}(user)/lamda)^2);
 Cofdma = Cofdma + bandwidths(area)*10^6*log2(1+SNR);
 if (user == 1)
 I = 0;
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 2)
 I = powers{area}(user-1);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 3)
 
 I = powers{area}(user-1) + powers{area}(user-2);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 4)    
 I = powers{area}(user-1) + powers{area}(user-2) + powers{area}(user-3);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 elseif (user == 5)
     
 I = powers{area}(user-1) + powers{area}(user-2) +powers{area}(user-3) + powers{area}(user-4);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 6)
 I = powers{area}(user-1) + powers{area}(user-2) + powers{area}(user-3) + powers{area}(user-4) + powers{area}(user-5);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 7)
 I = powers{area}(user-1) + powers{area}(user-2) +powers{area}(user-3) + powers{area}(user-4) + powers{area}(user-5) + powers{area}(user-6);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 8)
 I = powers{area}(user-1) + powers{area}(user-2) +powers{area}(user-3) + powers{area}(user-4) + powers{area}(user-5) + powers{area}(user-6) + powers{area}(user-7);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 
 elseif (user == 9)
 I = powers{area}(user-1) + powers{area}(user-2) +powers{area}(user-3) + powers{area}(user-4) + powers{area}(user-5) + powers{area}(user-6) + powers{area}(user-7) +powers{area}(user-8);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*log2(1+SNR));
 
 elseif (user == 10)
 I = powers{area}(user-1) + powers{area}(user-2) +powers{area}(user-3) + powers{area}(user-4) + powers{area}(user-5) + powers{area}(user-6) + powers{area}(user-7) +powers{area}(user-8) + powers{area}(user-9);
 SNR = (powers{area}(user)/(4*pi*rhosArea{area}(user)/lamda)^2)/(N+I/(4*pi*rhosArea{area}(user)/lamda)^2);
 Cnoma = Cnoma + (W*10^6*log2(1+SNR));
 end
 end
end

BaseTotalCapacityWithNOMA = double(Cnoma) 
BaseTotalCapacityWithOFDMA = double(Cofdma)


% Static power consumption (assumed to be 20 Watts as per your previous message)
P_static = 100; 

% Computing EE for NOMA
EE_NOMA = BaseTotalCapacityWithNOMA / (P_static + Pdis);  % Energy Efficiency for NOMA in bits/Joule

% Computing EE for OFDMA
EE_OFDMA = BaseTotalCapacityWithOFDMA / (P_static + Pdis);  % Energy Efficiency for OFDMA in bits/Joule

% Display results
fprintf('Energy Efficiency (EE) for NOMA: %e bits/Joule\n', EE_NOMA);
fprintf('Energy Efficiency (EE) for OFDMA: %e bits/Joule\n', EE_OFDMA);

% User count data
%power=20
%bandwidths=20
user_counts = [5, 10, 20, 35];

% Capacity data for NOMA and OFDMA
noma_capacity = [
    mean([3.0615e+09, 2.8598e+09, 2.8643e+09]),
    mean([2.8707e+09, 3.6496e+09, 2.8883e+09]),
    mean([3.8874e+09, 3.8286e+09, 4.0285e+09]),
    mean([3.9013e+09, 3.9132e+09, 4.2354e+09])
];
ofdma_capacity = [
    mean([3.0171e+09, 2.8511e+09, 2.8563e+09]),
    mean([2.8048e+09, 3.5271e+09, 2.8072e+09]),
    mean([3.6161e+09, 3.5435e+09, 3.6948e+09]),
    mean([3.6209e+09, 3.6376e+09, 3.7197e+09])
];

% Energy efficiency data for NOMA and OFDMA
noma_ee = [
    mean([2.551226e+07, 2.383190e+07, 2.386890e+07]),
    mean([2.392228e+07, 3.041352e+07, 2.406889e+07]),
    mean([3.239531e+07, 3.190500e+07, 3.357045e+07]),
    mean([3.251094e+07, 3.261000e+07, 3.529529e+07])
];
ofdma_ee = [
    mean([2.514267e+07, 2.375920e+07, 2.380282e+07]),
    mean([2.337369e+07, 2.939243e+07, 2.339312e+07]),
    mean([3.013380e+07, 2.952901e+07, 3.078960e+07]),
    mean([3.017399e+07, 3.031354e+07, 3.099747e+07])
];

% Plot the graphs
figure;

% Capacity vs. Number of Users
subplot(2, 1, 1);
plot(user_counts, noma_capacity, 'o-', 'LineWidth', 2);
hold on;
plot(user_counts, ofdma_capacity, 's-', 'LineWidth', 2);
xlabel('Number of Users');
ylabel('Capacity (bits/Joule)');
title('Capacity vs. Number of Users');
legend('NOMA Capacity', 'OFDMA Capacity');
grid on;

% Energy Efficiency vs. Number of Users
subplot(2, 1, 2);
plot(user_counts, noma_ee, 'o-', 'LineWidth', 2);
hold on;
plot(user_counts, ofdma_ee, 's-', 'LineWidth', 2);
xlabel('Number of Users');
ylabel('Energy Efficiency (bits/Joule)');
title('Energy Efficiency vs. Number of Users');
legend('NOMA EE', 'OFDMA EE');
grid on;

% Base Station power levels
%bandwidths=20
%user=20
powers = [1, 10, 50, 100];  

% Capacity data for NOMA and OFDMA
noma_capacity = [
    mean([3.5254e+09, 3.5002e+09, 3.5068e+09]),
    mean([3.7734e+09, 3.7644e+09, 3.9810e+09]),
    mean([4.1533e+09, 4.1703e+09, 4.1106e+09]),
    mean([4.2526e+09, 3.9837e+09, 4.4377e+09])
];

ofdma_capacity = [
    mean([3.1726e+09, 3.2609e+09, 3.1715e+09]),
    mean([3.5113e+09, 3.4990e+09, 3.5529e+09]),
    mean([3.7690e+09, 3.8798e+09, 3.8841e+09]),
    mean([3.9551e+09, 3.7413e+09, 4.0562e+09])
];

% Energy efficiency data for NOMA and OFDMA
noma_ee = [
    mean([3.490537e+07, 3.465567e+07, 3.472109e+07]),
    mean([3.430374e+07, 3.422205e+07, 3.619072e+07]),
    mean([2.768852e+07, 2.780220e+07, 2.740380e+07]),
    mean([2.126295e+07, 1.991841e+07, 2.218862e+07])
];

ofdma_ee = [
    mean([3.141195e+07, 3.228619e+07, 3.140091e+07]),
    mean([3.192135e+07, 3.180882e+07, 3.229866e+07]),
    mean([2.512652e+07, 2.586525e+07, 2.589414e+07]),
    mean([1.977573e+07, 1.870670e+07, 2.028096e+07])
];

% Plot the graphs
figure;

% Capacity vs. Power of Base Station
subplot(2,1,1);
plot(powers, noma_capacity, 'o-', 'LineWidth', 2);
hold on;
plot(powers, ofdma_capacity, 's-', 'LineWidth', 2);
xlabel('Power of BS (W)');
ylabel('Capacity (bits/Joule)');
title('Capacity vs. Power of BS');
legend('NOMA Capacity', 'OFDMA Capacity');
grid on;

% Energy Efficiency vs. Power of Base Station
subplot(2,1,2);
plot(powers, noma_ee, 'o-', 'LineWidth', 2);
hold on;
plot(powers, ofdma_ee, 's-', 'LineWidth', 2);
xlabel('Power of BS (W)');
ylabel('Energy Efficiency (bits/Joule)');
title('Energy Efficiency vs. Power of BS');
legend('NOMA EE', 'OFDMA EE');
grid on;

% Bandwidth data
%user=20
%power=20
bandwidths = [1.4, 10, 20, 50];

% Capacity data for NOMA and OFDMA
noma_capacity_data = [
    mean([2.7641e+08, 2.7925e+08, 2.6000e+08]),  % 1.4 MHz
    mean([1.9527e+09, 1.9169e+09, 1.8392e+09]),  % 10 MHz
    mean([3.9416e+09, 4.1281e+09, 3.9502e+09]),  % 20 MHz
    mean([9.8373e+09, 9.1528e+09, 1.0367e+10])   % 50 MHz
];
ofdma_capacity_data = [
    mean([2.7990e+08, 2.7719e+08, 2.7165e+08]),  % 1.4 MHz
    mean([1.8700e+09, 1.8498e+09, 1.8145e+09]),  % 10 MHz
    mean([3.6733e+09, 3.7090e+09, 3.6392e+09]),  % 20 MHz
    mean([8.6450e+09, 8.5781e+09, 8.8019e+09])   % 50 MHz
];

% Energy efficiency data for NOMA and OFDMA
noma_ee_data = [
    mean([2.303418e+06, 2.327102e+06, 2.166677e+06]),  % 1.4 MHz
    mean([1.627281e+07, 1.597386e+07, 1.532654e+07]),  % 10 MHz
    mean([3.284674e+07, 3.440051e+07, 3.291832e+07]),  % 20 MHz
    mean([8.197770e+07, 7.627293e+07, 8.639070e+07])   % 50 MHz
];
ofdma_ee_data = [
    mean([2.332466e+06, 2.309948e+06, 2.263779e+06]),  % 1.4 MHz
    mean([1.558320e+07, 1.541510e+07, 1.512046e+07]),  % 10 MHz
    mean([3.061085e+07, 3.090858e+07, 3.032667e+07]),  % 20 MHz
    mean([7.204201e+07, 7.148421e+07, 7.334955e+07])   % 50 MHz
];

% Plot the graphs for Capacity and Energy Efficiency
figure;

% Capacity vs. Bandwidth
subplot(2, 1, 1);
plot(bandwidths, noma_capacity_data, 'o-', 'LineWidth', 2);
hold on;
plot(bandwidths, ofdma_capacity_data, 's-', 'LineWidth', 2);
xlabel('Bandwidth (MHz)');
ylabel('Capacity (bits/Joule)');
title('Capacity vs. Bandwidth');
legend('NOMA Capacity', 'OFDMA Capacity');
grid on;

% Energy Efficiency vs. Bandwidth
subplot(2, 1, 2);
plot(bandwidths, noma_ee_data, 'o-', 'LineWidth', 2);
hold on;
plot(bandwidths, ofdma_ee_data, 's-', 'LineWidth', 2);
xlabel('Bandwidth (MHz)');
ylabel('Energy Efficiency (bits/Joule)');
title('Energy Efficiency vs. Bandwidth');
legend('NOMA EE', 'OFDMA EE');
grid on;

% Define the number of users and their corresponding improvements
users = [5, 10, 20, 35];
improvements = [0.7, 3, 8.3, 9.7];  % in percent

% Convert percentages to a decimal for plotting
improvements_decimal = improvements;

% Plot the improvements for NOMA over OFDMA
figure;
plot(users, improvements_decimal, '-o', 'DisplayName', 'Improvement of NOMA over OFDMA');
title('Improvement Number of Users');
xlabel('Number of Users');
ylabel('Improvement (bit/Joule)%');
legend show;
ylim([0 15]);

% Define the power and their corresponding improvements
power = [1, 10, 50, 100];
improvements = [9.4, 9, 7.8, 7.7];  % in percent

% Convert percentages to a decimal for plotting
improvements_decimal = improvements;

% Plot the improvements for NOMA over OFDMA
figure;
plot(power, improvements_decimal, '-o', 'DisplayName', 'Improvement of NOMA over OFDMA');
title('Improvement power');
xlabel('power');
ylabel('Improvement (bit/Joule)%');
legend show;
ylim([0 15]);
% Define the bandwidths and their corresponding improvements
bandwidths = [1.4, 10, 20, 50];
improvements = [-1.7, 2.79, 8.8, 12.7];  % in percent

% Convert percentages to a decimal for plotting
improvements_decimal = improvements;

% Plot the improvements for NOMA over OFDMA
figure;
plot(bandwidths, improvements_decimal, '-o', 'DisplayName', 'Improvement of NOMA over OFDMA');
title('Improvement bandwidths');
xlabel('bandwidths');
ylabel('Improvement (bit/Joule)%');
legend show;
ylim([-2 15]);