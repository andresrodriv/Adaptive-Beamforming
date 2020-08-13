%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Wireless Communication project: Adaptive Beamforming
%Politecnico of Milan
%M.S. in Telecommunication Engineering
%Professor: Luca Reggiani
%Lab Assistant: Davide Scazzoli
%Authors:
% - Giovanni Vacirca
% - Andrés Felipe Rodrìguez Vanegas
% - Gianpaolo Alesiani
%A.Y.: 2019/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%% Parameters
%we define the parameters useful for all the script
Pars.c = physconst('LightSpeed');
Pars.fc_carrier = 1e9;                  %Carrier frequency
Pars.lambda = Pars.c/Pars.fc_carrier;   %Wavelenght

Pars.NInter = 2;%;randi([2 10]);        %Number of Interferences
for n = 1:Pars.NInter                   %Interference
    Pars.fc_Inter(n) = 2e9;%randi([1e9 6e9]);
    Pars.lambda_Inter(n) = Pars.c/Pars.fc_Inter(n); %Wavelenght
end

Pars.NTx_ant = 5;%randi([2 10]);        %Number of Antennas in the Base station
Pars.antenna = phased.IsotropicAntennaElement('BackBaffled',false);%here we create an antenna

Pars.BS_power = 10;
Pars.txPower = 90;              % Transmitter power [watt]
Pars.interfPower = 90;          % Interferences power [watt]
Pars.SNR = 20;                  % SNR awgn noise [dB]
%% Geometry
%Base station
Geometry.BS_array = phased.URA('Element',Pars.antenna',...
    'Size',[Pars.NTx_ant Pars.NTx_ant],...
    'ElementSpacing',[Pars.lambda/2,Pars.lambda/2], 'ArrayNormal', 'x');% creation of the URA of antennas
Geometry.BS_Antpos = getElementPosition(Geometry.BS_array);
Geometry.Rx_Pos = [0,0,25]; 

%Vehicles
%%Vehicle 1
%Geometry.V1_initPos = [randi([50 100]),randi([-100 -50]),1.5];
Geometry.V1_initPos = [50,-100,25];
Geometry.V1_finalPos = [Geometry.V1_initPos(1),Geometry.V1_initPos(2)+200,25];
Geometry.TrackV1 = sqrt(sum((Geometry.V1_finalPos(1,:)-Geometry.V1_initPos(1,:)).^2));%Trajectory distance
Geometry.Dist_V1(1) = sqrt(sum((Geometry.V1_initPos(1:2)-Geometry.Rx_Pos(1:2)).^2));%Distance of the Vehicule to the BS
%Caculate DoA for V1 at start
Geometry.AOA_V1(1) = atan2(Geometry.V1_initPos(1,2)-Geometry.Rx_Pos(1,2),...
    Geometry.V1_initPos(1,1)-Geometry.Rx_Pos(1,1))*180/pi;
Geometry.ZOA_V1(1) = -atan2(Geometry.Rx_Pos(1,3) - Geometry.V1_initPos(1,3),Geometry.Dist_V1)*180/pi;
Geometry.DOA_V1(1,:) = [Geometry.AOA_V1(1) Geometry.ZOA_V1(1)];
%DoA at each point of the trajectory
for i=1:Geometry.TrackV1
    PosV1 = Geometry.V1_initPos(1,:)+[0,i,0];
    Geometry.Dist_V1(i+1) = sqrt(sum((PosV1(1:2)-Geometry.Rx_Pos(1:2)).^2));
    Geometry.AOA_V1(i+1) = atan2(PosV1(1,2)-Geometry.Rx_Pos(1,2),PosV1(1,1)-Geometry.Rx_Pos(1,1))*180/pi;
    Geometry.ZOA_V1(i+1) = -atan2(Geometry.Rx_Pos(1,3) - PosV1(1,3),Geometry.Dist_V1(i+1))*180/pi;
    Geometry.DOA_V1(i+1,:) = [Geometry.AOA_V1(i+1) Geometry.ZOA_V1(i+1)];
end

%%Vehicle 2
%Geometry.V2_initPos = [randi([100 150]),randi([-50 0]),1.5];
Geometry.V2_initPos = [100,-50,25];
Geometry.V2_finalPos = [Geometry.V2_initPos(1)-200,Geometry.V2_initPos(2),25];
Geometry.TrackV2 = sqrt(sum((Geometry.V2_finalPos(1,:)-Geometry.V2_initPos(1,:)).^2));%Trajectory distance
Geometry.Dist_V2(1) = sqrt(sum((Geometry.V2_initPos(1:2)-Geometry.Rx_Pos(1:2)).^2));%Distance of the Vehicule to the BS
%Caculate DoA for V2 at start
Geometry.AOA_V2(1) = atan2(Geometry.V2_initPos(1,2)-Geometry.Rx_Pos(1,2),...
    Geometry.V2_initPos(1,1)-Geometry.Rx_Pos(1,1))*180/pi;
Geometry.ZOA_V2(1) = -atan2(Geometry.Rx_Pos(1,3) - Geometry.V2_initPos(1,3),Geometry.Dist_V2)*180/pi;
Geometry.DOA_V2(1,:) = [Geometry.AOA_V2(1) Geometry.ZOA_V2(1)];
%DoA at each point of the trajectory
for i=1:Geometry.TrackV2
    PosV2 = Geometry.V2_initPos(1,:)-[i,0,0];
    Geometry.Dist_V2(i+1) = sqrt(sum((PosV2(1:2)-Geometry.Rx_Pos(1:2)).^2));
    Geometry.AOA_V2(i+1) = atan2(PosV2(1,2)-Geometry.Rx_Pos(1,2),PosV2(1,1)-Geometry.Rx_Pos(1,1))*180/pi;
    Geometry.ZOA_V2(i+1) = -atan2(Geometry.Rx_Pos(1,3) - PosV2(1,3),Geometry.Dist_V2(i+1))*180/pi;
    Geometry.DOA_V2(i+1,:) = [Geometry.AOA_V2(i+1) Geometry.ZOA_V2(i+1)];
end

%Interferences
Geometry.IPos(1,:) = [100,-100,25];
Geometry.IPos(2,:) = [100,100,25];
for n = 1:Pars.NInter
    %Geometry.IPos(n,:) = [randi([-150 200]),randi([-200 100]),1.5];
    PosInterf = Geometry.IPos(n,:);
    Geometry.DistInterf(n) = sqrt(sum((PosInterf(1:2)-Geometry.Rx_Pos(1:2)).^2));
    Geometry.AOA_Interf(n) = atan2(PosInterf(1,2)-Geometry.Rx_Pos(1,2),PosInterf(1,1)-Geometry.Rx_Pos(1,1))*180/pi;
    Geometry.ZOA_Interf(n) = -atan2(Geometry.Rx_Pos(1,3) - PosInterf(1,3),Geometry.DistInterf(n))*180/pi;
    Geometry.DOA_Interf(n,:) = [Geometry.AOA_Interf(n) Geometry.ZOA_Interf(n)];
end
%% Level D
%% Create the geometry scenario
% Create the layout
l = qd_layout;
%Set scenario parameters
l.set_scenario('QuaDRiGa_UD2D_LOS');
%Define tx and rx
txArr = qd_arrayant('omni');
rxArr = qd_arrayant('omni');
rxARR.no_elements = Pars.NTx_ant^2;
rxArr.element_position = Geometry.BS_Antpos;
%Pssing the Geometry data to the layaout
l.tx_array = txArr;
l.rx_array = rxArr;
l.no_tx = 2 + Pars.NInter;
l.no_rx = 1;
tx_track1 = qd_track('linear', Geometry.TrackV1, pi/2);
tx_track1.name = 'trackV1';
tx_track2 = qd_track('linear', Geometry.TrackV2, pi);
tx_track2.name = 'trackV2';
l.tx_track(1,1)=copy(tx_track1);
l.tx_track(1,2)=copy(tx_track2);
l.tx_position=[Geometry.V1_initPos', Geometry.V2_initPos', Geometry.IPos(:,:)'];
l.rx_position=Geometry.Rx_Pos';
%Visualize the Geometry
l.visualize();
%% Creation Of the Narrow Band Signal
%Signal Narrow Band
Fsin_V1 = 500;          %frequency of sinusoidal signal of vehicle 1 
Fsin_V2 = 1000;         %frequency of sinusoidal signal of vehicle 2
Fsin_Interf = 1500;     %frequency of sinusoidal signal of interferences

Ts = 1e-4;              %Sample time 100 us
Fsample = 1/Ts;         %Sample Frequency of 100 Khz
TsVect = 0:Ts:5/Fsin_V1;%Time vector for the sinusoidals

%Create a Narrow Band Signal
waveform_V1 = sqrt(Pars.txPower*2)*sin(2*pi*Fsin_V1*TsVect);%sinusoidal waveform.
waveform_V2 = sqrt(Pars.txPower*2)*sin(2*pi*Fsin_V2*TsVect);%sinusoidal waveform.
for n = 1:Pars.NInter
    waveform_Interf(n,:) = sqrt(Pars.interfPower*2)*sin(2*pi*Fsin_Interf*TsVect);%sinusoidal waveform.
    waveform_Interf_rx(n,:) = waveform_Interf(n,:)/(4*pi*Geometry.DistInterf(n)/Pars.lambda_Inter(n));
end

%Calculate the received waveform with a Free Space Attenution, at each point
%of the trajectory
for i=1:Geometry.TrackV1+1
    waveform_V1_rx(i,:) = waveform_V1/(4*pi*Geometry.Dist_V1(i)/Pars.lambda);
    waveform_V2_rx(i,:) = waveform_V2/(4*pi*Geometry.Dist_V2(i)/Pars.lambda);
end


%Plot of the recived an transmitted signal, considering just FSPL.
figure;
subplot(221)
plot(TsVect,waveform_V1,'b');
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('V1 Tx Signal'); 
grid on;

subplot(223)
plot(TsVect,waveform_V1_rx(i,:),'b');
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('V1 Rx Signal'); 
grid on;

subplot(222)
plot(TsVect,waveform_V2,'r');
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('V2 Tx Signal');
grid on;

subplot(224)
plot(TsVect,waveform_V2_rx(i,:),'r');
xlabel('Time [s]'); 
ylabel('Amplitude [V]');
title('V2 Rx Signal');
grid on;
%% Power of the Transmitted and Received Signal
%Calculate transmitted and received power (due to FS) of vehicle 1
Power_V1_tx = sum(waveform_V1.^2)/length(waveform_V1);
Power_V1_rx = Power_V1_tx./(4*pi*Geometry.Dist_V1/Pars.lambda).^2;

%Calculate transmitted and received power (due to FS) of vehicle 2
Power_V2_tx = sum(waveform_V2.^2)/length(waveform_V2);
Power_V2_rx = Power_V2_tx./((4*pi*Geometry.Dist_V2/Pars.lambda).^2);

%Calculate transmitted and received power (due to FS) of interferences
Power_Interf_tx = (sum(waveform_Interf.^2,2)./length(waveform_Interf))';
Power_Interf_rx = Power_Interf_tx ./((4*pi*Geometry.DistInterf./Pars.lambda_Inter).^2);

%SINR 
Power_Noise_V1 = Power_V1_rx / db2pow(Pars.SNR); %SNR = Pr/Pn --> Pn = Pr/SNR 
SINR_1 = pow2db(Power_V1_rx ./ (Power_Noise_V1 + sum(Power_Interf_rx) + Power_V2_rx));
Power_Noise_V2 = Power_V2_rx / db2pow(Pars.SNR);
SINR_2 = pow2db(Power_V2_rx ./ (Power_Noise_V2 + sum(Power_Interf_rx) + Power_V1_rx));

%Plot the Received Power at each point of the trajectory.
figure;
subplot(211);
plot(0:1:Geometry.TrackV1,pow2db(Power_V1_rx),'b');
hold on;
plot(0:1:Geometry.TrackV2,pow2db(Power_V2_rx),'r');
xlabel('Distance [m]'); 
ylabel('Received Power [dBW]'); 
title('Received Power'); 
grid on;

%Plot SINR
subplot(212);
plot(0:1:Geometry.TrackV1,SINR_1,'b');
hold on;
plot(0:1:Geometry.TrackV2,SINR_2,'r');
xlabel('Distance [m]'); 
ylabel('SINR [dB]'); 
title('SINR'); 
grid on;
%% Received Signal in the Base Station
for i=1:Geometry.TrackV1+1
    receivedW(:,:,i) = collectPlaneWave(Geometry.BS_array,...
        [waveform_V1_rx(i,:)' waveform_V2_rx(i,:)' waveform_Interf_rx'],...
        [Geometry.AOA_V1(i) Geometry.AOA_V2(i) Geometry.AOA_Interf],...
        Pars.fc_carrier, Pars.c);
    chOut(:,:,i)=awgn(receivedW(:,:,i), Pars.SNR, 'measured');% Add AWGN
end
% Plot the Signal in the 4 first elements of the Base Station (initial position of the vehicles).
figure;
subplot(2,2,1)
plot(TsVect, chOut(:,1,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]');
title('Received Signal, 1st Element'); 
grid on;

subplot(2,2,2)
plot(TsVect, chOut(:,2,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]');
title('Received Signal, 2nd Element'); 
grid on;

subplot(2,2,3)
plot(TsVect, chOut(:,3,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]');
title('Received Signal, 3rd Element'); 
grid on;

subplot(2,2,4)
plot(TsVect, chOut(:,4,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]');
title('Received Signal, 4th Element'); 
grid on;
%% Calculated Steered Vector
%Calculate wave vector
xi_V1 = 2*pi/Pars.lambda;
xi_V2 = 2*pi/Pars.lambda;
%Calculate distance between antennas
d_V1 = Pars.lambda/2;
d_V2 = Pars.lambda/2;
%Assign azimuth angles
phi_V1 = Geometry.AOA_V1;
phi_V2 = Geometry.AOA_V2;
%Calculate the phase for steering vector
alpha_V1 =  - xi_V1 * d_V1 * sind(phi_V1);
alpha_V2 =  - xi_V2 * d_V2 * sind(phi_V2);
%Assign the positional indeces of the antennas
if mod(Pars.NTx_ant,2) == 0
    n = -Pars.NTx_ant/2+1/2:Pars.NTx_ant/2-1/2;
else
    n = -floor(Pars.NTx_ant/2):floor(Pars.NTx_ant/2);
end
rho = 300;      %Amplitude factor
%Steered vector for the first row
Steered_direction_V1 = exp(1j*n.*alpha_V1')';
Steered_direction_V2 = exp(1j*n.*alpha_V2')';
%Steered matrix for all the URA
Pattern_V1 = kron(Steered_direction_V1,ones(Pars.NTx_ant,1));
Pattern_V2 = kron(Steered_direction_V2,ones(Pars.NTx_ant,1));

% Steered Interference, useful for Null Beamforming
% xi_Interf = Pars.lambda_Inter.\2*pi;
% d = Pars.lambda/2;
% phi_Interf = Geometry.AOA_Interf;
% alpha_Interf = - xi_V2 * d_V2 * sind(phi_Interf);
% if mod(Pars.NTx_ant,2) == 0
%     n = -Pars.NTx_ant/2+1/2:Pars.NTx_ant/2-1/2;
% else
%     n = -floor(Pars.NTx_ant/2):floor(Pars.NTx_ant/2);
% end
% Steered_direction_Interf = exp(1j*n.*alpha_Interf')';
%Steered_direction_Interf = kron(Steered_direction_Interf,ones(Param.NTx_ant,1));
%% Null Beamforming
% g_V1 = [ 1 0 zeros(1,Param.NInter)];
% g_V2 = [ 0 1 zeros(1,Param.NInter)];
% for i=1:Geometry.TrackV1+1
% for i=1:1
%     S(:,:,i) = [Steered_direction_V1(:,i) Steered_direction_V2(:,i) Steered_direction_Interf];   
%     
%     inverse = (S(:,:,i)*S(:,:,i)')^(-1);
%     
%     aux1_V1 = g_V1 * S(:,:,i)' * inverse;
%     aux1_V2 = g_V2 * S(:,:,i)' * inverse;
%     
%     aux2_V1 = kron(aux1_V1,ones(Param.NTx_ant,1));
%     aux2_V2 = kron(aux1_V2,ones(Param.NTx_ant,1));
%     
%     w_V1(:,i) = reshape(aux2_V1,[Param.NTx_ant*Param.NTx_ant,1]);
%     w_V2(:,i) = reshape(aux2_V2,[Param.NTx_ant*Param.NTx_ant,1]);
% 
%     nullbeam_signal_V1(:,i) = chOut(:,:,1) * w_V1(:,i);
%     nullbeam_signal_V2(:,i) = chOut(:,:,1) * w_V2(:,i);     
% end
%% Steering Vector Beamformer (Matlab function)
%Calculation of steering vector with the Matlab function
steervector = phased.SteeringVector('SensorArray',Geometry.BS_array);
sv = steervector(Pars.fc_carrier,[Geometry.AOA_V1; zeros(1,Geometry.TrackV1+1)]); 
% figure;
% for i=1:10:Geometry.TrackV1
%     subplot(121)
%     pattern(Geometry.BS_array,Pars.fc_carrier,-180:180,0,'PropagationSpeed',Pars.c,'Weights',Pattern_V1(:,i),...
%         'CoordinateSystem','polar', 'type', 'power'); 
%     title('Steering vector calculated by hand');
%     subplot(122)
%     pattern(Geometry.BS_array,Pars.fc_carrier,-180:180,0,'PropagationSpeed',Pars.c,'Weights',sv(:,i),...
%         'CoordinateSystem','polar', 'type', 'power'); 
%     title('Steering vector calculated by Matlab function');
%     pause(0.01);
% end
%% Beamforming Narrow Band;
%Find the final signal at the receiver after beamforming
for i=1:Geometry.TrackV1+1
    beamformed_signal_V1(i,:) = Pattern_V1(:,i)' * chOut(:,:,i).';
    beamformed_signal_V2(i,:) = Pattern_V2(:,i)' * chOut(:,:,i).';
end
%% Plot Environment Vehicle 1
arrayresp = phased.ArrayResponse('SensorArray',Geometry.BS_array,'WeightsInputPort',true);
figure;
for i=1:1:Geometry.TrackV1+1
%for i=1:1
%%%%%%%%%%%%%% Time Signal %%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,1)
    %plot(TsVect, nullbeam_signal_V1(:,i),'b');
    plot(TsVect, beamformed_signal_V1(i,:),'b');
    ylim([-0.25 0.25]);
    xlabel('Time [s]'); 
    ylabel('Amplitude [V]');
    title('Received Signal Vehicle 1'); 
    grid on;
%%%%%%%%%%%%%%%%% FFT %%%%%%%%%%%%%%%%%%%%%
%     FFT_V1 = fft(nullbeam_signal_V1(:,i));
%     L_V1 = length(nullbeam_signal_V1(:,i));
    %Calculation of FFT
    FFT_V1 = fft(beamformed_signal_V1(i,:));
    L_V1 = length(beamformed_signal_V1(i,:));
    
    P2_V1 = abs(FFT_V1/L_V1);
    P1_V1 = P2_V1(1:L_V1/2+1);
    P1_V1(2:end-1) = 2*P1_V1(2:end-1);

    f = Fsample*(0:(L_V1/2))/L_V1;
    SINR_V1(i) = (P1_V1(floor(Fsin_V1/(Fsample/L_V1)))).^2/((P1_V1(floor(Fsin_Interf/(Fsample/L_V1)))).^2+...
        (P1_V1(floor(Fsin_V2/(Fsample/L_V1)))).^2);
    subplot(2,2,3)
    plot(f,P1_V1,'b')
    ylim([0 0.25]);
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f [Hz]')
    ylabel('|P1(f)| [V]')
    grid on;

%%%%%%%%%%%%%% Polar Environment %%%%%%%%%%%%%%%%%%%%%%

    arrayResponse(:,1) = arrayresp(Pars.fc_carrier,[-180:180],Pattern_V1(:,i));
    %arrayResponse(:,1) = arrayresp(Param.fc_carrier,[-180:180],w_V1(:,i));
    
    %Graphic adjustements
    subplot(2,2,[2 4])
    hNull = polar(pi, rho, 'sw');
    hPlotAxes = get(hNull, 'Parent');
    
    hold(hPlotAxes, 'on')
    
    %plot vehicle 1 position
    h_V1 = polar(deg2rad(Geometry.AOA_V1(i)), Geometry.Dist_V1(i), 'sg');
    set(h_V1, 'MarkerFaceColor', 'b')
    
    %plot vehicle 2 position
    h_V2= polar(deg2rad(Geometry.AOA_V2(i)), Geometry.Dist_V2(i), 'sg');
    set(h_V2, 'MarkerFaceColor', 'r')
    
    %plot interferences position
    h_Inter = polar(deg2rad(Geometry.AOA_Interf), Geometry.DistInterf, 'sg');
    set(h_Inter, 'MarkerFaceColor', 'g')
    
    %plot beamforming
    hTxArray_V1 = polar(deg2rad([-180:180]),(rho/(Pars.NTx_ant^2))*abs(arrayResponse(:,1)'), '-b');
    
    hold(hPlotAxes, 'off')
    title('Radiated Pattern'); 

    pause(0.01);
end 
hold on;
%% Plot Environment Vehicle 2
arrayresp = phased.ArrayResponse('SensorArray',Geometry.BS_array,'WeightsInputPort',true);
figure;
for i=1:Geometry.TrackV2+1
%for i=1:1
%%%%%%%%%%%%%% Time Signal %%%%%%%%%%%%%%%%%%%%%%
    
    subplot(2,2,1)
    plot(TsVect, beamformed_signal_V2(i,:),'r');
    ylim([-0.5 0.5]);
    xlabel('Time [s]'); 
    ylabel('Amplitude [V]');
    title('Received Signal Vehicle 2'); 
    grid on;
    
%%%%%%%%%%%%%%%%% FFT %%%%%%%%%%%%%%%%%%%%%
    Pos_Vehicule_V2 = i;
    FFT_V2 = fft(beamformed_signal_V2(Pos_Vehicule_V2,:));
    L_V2 = length(beamformed_signal_V2(Pos_Vehicule_V2,:));
    P2_V2 = abs(FFT_V2/L_V2);
    P1_V2 = P2_V2(1:L_V2/2+1);
    P1_V2(2:end-1) = 2*P1_V2(2:end-1);

    f = Fsample*(0:(L_V2/2))/L_V2;
    SINR_V2(i) = (P1_V2(floor(Fsin_V2/(Fsample/L_V2)))).^2/((P1_V2(floor(Fsin_Interf/(Fsample/L_V2)))).^2+...
        (P1_V2(floor(Fsin_V1/(Fsample/L_V2)))).^2);
    subplot(2,2,3)
    plot(f,P1_V2,'r') 
    ylim([0 0.25]);
    title('Single-Sided Amplitude Spectrum of X(t)')
    xlabel('f [Hz]')
    ylabel('|P2(f)| [V]')
    grid on;
    
%%%%%%%%%%%%%% Polar Environment %%%%%%%%%%%%%%%%%%%%%%

    arrayResponse(:,2) = arrayresp(Pars.fc_carrier,[-180:180],Pattern_V2(:,i));
    %arrayResponse(:,2) = arrayresp(Param.fc_carrier,[-180:180],w_V2(:,i));
    
    subplot(2,2,[2 4])
    hNull = polar(pi, rho, 'sw');
    hPlotAxes = get(hNull, 'Parent');
    
    hold(hPlotAxes, 'on')
    
    h_V1 = polar(deg2rad(Geometry.AOA_V1(i)), Geometry.Dist_V1(i), 'sg');
    set(h_V1, 'MarkerFaceColor', 'b')
    
    h_V2= polar(deg2rad(Geometry.AOA_V2(i)), Geometry.Dist_V2(i), 'sg');
    set(h_V2, 'MarkerFaceColor', 'r')
    
    h_Inter = polar(deg2rad(Geometry.AOA_Interf), Geometry.DistInterf, 'sg');
    set(h_Inter, 'MarkerFaceColor', 'g')
    
    hTxArray_V2 = polar(deg2rad([-180:180]),(rho/(Pars.NTx_ant^2))*abs(arrayResponse(:,2)'), '-r');
    
    hold(hPlotAxes, 'off')
    title('Radiated Pattern'); 

    pause(0.01);
end 
hold on;
%% SINR
figure,
subplot(211)
plot(1:Geometry.TrackV1+1,pow2db(SINR_V1));
xlabel('Distance [m]'); 
ylabel('SINR [dB]'); 
title('SINR V1'); 
axis tight;
grid on;

subplot(212)
plot(1:Geometry.TrackV2+1,pow2db(SINR_V2),'r');
xlabel('Distance [m]'); 
ylabel('SINR [dB]'); 
title('SINR V2'); 
axis tight;
grid on;
%% Phase Shift Beamformer
%Use of phase shift beamformer function to double check
beamformer = phased.PhaseShiftBeamformer('SensorArray',Geometry.BS_array,...
    'OperatingFrequency',Pars.fc_carrier,'Direction',[Geometry.AOA_V1; zeros(1,Geometry.TrackV1+1)],...
    'WeightsOutputPort',true);
[y,w] = beamformer(chOut(:,:,1)); %Here we apply the beamformer to the signal

figure;
for i=1:5:Geometry.TrackV1
    subplot(311)
    plot(TsVect,waveform_V1_rx(i,:))
    axis tight
    title('Original Signal')
    ylabel('Magnitude')

    subplot(312)
    plot(TsVect,y(:,i))
    axis tight
    title('Received Signal with Phase Shift Beamforming')
    ylabel('Magnitude')
    xlabel('Time [s]')

    subplot(313)
    plot(TsVect, beamformed_signal_V1(i,:));
    axis tight
    title('Received Signal with Steered Vector Beamforming')
    ylabel('Magnitude')
    xlabel('Time [s]')
    
    pause(0.02);
end