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
Pars.fc = 1e9;                      %carrier frequency
Pars.fc_Inter = 1e9;
Pars.c = physconst('Lightspeed');
Pars.lambda = Pars.c/Pars.fc;       %Wavelenght 
Pars.lambda_Inter = Pars.c/Pars.fc_Inter;   %Wavelenght interferences
Pars.SNR=20; %SNR = 20db
Pars.Nb_TxAnt = 5;                  %Number of antennas
Pars.NbInterferences = 2;           %Number of interferences
%% Geometry (Uplink)
Geometry.BSPos=[0,0,25];                    %height of Base Station
Geometry.V1PosStart = [50,-100,1.5];
Geometry.V1PosEnd = [50,0,1.5];
Geometry.V2PosStart = [100,-50,1.5];
Geometry.V2PosEnd = [0,-50,1.5];
Geometry.PosInterf = zeros(Pars.NbInterferences,3);
Geometry.PosInterf(1,:) = [10, -210, 1.5];
Geometry.PosInterf(2,:) = [-150, 100, 1.5];
% Geometry.PosInterf(3,:) = [20, -100, 1.5];

%Trajectory of the vehicles
Geometry.T1=sqrt(sum((Geometry.V1PosEnd(1,1:2)-Geometry.V1PosStart(1,1:2)).^2));
Geometry.T2=sqrt(sum((Geometry.V2PosEnd(1,1:2)-Geometry.V2PosStart(1,1:2)).^2));
%Initial distance Base Station - Vehicles 
Geometry.DistV1Start=sqrt(sum((Geometry.V1PosStart(1,1:2)-Geometry.BSPos(1,1:2)).^2));
Geometry.DistV2Start=sqrt(sum((Geometry.V2PosStart(1,1:2)-Geometry.BSPos(1,1:2)).^2));

%DOA for the two vehicles and interferences for all their positions
Geometry.DOAV1=zeros(abs(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2))+1, 2);
Geometry.DOAV2=zeros(abs(Geometry.V2PosEnd(1)-Geometry.V2PosStart(1))+1, 2);
% DoA for V1
Geometry.DistV1 = zeros(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2),1);
for j=1:abs(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2))+1
    Geometry.V1Pos(j,:)=[Geometry.V1PosStart] + [0 j-1 0];
    CurrentDistV1=sqrt(sum((Geometry.V1Pos(j,1:2)-Geometry.BSPos(1,1:2)).^2));
    Geometry.DistV1(j) = CurrentDistV1;
    Geometry.AOAV1=atan2(Geometry.V1Pos(j,2)-Geometry.BSPos(1,2),Geometry.V1Pos(j,1)-Geometry.BSPos(1,1))*180/pi; %check this
    Geometry.ZOAV1=atan2(Geometry.V1Pos(j,3)-Geometry.BSPos(1,3), CurrentDistV1)*180/pi;
    Geometry.DOAV1(j,:)= [Geometry.AOAV1 Geometry.ZOAV1];   
end
%DoA for V2
Geometry.DistV2 = zeros(abs(Geometry.V2PosEnd(1)-Geometry.V2PosStart(1)),1);
for j = 1:abs(Geometry.V2PosEnd(1)-Geometry.V2PosStart(1))+1
    Geometry.V2Pos(j,:) = [Geometry.V2PosStart] - [(j-1) 0 0];
    CurrentDistV2 = sqrt(sum((Geometry.V2Pos(j,1:2)-Geometry.BSPos(1,1:2)).^2));
    Geometry.DistV2(j) = CurrentDistV2;
    Geometry.AOAV2=atan2(Geometry.V2Pos(j,2)-Geometry.BSPos(1,2), Geometry.V2Pos(j,1)-Geometry.BSPos(1,1))*180/pi;
    Geometry.ZOAV2=atan2(Geometry.V2Pos(j,3)-Geometry.BSPos(1,3), CurrentDistV2)*180/pi;
    Geometry.DOAV2(j,:)= [Geometry.AOAV2 Geometry.ZOAV2];
end
%DOA for all interferences
CurrentPosInterf = zeros(Pars.NbInterferences,3);
for i = 1:Pars.NbInterferences
    CurrentPosInterf = Geometry.PosInterf(i,:);
    Geometry.DistInterf(i) = sqrt(sum((CurrentPosInterf(1:2)-Geometry.BSPos(1,1:2)).^2));
    Geometry.AOA_Interf(i) = atan2(CurrentPosInterf(1,2)-Geometry.BSPos(1,2),CurrentPosInterf(1,1)-Geometry.BSPos(1,1))*180/pi;
    Geometry.ZOA_Interf(i) = atan2(CurrentPosInterf(1,3)-Geometry.BSPos(1,3), Geometry.DistInterf(i))*180/pi;
    Geometry.DOA_Interf(i,:) = [Geometry.AOA_Interf(i) Geometry.ZOA_Interf(i)];
end

Geometry.BSarray = phased.URA('Size',[Pars.Nb_TxAnt Pars.Nb_TxAnt],...    
    'ElementSpacing',[Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');
Geometry.BSAntennaPos = getElementPosition(Geometry.BSarray);
%% Creation of the geometry scenario
%Create the layout
l = qd_layout;
%Set scenario parameters
l.set_scenario('QuaDRiGa_UD2D_LOS');
%Define tx and rx
txArr = qd_arrayant('omni');
rxArr = qd_arrayant('omni');
rxARR.no_elements = Pars.Nb_TxAnt^2;
rxArr.element_position = Geometry.BSAntennaPos;
%Pssing the Geometry data to the layaout
l.tx_array = txArr;
l.rx_array = rxArr;
l.no_tx = 2 + Pars.NbInterferences;
l.no_rx = 1;
tx_track1 = qd_track('linear', Geometry.T1, pi/2);
tx_track1.name = 'trackV1';
tx_track2 = qd_track('linear', Geometry.T2, pi);
tx_track2.name = 'trackV2';
l.tx_track(1,1)=copy(tx_track1);
l.tx_track(1,2)=copy(tx_track2);
l.tx_position=[Geometry.V1PosStart', Geometry.V2PosStart', Geometry.PosInterf(:,:)'];
l.rx_position=Geometry.BSPos';
%Visualize the Geometry
l.visualize();
%% Creation of the OFDM signal
%creation of OFDM signal of vehicle 1
ofdmMod1 = comm.OFDMModulator('FFTLength', 12, ...
    'NumGuardBandCarriers', [0;0], ...
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [2], ...
    'Windowing', false, ...
    'NumSymbols', 14, ...
    'NumTransmitAntennas', 1, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [11]);

M = 4; 	 % Modulation order
ofdmInfo_V1 = info(ofdmMod1);
NumSymbols_V1 = ofdmInfo_V1.DataInputSize(2);
FFTLength_Data_V1 = ofdmInfo_V1.DataInputSize(1)-ofdmInfo_V1.PilotInputSize(1);

% input bit source:
in_V1 = randi([0 1], (NumSymbols_V1*ofdmInfo_V1.DataInputSize(1))*2, 1);

dataInput_V1 = qammod(in_V1, M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
ofdmSize_V1 = ofdmInfo_V1.DataInputSize;
dataInput_V1 = reshape(dataInput_V1, ofdmSize_V1);

% waveform generation:
pilotInput_V1 = ones(ofdmInfo_V1.PilotInputSize(1), NumSymbols_V1);
waveform_OFDM_V1 = ofdmMod1(dataInput_V1, pilotInput_V1);

Fs_V1 = 180000; 				%sample rate of waveform
TsV1 = 1/Fs_V1;                 %sampling time for sinusoid
TsVectV1 = 0:TsV1:TsV1*195;     %Time vector for the sinusoidals
%figure;
%plot(TsVect,waveform_OFDM_V1);
%grid on;
Power_OFDM_V1 = mean(sum(waveform_OFDM_V1.^2));

%creation of OFDM signal of vehicle 2
ofdmMod2 = comm.OFDMModulator('FFTLength', 12, ...
    'NumGuardBandCarriers', [0;0], ...
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [2], ...
    'Windowing', false, ...
    'NumSymbols', 14, ...
    'NumTransmitAntennas', 1, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [2]);

M = 4; 	 % Modulation order
ofdmInfo_V2 = info(ofdmMod1);
NumSymbols_V2 = ofdmInfo_V2.DataInputSize(2);
FFTLength_Data_V2 = ofdmInfo_V2.DataInputSize(1)-ofdmInfo_V2.PilotInputSize(1);

% input bit source:
in_V2 = randi([0 1], (NumSymbols_V2*ofdmInfo_V2.DataInputSize(1))*2, 1);

dataInput_V2 = qammod(in_V2, M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
ofdmSize_V2 = ofdmInfo_V2.DataInputSize;
dataInput_V2 = reshape(dataInput_V2, ofdmSize_V2);

% waveform generation:
pilotInput_V2 = ones(ofdmInfo_V2.PilotInputSize(1), NumSymbols_V2);
waveform_OFDM_V2 = ofdmMod2(dataInput_V2, pilotInput_V2);

Fs_V2 = 180000;                     %sample rate of waveform
TsV2 = 1/Fs_V2;                     %sampling time for sinusoid
TsVectV2 = 0:TsV2:TsV2*195;         %Time vector for the sinusoidals
%figure;
%plot(TsVect,waveform_OFDM_V1);
%grid on;
Power_OFDM_V2 = mean(sum(waveform_OFDM_V2.^2));

%Creation of OFDM signals of interfereces
ofdmModInterf = comm.OFDMModulator('FFTLength', 12, ...
    'NumGuardBandCarriers', [0;0], ...
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [2], ...
    'Windowing', false, ...
    'NumSymbols', 14, ...
    'NumTransmitAntennas', 1, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [6]);

M = 4; 	 % Modulation order
ofdmInfo_Interf = info(ofdmModInterf);
NumSymbols_Interf = ofdmInfo_Interf.DataInputSize(2);
FFTLength_Data_Interf = ofdmInfo_Interf.DataInputSize(1)-ofdmInfo_Interf.PilotInputSize(1);

for i = 1:Pars.NbInterferences
    % input bit source:
    in_Interf = randi([0 1], (NumSymbols_Interf*ofdmInfo_Interf.DataInputSize(1))*2, 1);

    dataInput_Interf = qammod(in_Interf, M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
    ofdmSize_Interf = ofdmInfo_Interf.DataInputSize;
    dataInput_Interf = reshape(dataInput_Interf, ofdmSize_Interf);

    % waveform generation:
    pilotInput_Interf = ones(ofdmInfo_Interf.PilotInputSize(1), NumSymbols_Interf);
    waveform_OFDM_Interf(:,i) = ofdmModInterf(dataInput_Interf, pilotInput_Interf);
end
Fs_Interf = 180000;                         %sample rate of waveform
Ts_Interf = 1/Fs_Interf;                    %sampling time for sinusoid
TsVect_Interf = 0:Ts_Interf:Ts_Interf*195;  %Time vector for the sinusoidals
%figure;
%plot(TsVect,waveform_OFDM_V1);
%grid on;
Power_OFDM_Interf = mean(sum(waveform_OFDM_Interf.^2));
%% Visualize OFDM signals 
% Visualize OFDM Waveformes
% Spectrum Analyzer Vehicule 1 --> need of 1536 signals at least
spectrum = dsp.SpectrumAnalyzer('SampleRate', Fs_V1);
spectrum(waveform_OFDM_V1);
release(spectrum);
% OFDM Subcarrier Mapping Vehicule 1
showResourceMapping(ofdmMod1);

%Spectrum Analyzer Vehicule 2 --> need of 1536 signals at least
spectrum = dsp.SpectrumAnalyzer('SampleRate', Fs_V2);
spectrum(waveform_OFDM_V2);
release(spectrum);
% OFDM Subcarrier Mapping Vehicule 2
showResourceMapping(ofdmMod2);
%% Position assignment (Uplink)
% Vehicle 1
txPos_V1 = Geometry.V1Pos;
htx_V1 = txPos_V1(:,3);

% Vehicle 2
txPos_V2 = Geometry.V2Pos;
htx_V2 = txPos_V2(:,3);

% Interferences
txPos_Interf = Geometry.PosInterf;
htx_Interf = txPos_Interf(:,3);

% Base Station
rxPos = Geometry.BSPos;
hrx = rxPos(3);
rxPosR = rxPos;
rxPosR(3) = - rxPosR(3);
%% Free Space
pos2 = rxPos';
vel1 = [0;0;0];
vel2 = [0;0;0];

%Free Space with the function
channel = phased.FreeSpace('SampleRate',Fs_V1,'OperatingFrequency',Pars.fc);
%sets sample rate with the one fixed on the signal created at the start
%allows the output to mix the rays from the LoS and the reflected beams
figure

% Vehicle 1
pos1 = txPos_V1';
OFDMwav = waveform_OFDM_V1;
%for with function
% for j=1:abs(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2))+1
%     prop_signal_V1(:,j) = channel([swav],pos1(:,j),pos2,vel1,vel2);
%     RxPow_V1(j) = mean(abs(prop_signal_V1(:,j)).^2);
% end
for j=1:abs(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2))+1
    prop_signal_V1(:,j) = OFDMwav ./ (4*pi*Geometry.DistV1(j)/Pars.lambda);
    RxPow_V1(j) = pow2db(mean(abs(prop_signal_V1(:,j)).^2));
end

subplot(211)
semilogy(1:abs(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2))+1,RxPow_V1);
xlabel('Distance [m]'); 
ylabel('Received Power [dBW]'); 
axis tight
title('Received Power Vehicle 1'); 
grid on;

% Vehicle 2
pos1 = txPos_V2';
OFDMwav = waveform_OFDM_V2;
% for j=1:abs(Geometry.V2PosEnd(1)-Geometry.V2PosStart(1))+1
%     prop_signal_V2(:,j) = channel([swav],pos1(:,j),pos2,vel1,vel2);
%     RxPow_V2(j) = mean(abs(prop_signal_V2(:,j)).^2);
% end
for j=1:abs(Geometry.V2PosEnd(1)-Geometry.V2PosStart(1))+1
    prop_signal_V2(:,j) = OFDMwav ./ (4*pi*Geometry.DistV2(j)/Pars.lambda);
    RxPow_V2(j) = pow2db(mean(abs(prop_signal_V2(:,j)).^2));
end
subplot(212)
semilogy(1:abs(Geometry.V2PosEnd(1)-Geometry.V2PosStart(1))+1,RxPow_V2,'r');
xlabel('Distance [m]'); 
ylabel('Received Power [dBW]');
axis tight
title('Received Power Vehicle 2'); 
grid on;

% Interferences
pos1 = txPos_Interf';
OFDMwav = waveform_OFDM_Interf;
% for i = 1:Pars.NbInterferences
%     prop_signal_Interf(:,i) = channel([swav(:,i)],pos1(:,i),pos2,vel1,vel2);
%     RxPow_Interf(i) = mean(abs(prop_signal_Interf(:,i)).^2);
% end
for i = 1:Pars.NbInterferences
    prop_signal_Interf(:,i) = OFDMwav(:,i) ./ (4*pi*Geometry.DistInterf(i)/Pars.lambda_Inter);
    RxPow_Interf(i) = pow2db(mean(abs(prop_signal_Interf(:,i)).^2));
end
%% Recived Signal in the Base Station
for i=1:abs(Geometry.V1PosEnd(2)-Geometry.V1PosStart(2))+1
    receivedW(:,:,i) = collectPlaneWave(Geometry.BSarray,...
        [prop_signal_V1(:,i) prop_signal_V2(:,i) prop_signal_Interf],...
        [Geometry.DOAV1(i,:)' Geometry.DOAV2(i,:)' Geometry.DOA_Interf'],...
        Pars.fc);
    chOut(:,:,i)=awgn(receivedW(:,:,i), Pars.SNR, 'measured');% Add AWGN
end
% Plot the Signal in the 4 first elements of the Base Station array (initial position of the vehicules).
figure;
subplot(2,2,1)
plot(TsVectV1, chOut(:,1,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('Received Signal, 1st Element'); 
grid on;

subplot(2,2,2)
plot(TsVectV1, chOut(:,2,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('Received Signal, 2nd Element'); 
grid on;

subplot(2,2,3)
plot(TsVectV1, chOut(:,3,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('Received Signal, 3rd Element'); 
grid on;

subplot(2,2,4)
plot(TsVectV1, chOut(:,4,1));
xlabel('Time [s]'); 
ylabel('Amplitude [V]'); 
title('Received Signal, 4th Element'); 
grid on;
%% Calculated Steered Vector
%for Vehicle 1 and Vehicle 2
%Calculate wave vector
xi_V1 = 2*pi/Pars.lambda;
xi_V2 = 2*pi/Pars.lambda;
%Calculate distance between antennas
d_V1 = Pars.lambda/2;
d_V2 = Pars.lambda/2;
%Assign azimuth angles
phi_V1 = Geometry.DOAV1(:,1);
phi_V2 = Geometry.DOAV2(:,1);
%Assign elevation angles
teta_V1 = Geometry.DOAV1(:,2);
teta_V2 = Geometry.DOAV1(:,2);
%Calculate the phase for steering vector for vehicle 1
alpha_V1_azimuth = xi_V1 * d_V1 .* sind(phi_V1) .* cosd(teta_V1);
alpha_V1_elevation = xi_V1 * d_V1 .* sind(teta_V1) .* cosd(phi_V1);
%Calculate the phase for steering vector for vehicle 2
alpha_V2_azimuth = xi_V2 * d_V2 .* sind(phi_V2) .* cosd(teta_V2);
alpha_V2_elevation = xi_V2 * d_V2 .* sind(teta_V2) .* cosd(phi_V2);
%Assign the positional indeces of the antennas
if mod(Pars.Nb_TxAnt,2) == 0
    n = -Pars.Nb_TxAnt/2+1/2:Pars.Nb_TxAnt/2-1/2;
else
    n = -floor(Pars.Nb_TxAnt/2):floor(Pars.Nb_TxAnt/2);
end
rho=300;            %Amplitude factor
%Steered_direction_V1 = (rho/(Pars.NTx_ant^2))*exp(1j*n.*alpha_V1')';
%Steered_direction_V2 = (rho/(Pars.NTx_ant^2))*exp(1j*n.*alpha_V2')';

%Steered vector for the first row of antennas for vehicle 1
Steered_direction_V1_azimuth = exp(1j* n.* alpha_V1_azimuth);
Steered_direction_V1_elevation = exp(1j* -n.* alpha_V1_elevation);
%Steered vector for the first row of antennas for vehicle 2
Steered_direction_V2_azimuth = exp(1j* n.* alpha_V2_azimuth);
Steered_direction_V2_elevation = exp(1j* -n.* alpha_V2_elevation);

%initialization
Steer_Kron_V1 = zeros (Pars.Nb_TxAnt*Pars.Nb_TxAnt,Geometry.T1+1);
Steer_Kron_V2 = zeros (Pars.Nb_TxAnt*Pars.Nb_TxAnt,Geometry.T2+1);
% Steer_Matrix_V1 = zeros (Pars.Nb_TxAnt,Pars.Nb_TxAnt,Geometry.T1+1);
% Steer_Matrix_V2 = zeros (Pars.Nb_TxAnt,Pars.Nb_TxAnt,Geometry.T2+1);

for i=1:Geometry.T1+1
    %creation of steering matrix
    Steer_Kron_V1(:,i) = kron(Steered_direction_V1_azimuth(i,:),Steered_direction_V1_elevation(i,:));
    Steer_Kron_V2(:,i) = kron(Steered_direction_V2_azimuth(i,:),Steered_direction_V2_elevation(i,:));
    
%     Steer_Matrix_V1(:,:,i) = reshape(Steer_Kron_V1(:,i),[Pars.Nb_TxAnt Pars.Nb_TxAnt]);
%     Steer_Matrix_V2(:,:,i) = reshape(Steer_Kron_V2(:,i),[Pars.Nb_TxAnt Pars.Nb_TxAnt]);
end
%% Steering Vector (Matlab Function)
steervector = phased.SteeringVector('SensorArray',Geometry.BSarray);
sv = steervector(Pars.fc,Geometry.DOAV1'); 
% sv_matrix=reshape(sv,[Pars.Nb_TxAnt Pars.Nb_TxAnt Geometry.T1+1]);
%Comparisons between calculated weights by hand and calculated weights with
%Matlab function
figure,
for i=1:10:Geometry.T1
%for i=1:1
    subplot(221)
    pattern(Geometry.BSarray,Pars.fc,0,-90:90,'PropagationSpeed',Pars.c,'Weights',sv(:,i),...
        'CoordinateSystem','polar', 'type', 'power');    
    subplot(222)
    pattern(Geometry.BSarray,Pars.fc,0,-90:90,'PropagationSpeed',Pars.c,'Weights',Steer_Kron_V1(:,i),...
        'CoordinateSystem','polar', 'type', 'power'); 
    subplot(223)
    pattern(Geometry.BSarray,Pars.fc,-180:180,0,'PropagationSpeed',Pars.c,'Weights',sv(:,i),...
        'CoordinateSystem','polar', 'type', 'power'); 
    subplot(224)
    pattern(Geometry.BSarray,Pars.fc,-180:180,0,'PropagationSpeed',Pars.c,'Weights',Steer_Kron_V1(:,i),...
        'CoordinateSystem','polar', 'type', 'power');
    pause(0.1);
end
%% Demodulator
ofdmDemod1 = comm.OFDMDemodulator(ofdmMod1);
ofdmDemod2 = comm.OFDMDemodulator(ofdmMod2);

out1 = ofdmDemod1(chOut(:,1,1));
out2 = ofdmDemod2(chOut(:,1,1));
%Scatter plot of the first antenna in the first position of the trajectory
figure,
x=real(out1);
x=reshape(x,[ofdmInfo_V1.DataInputSize(1)*ofdmInfo_V1.DataInputSize(2),1]);
y=imag(out1);
y=reshape(y,[ofdmInfo_V1.DataInputSize(1)*ofdmInfo_V1.DataInputSize(2),1]);
scatter(x,y)
title('Vehicle 1 without Beamforming')
grid on;

figure,
x=real(out2);
x=reshape(x,[ofdmInfo_V2.DataInputSize(1)*ofdmInfo_V2.DataInputSize(2),1]);
y=imag(out2);
y=reshape(y,[ofdmInfo_V2.DataInputSize(1)*ofdmInfo_V2.DataInputSize(2),1]);
scatter(x,y, 'r')
title('Vehicle 2 without Beamforming')
grid on;
%% Phased shifted Beamformer (Matlab Function)
beamformer_V1 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'Direction',Geometry.DOAV1(100,:)',...
    'WeightsOutputPort',true);

beamformer_V2 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'Direction',Geometry.DOAV2(100,:)',...
    'WeightsOutputPort',true);
%% Plot Environment Vehicule 1
arrayresp = phased.ArrayResponse('SensorArray',Geometry.BSarray,'WeightsInputPort',true);
figure,
for i=1:Geometry.T1+1
    beamformed_signal_V1 = chOut(:,:,i) * conj(Steer_Kron_V1(:,i));
    %[ps1,w1] = beamformer_V1(chOut(:,:,i)); %Here we apply the beamformr to the signal
    out1_Beam(:,:,i) = ofdmDemod1(beamformed_signal_V1);    %output signal of beamformed signal
    
    subplot(211)
    x=real(out1_Beam(:,:,i));
    x=reshape(x,[ofdmInfo_V1.DataInputSize(1)*ofdmInfo_V1.DataInputSize(2),1]);
    y=imag(out1_Beam(:,:,i));
    y=reshape(y,[ofdmInfo_V1.DataInputSize(1)*ofdmInfo_V1.DataInputSize(2),1]);
    scatter(x,y)
    title('Vehicle 1 with Beamforming')
    grid on;
   
    %%%%%%%%%%%%%% Polar Environment %%%%%%%%%%%%%%%%%%%%%%

    arrayResponse(:,1) = arrayresp(Pars.fc,[-180:180],Steer_Kron_V1(:,i));
    
    subplot(212)
    hNull = polar(pi, rho, 'sw');
    hPlotAxes = get(hNull, 'Parent');
    
    hold(hPlotAxes, 'on')
    
    h_V1 = polar(deg2rad(Geometry.DOAV1(i,1)), Geometry.DistV1(i), 'sg');
    set(h_V1, 'MarkerFaceColor', 'b')
    
    h_V2= polar(deg2rad(Geometry.DOAV2(i,1)), Geometry.DistV2(i), 'sg');
    set(h_V2, 'MarkerFaceColor', 'r')
    
    h_Inter = polar(deg2rad(Geometry.DOA_Interf(:,1).'), Geometry.DistInterf, 'sg');
    set(h_Inter, 'MarkerFaceColor', 'g')
    
    hTxArray_V1 = polar(deg2rad([-180:180]),((rho/(Pars.Nb_TxAnt^2))*abs(arrayResponse(:,1)')), '-b');
    
    hold(hPlotAxes, 'off')
    title('Radiated Pattern'); 

    pause(0.1);
end
%% Plot Environment Vehicule 2
figure,
for i=1:Geometry.T2+1
    
    beamformed_signal_V2 = chOut(:,:,i) * conj(Steer_Kron_V2(:,i));
    %[ps2,w2] = beamformer_V2(chOut(:,:,i)); %Here we apply the beamformr to the signal
    out2_Beam(:,:,i) = ofdmDemod2(beamformed_signal_V2);
    
    subplot(211)
    x=real(out2_Beam(:,:,i));
    x=reshape(x,[ofdmInfo_V2.DataInputSize(1)*ofdmInfo_V2.DataInputSize(2),1]);
    y=imag(out2_Beam(:,:,i));
    y=reshape(y,[ofdmInfo_V2.DataInputSize(1)*ofdmInfo_V2.DataInputSize(2),1]);
    scatter(x,y,'r')
    title('Vehicle 2 with Beamforming')
    grid on;

    %%%%%%%%%%%%%% Polar Environment %%%%%%%%%%%%%%%%%%%%%%

    arrayResponse(:,2) = arrayresp(Pars.fc,[-180:180],Steer_Kron_V2(:,i));
    
    subplot(212)
    hNull = polar(pi, rho, 'sw');
    hPlotAxes = get(hNull, 'Parent');
    
    hold(hPlotAxes, 'on')
    
    h_V1 = polar(deg2rad(Geometry.DOAV1(i,1)), Geometry.DistV1(i), 'sg');
    set(h_V1, 'MarkerFaceColor', 'b')
    
    h_V2= polar(deg2rad(Geometry.DOAV2(i,1)), Geometry.DistV2(i), 'sg');
    set(h_V2, 'MarkerFaceColor', 'r')
    
    h_Inter = polar(deg2rad(Geometry.DOA_Interf(:,1).'), Geometry.DistInterf, 'sg');
    set(h_Inter, 'MarkerFaceColor', 'g')
    
    hTxArray_V2 = polar(deg2rad([-180:180]),((rho/(Pars.Nb_TxAnt^2))*abs(arrayResponse(:,2)')), '-r');
    
    hold(hPlotAxes, 'off')
    title('Radiated Pattern'); 

    pause(0.1);
end
%% BER
for i=1:Geometry.T1+1
    out1 = ofdmDemod1(chOut(:,1,i));
    %BER per V1
    bit_rx_V1 = qamdemod(out1, M, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true); %QAM demod
    out_bit_rx_V1 = reshape(bit_rx_V1,[length(in_V1),1]);
    nrBit_Rx_V1_wrong = sum(out_bit_rx_V1 ~= in_V1); %number of wrong bits
    nrBit_Tx_V1 = length(in_V1);
    BER_V1(i)=nrBit_Rx_V1_wrong/nrBit_Tx_V1; %wrong bits over sent bit without beamforming

    bit_rx_V1 = qamdemod(out1_Beam(:,:,i), M, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true); %QAM demod
    out_bit_rx_V1 = reshape(bit_rx_V1,[length(in_V1),1]); 
    nrBit_Rx_V1_wrong = sum(out_bit_rx_V1 ~= in_V1); %number of wrong bits
    nrBit_Tx_V1 = length(in_V1);
    BER_V1_Beam(i)=nrBit_Rx_V1_wrong/nrBit_Tx_V1; %wrong bits over sent bit with beamforming
end
for i=1:Geometry.T2+1
    out2 = ofdmDemod2(chOut(:,1,i));
    %BER per V2
    bit_rx_V2 = qamdemod(out2, M, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true); %QAM demod
    out_bit_rx_V2 = reshape(bit_rx_V2,[length(in_V2),1]);
    nrBit_Rx_V2_wrong = sum(out_bit_rx_V2 ~= in_V2); %number of wrong bits
    nrBit_Tx_V2 = length(in_V2);
    BER_V2(i)=nrBit_Rx_V2_wrong/nrBit_Tx_V2; %wrong bits over sent bit without beamforming

    bit_rx_V2 = qamdemod(out2_Beam(:,:,i), M, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true); %QAM demod
    out_bit_rx_V2 = reshape(bit_rx_V2,[length(in_V2),1]);
    nrBit_Rx_V2_wrong = sum(out_bit_rx_V2 ~= in_V2); %number of wrong bits
    nrBit_Tx_V2 = length(in_V2);
    BER_V2_Beam(i)=nrBit_Rx_V2_wrong/nrBit_Tx_V2; %wrong bits over sent bit with beamforming
end

%Plot of the BER
figure,
subplot(211)
plot(1:Geometry.T1+1,BER_V1,'b');
xlabel('Distance [m]'); 
ylabel('BER'); 
axis ([0 Geometry.T1+1 0 1])
title('Vehicle 1 without Beamforming');
grid on

subplot(212)
plot(1:Geometry.T1+1,BER_V1_Beam,'b');
xlabel('Distance [m]'); 
ylabel('BER'); 
axis ([0 Geometry.T1+1 0 1])
title('Vehicle 1 with Beamforming');
grid on

figure,
subplot(211)
plot(1:Geometry.T2+1,BER_V2, 'r');
xlabel('Distance [m]'); 
ylabel('BER'); 
axis ([0 Geometry.T2+1 0 1])
title('Vehicle 2 without Beamforming')
grid on

subplot(212)
plot(1:Geometry.T2+1,BER_V2_Beam, 'r');
xlabel('Distance [m]'); 
ylabel('BER'); 
axis ([0 Geometry.T2+1 0 1])
title('Vehicle 2 with Beamforming')
grid on