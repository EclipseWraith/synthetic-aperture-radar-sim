% =========================================================================
% ADVANCED STRIPMAP SAR - NAVAL STRIKE GROUP (OPTIMIZED FOR SPEED)
% Description: Presentation-ready script utilizing a 3D Aircraft Carrier 
%              and multiple Fighter Jets. Optimized for live-demo speed.
% =========================================================================
clear; clc; close all;

%% 1. Radar Hardware & Physics Parameters
speedOfLight = physconst('LightSpeed');
carrierFreq  = 4e9;                       % 4 GHz (C-Band)
lambda       = speedOfLight / carrierFreq;

% Resolution Requirements
rangeRes      = 3;                        % Target: 3m Range Resolution
bandwidth     = speedOfLight / (2 * rangeRes);

% Timing & Sampling
pulseRepFreq = 1000;                     
pulseWidth   = 3e-6;                      
samplingRate = 120e6;                     
antennaAperture = 4;                      % 4-meter physical antenna

% Kinematics (OPTIMIZED: Reduced flight time from 4s to 2.5s for speed)
platformSpeed  = 100;  
flightDuration = 4; 
slowTimeStep   = 1 / pulseRepFreq;
totalPulses    = (flightDuration / slowTimeStep) + 1;

% Range Limits 
maxSlantRange = 2500;
truncSamples  = ceil((2 * maxSlantRange / speedOfLight) * samplingRate);
fastTimeAxis  = (0:1/samplingRate:(truncSamples-1)/samplingRate);
referenceRng  = 1000;                     

%% 2. Command Window Presentation Output
beamwidth_rad = lambda / antennaAperture;
normalRadarAzRes = referenceRng * beamwidth_rad;
sarAzRes = antennaAperture / 2;
fprintf('====================================================\n');
fprintf('           NAVAL SCENE PHYSICS LIMITS               \n');
fprintf('====================================================\n');
fprintf('Frequency: %.2f GHz (Wavelength: %.4f m)\n', carrierFreq/1e9, lambda);
fprintf('Theoretical Slant Range Resolution (c/2B): %.2f m\n', rangeRes);
fprintf('SAR (Synthetic Aperture) Azimuth Resolution: %.2f meters\n', sarAzRes);
fprintf('====================================================\n\n');

%% 3. Phased Array Toolbox Initialization
fmWaveform = phased.LinearFMWaveform('SampleRate', samplingRate, ...
    'PulseWidth', pulseWidth, 'PRF', pulseRepFreq, 'SweepBandwidth', bandwidth);
platformRadar = phased.Platform('InitialPosition', [0; -125; 500], ... % Adjusted start for 2.5s flight
    'Velocity', [0; platformSpeed; 0]); 

antenna     = phased.CosineAntennaElement('FrequencyRange', [1e9 6e9]);
antennaGain = aperture2gain(antennaAperture, lambda); 
transmitter = phased.Transmitter('PeakPower', 50e3, 'Gain', antennaGain);
radiator    = phased.Radiator('Sensor', antenna, 'OperatingFrequency', carrierFreq, 'PropagationSpeed', speedOfLight);
collector   = phased.Collector('Sensor', antenna, 'PropagationSpeed', speedOfLight, 'OperatingFrequency', carrierFreq);
receiver    = phased.ReceiverPreamp('SampleRate', samplingRate, 'NoiseFigure', 30); 
propChannel = phased.FreeSpace('PropagationSpeed', speedOfLight, ...
    'OperatingFrequency', carrierFreq, 'SampleRate', samplingRate, 'TwoWayPropagation', true);

%% 4. Target Scene Geometry (Naval Strike Group & Formations)
carrierDeckHeight = 15; % Flight deck is 15 meters above the water

% 1. The Aircraft Carrier (Center: 1000m Range, 0m Azimuth)
[c1_x, c1_y, c1_z, c1_rcs] = createAircraftCarrier3D(1000, 0, 0); 

% 2. DECK: Two Jets Parked Side-by-Side (Z = 15m)
[j1_x, j1_y, j1_z, j1_rcs] = createFighterJet3D(1040, 10, carrierDeckHeight); 
[j2_x, j2_y, j2_z, j2_rcs] = createFighterJet3D(1040, -10, carrierDeckHeight); 

% 3. AIR: Jet flying directly ABOVE the carrier (Z = 150m Altitude)
[j3_x, j3_y, j3_z, j3_rcs] = createFighterJet3D(1000, 0, 150); 

% 4. AIR: Escort V-Formation flying alongside (Z = 100m Altitude)
[j4_x, j4_y, j4_z, j4_rcs] = createFighterJet3D(1000, 70, 100);  % Lead Jet
[j5_x, j5_y, j5_z, j5_rcs] = createFighterJet3D(1020, 85, 100);  % Wingman Right
[j6_x, j6_y, j6_z, j6_rcs] = createFighterJet3D(1020, 55, 100);  % Wingman Left

% >>> ADD NEW OBJECTS TO THESE MATRICES TO AUTO-UPDATE GROUND TRUTH <<<
targetCoords = [c1_x, j1_x, j2_x, j3_x, j4_x, j5_x, j6_x; ... % Append new X arrays here
                c1_y, j1_y, j2_y, j3_y, j4_y, j5_y, j6_y; ... % Append new Y arrays here
                c1_z, j1_z, j2_z, j3_z, j4_z, j5_z, j6_z];    % Append new Z arrays here
targetRCS    = [c1_rcs, j1_rcs, j2_rcs, j3_rcs, j4_rcs, j5_rcs, j6_rcs]; % Append new RCS arrays here
targetVels   = zeros(3, length(targetRCS));

% Swerling 1 Scintillation Applied
targetScatterers = phased.RadarTarget('OperatingFrequency', carrierFreq, ...
    'MeanRCS', targetRCS, 'Model', 'Swerling1');
platformTargets  = phased.Platform('InitialPosition', targetCoords, 'Velocity', targetVels);

%% 5. Forward Model: SAR Data Collection
radarBoresight = zeros(1, size(targetCoords, 2));
rawDataMatrix  = complex(zeros(truncSamples, totalPulses)); 

disp('Simulating Complex Naval Geometry (Collecting Pulse Data)...');
for i = 1:totalPulses
    [radarPos, radarVel] = platformRadar(slowTimeStep);
    [tgtPos, tgtVel]     = platformTargets(slowTimeStep);
    [tgtRange, tgtAngle] = rangeangle(tgtPos, radarPos);
    
    sigTx = fmWaveform();
    sigTx = sigTx(1:truncSamples);
    sigTx = transmitter(sigTx);
    
    tgtAngle(1,:) = radarBoresight; 
    
    sigProp = radiator(sigTx, tgtAngle);
    sigProp = propChannel(sigProp, radarPos, tgtPos, radarVel, tgtVel);
    sigRefl = targetScatterers(sigProp, true); 
    
    sigRx = collector(sigRefl, tgtAngle);
    rawDataMatrix(:, i) = receiver(sigRx);
end

disp('Injecting Multiplicative Speckle and Thermal Noise into Phase History...');
speckle_variance = 0.05;
speckle_noise = 1 + sqrt(speckle_variance/2) * (randn(size(rawDataMatrix)) + 1i*randn(size(rawDataMatrix)));
rawDataMatrix = rawDataMatrix .* speckle_noise; 

thermal_noise = 0.01 * max(abs(rawDataMatrix(:))) * (randn(size(rawDataMatrix)) + 1i*randn(size(rawDataMatrix)));
rawDataMatrix = rawDataMatrix + thermal_noise;
disp('Data Collection Complete.');

%% 6. DSP Pipeline: Range & Azimuth Compression (w-k Algorithm)
disp('Executing Range Compression...');
pulseCompression = phased.RangeResponse('RangeMethod', 'Matched filter', ...
    'PropagationSpeed', speedOfLight, 'SampleRate', samplingRate);
matchingCoeffs   = getMatchedFilter(fmWaveform);
[rangeCompData, rangeGrid] = pulseCompression(rawDataMatrix, matchingCoeffs);

disp('Executing Range Migration Algorithm (Stolt Interpolation)...');
rmaImage = applyRMA(rangeCompData, fastTimeAxis, carrierFreq, samplingRate, ...
    pulseRepFreq, platformSpeed, totalPulses, speedOfLight, referenceRng);

%% 7. VISUALIZATION 1: MACRO FLIGHT PATH & MICRO GROUND TRUTH
disp('Generating Dark Mode 3D Ground Truth Visualization...');
figure('Name', 'Naval Strike Group (Macro & Micro)', 'Position', [100, 200, 1200, 500], 'Color', 'k');

% Colors: Cyan for sea-level/deck, Neon Yellow for high-altitude/island
colors = repmat([0 1 1], length(targetCoords(3,:)), 1); 
colors(targetCoords(3,:) > 15, :) = repmat([1 1 0], sum(targetCoords(3,:) > 15), 1); 

% ==========================================
% LEFT PANEL: MACRO 3D VIEW (FLIGHT PATH)
% ==========================================
subplot(1, 2, 1);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); 
hold on; grid on; view([-25, 25]); 

scatter3(targetCoords(1,:), targetCoords(2,:), targetCoords(3,:), 10, colors, 'filled');
flight_y = linspace(-125, 125, 100);
plot3(zeros(1, 100), flight_y, 500 * ones(1, 100), 'w-', 'LineWidth', 2); 
scatter3(0, -125, 500, 80, 'w', '^', 'filled');
scatter3(0, 0, 500, 120, 'w', 'square', 'filled');

plot3([0, 1000], [0, 0], [500, 0], 'r-', 'LineWidth', 1.5, 'Color', [1 0 0 0.6]);
text(0, 0, 550, 'SAR Platform (500m)', 'Color', 'w', 'FontWeight', 'bold');

xlabel('Ground Range (m)', 'FontWeight', 'bold');
ylabel('Azimuth (m)', 'FontWeight', 'bold');
zlabel('Altitude (m)', 'FontWeight', 'bold');
title('Macro View: Radar Geometry', 'FontSize', 14, 'Color', 'w');
xlim([-100 1100]); ylim([-250 250]); zlim([0 600]); 
rotate3d on; 
hold off;

% ==========================================
% RIGHT PANEL: MICRO TOP-DOWN VIEW
% ==========================================
subplot(1, 2, 2);
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w'); 
hold on; grid on; view(2); 

scatter3(targetCoords(1,:), targetCoords(2,:), targetCoords(3,:), 15, colors, 'filled');

text(960, -25, 15, 'Carrier Island', 'FontSize', 9, 'Color', 'w');
text(1040, 0, 15, 'Deck Jets', 'FontSize', 9, 'Color', 'w');
text(1000, 0, 150, 'Jet (Above Carrier)', 'FontSize', 9, 'Color', 'y', 'FontWeight', 'bold');
text(980, 70, 100, 'Escort Formation', 'FontSize', 9, 'Color', 'y', 'FontWeight', 'bold');

xlabel('Ground Range (m)', 'FontWeight', 'bold');
ylabel('Azimuth (m)', 'FontWeight', 'bold');
title('Micro View: Target Positions', 'FontSize', 14, 'Color', 'w');
xlim([900 1100]); ylim([-80 120]); 
axis equal; 
hold off;

%% 8. VISUALIZATION 2: SAR IMAGE
disp('Generating Photorealistic SAR Image Results...');
figure('Name', 'SAR Reconnaissance Image', 'Position', [150, 50, 700, 700], 'Color', 'k');
azimuthGrid = linspace(-125, 125, totalPulses); 
azCrop = 250:2250; 
rgCrop = 700:1500;

raw_complex_crop = rmaImage(azCrop, rgCrop).';
sar_dB = 20 * log10(abs(raw_complex_crop));

imagesc(azimuthGrid(azCrop), rangeGrid(rgCrop), sar_dB);
set(gca, 'YDir', 'normal'); 
colormap(gca, 'gray'); 
clim([max(sar_dB(:)) - 35, max(sar_dB(:))]); 
axis image; 

% Expanded view to fit the Layover effect and Escort Formation
xlim([-50 110]);     % Horizontal Axis (Azimuth)
ylim([1030 1170]);   % Vertical Axis (Slant Range)

title('SAR Image (Native Slant Range)', 'Color', 'w', 'FontSize', 14);
xlabel('Azimuth / Cross-Range (Meters)', 'Color', 'w', 'FontWeight', 'bold');
ylabel('Slant Range (Meters)', 'Color', 'w', 'FontWeight', 'bold');

%% =========================================================================
% FUNCTION DEFINITIONS
% =========================================================================

% --- True 3D Mesh Generator for Aircraft Carrier ---
function [X, Y, Z, RCS] = createAircraftCarrier3D(centerX, centerY, baseZ)
    % Carrier Dimensions: Length 150m (X), Width 40m (Y), Height 15m (Z)
    
    % OPTIMIZED: Reduced point density to speed up the loop (400 points total)
    x_d = linspace(centerX-75, centerX+75, 40);
    y_d = linspace(centerY-20, centerY+20, 10);
    [X_D, Y_D] = meshgrid(x_d, y_d);
    
    x_deck = X_D(:).';
    y_deck = Y_D(:).';
    z_deck = repmat(baseZ + 15, 1, length(x_deck));
    rcs_deck = 5 * ones(1, length(x_deck)); 
    
    % Island Superstructure
    x_i = linspace(centerX-15, centerX+5, 8);
    y_i = linspace(centerY-10, centerY-20, 4);
    z_i = linspace(baseZ+15, baseZ+30, 4);
    [X_I, Y_I, Z_I] = ndgrid(x_i, y_i, z_i);
    
    x_isl = X_I(:).';
    y_isl = Y_I(:).';
    z_isl = Z_I(:).';
    rcs_isl = 50 * ones(1, length(x_isl)); 
    
    X = [x_deck, x_isl]; Y = [y_deck, y_isl]; Z = [z_deck, z_isl]; RCS = [rcs_deck, rcs_isl];
end

% --- True 3D Mesh Generator for Fighter Jet ---
function [X, Y, Z, RCS] = createFighterJet3D(centerX, centerY, baseZ)
    % Jet Dimensions: Length 15m (X), Wingspan 10m (Y), Height 3m (Z)
    
    x_fuse = linspace(centerX-7.5, centerX+7.5, 12);
    y_fuse = zeros(1, 12) + centerY;
    z_fuse = zeros(1, 12) + baseZ;
    rcs_fuse = 10 * ones(1, 12);
    
    x_canopy = linspace(centerX+1, centerX+3, 4);
    y_canopy = zeros(1, 4) + centerY;
    z_canopy = linspace(baseZ+0.5, baseZ+1.5, 4);
    rcs_canopy = 15 * ones(1, 4);
    
    x_wing = linspace(centerX+1, centerX-4, 8);
    y_wing_port = linspace(centerY+0.5, centerY+5, 8);
    y_wing_stbd = linspace(centerY-0.5, centerY-5, 8);
    z_wing = zeros(1, 8) + baseZ;
    rcs_wing = 20 * ones(1, 8); 
    
    x_htail = linspace(centerX-6, centerX-7.5, 4);
    y_htail_port = linspace(centerY+0.5, centerY+2.5, 4);
    y_htail_stbd = linspace(centerY-0.5, centerY-2.5, 4);
    z_htail = zeros(1, 4) + baseZ;
    rcs_htail = 5 * ones(1, 4);
    
    x_vtail = linspace(centerX-6, centerX-7.5, 4);
    y_vtail = zeros(1, 4) + centerY;
    z_vtail = linspace(baseZ, baseZ+3, 4); 
    rcs_vtail = 50 * ones(1, 4); 
    
    X = [x_fuse, x_canopy, x_wing, x_wing, x_htail, x_htail, x_vtail];
    Y = [y_fuse, y_canopy, y_wing_port, y_wing_stbd, y_htail_port, y_htail_stbd, y_vtail];
    Z = [z_fuse, z_canopy, z_wing, z_wing, z_htail, z_htail, z_vtail];
    RCS = [rcs_fuse, rcs_canopy, rcs_wing, rcs_wing, rcs_htail, rcs_htail, rcs_vtail];
end

function azCompressedData = applyRMA(sigData, fastTime, fc, fs, prf, speed, numPulses, c, Rc)
    freqRange = linspace(fc - fs/2, fc + fs/2, length(fastTime));
    kRange    = 2 * (2 * pi * freqRange) / c; 
    kAz       = 2 * pi * linspace(-prf/2, prf/2, numPulses) ./ speed; 
    
    kAzimuth = kAz.';
    kx = kRange.^2 - kAzimuth.^2;
    kx = sqrt(kx .* (kx > 0)); 
    
    kFinal = exp(1i * kx .* Rc);
    sData  = fftshift(fft(fftshift(fft(sigData, [], 1), 1), [], 2), 2);
    fsmPol = (sData.') .* kFinal;
    
    stoltPol = fsmPol;
    for i = 1:size((fsmPol), 1)
        stoltPol(i, :) = interp1(kx(i, :), fsmPol(i, :), kRange(1, :));
    end
    
    stoltPol(isnan(stoltPol)) = 1e-30;
    stoltPol = stoltPol .* exp(-1i * kRange .* Rc);
    azCompressedData = ifft2(stoltPol);
end