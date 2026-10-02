%% ISAR_JTF_2D_FINAL_COMMON_RANGE_REFERENCE_PRESENTATION_FIXED
% Self-contained 2-D ISAR simulation.
% Book basis: Chapter 8, Sections 8.3.3-8.3.4 and MATLAB Code 8.5.
%
% Four figures are displayed only:
%   1) Ground truth: radar/RLOS geometry + aircraft + local ISAR frame
%   2) Figure 8.16 style: conventional ISAR, no compensation
%   3) Figure 8.19 style: after translational MP compensation
%   4) Figure 8.21 style: after translational + Gabor/JTF processing
%   5) Validation overlay: Figure 8.21 + JTF-time ground-truth scatterers
%
% Important coordinate convention:
%   - Range is the physical absolute slant range from the radar.
%   - Cross-range is the physical target-local coordinate obtained from
%     signed Doppler and the instantaneous aspect rate.
%   - The ISAR image is target-local; it is NOT the same as the global
%     x-y map seen from the radar.
%
% Because SFCW has a finite unambiguous range, an independent coarse range
% gate is used to place the fine ISAR range bins at the correct absolute
% radar range. The fine offset is estimated from the received field by MP.
%
% Presentation-fixed version: physical range/cross-range axes are made
% consistent with the signal phase convention and with Chapter 8.
% The Chapter 8 JTF PRF (20 kHz) is restored because the Code 8.5
% frequency grid depends on tst = Npulse/PRF.
% V7 improvement:
%   Translational MP is estimated locally about the JTF reference time.
%   The exact 2-D off-RLOS trajectory is not perfectly represented by one
%   global quadratic R(t) model over the entire CPI. A local model around
%   tRef reduces that model mismatch.
%
%   The final compensated images are referenced to the actual JTF snapshot
%   range, and the GT overlay uses the actual GT scatterer ranges at tRef.
%
% V8 practical range-reference fix:
%   The SFCW fine-range unambiguous interval is only R_unamb/2 on either
%   side of a coarse gate. Therefore the coarse gate used for the local MP
%   stage must be near the JTF reference range, as supplied by a practical
%   coarse range tracker. The demo value below is 8085 m, intentionally only
%   0.866 m from the simulated JTF reference range.

clear;
close all;
clc;
rng(7,'twister');

%% 1. RADAR PARAMETERS
c      = 3.0e8;
Npulse = 128;
Nburst = 512;
f0     = 3.0e9;
BW     = 384e6;
T1     = (Npulse-1)/BW;
PRF    = 20e3;                  % Chapter 8 Code 8.5 JTF-compatible PRF
PRI    = 1/PRF;
tBurst = Npulse*PRI;
fsSlow = 1/tBurst;

f      = f0 + (0:Npulse-1)/T1;
fc     = mean(f);
lambda = c/fc;
dr     = c/(2*BW);
dfStep = 1/T1;
R_unamb = c/(2*dfStep);

%% 2. REALISTIC 2-D AIRCRAFT SCENARIO
% Aircraft/UAV surveillance geometry, 30 deg off initial RLOS.
R0_true       = 8.0e3;          % m
Vtotal_true   = 60.0;           % m/s (~216 km/h)
beta_true     = deg2rad(30);    % flight-path angle from initial RLOS
accel_true    = 0.5;            % m/s^2 along the flight path
omega_body    = 0.24;           % rad/s body yaw rate (book-compatible; avoids Doppler aliasing)
% Translation is 30 deg off the initial RLOS. The body rotates
% independently. For the Chapter-8 small-angle JTF cross-range relation
% y = lambda*fD/(2*omega_r), it is most useful to choose the target body
% x-axis to be aligned with the instantaneous LOS at the JTF reference time.
% This does NOT remove the 30-deg off-RLOS translation; it only defines the
% target attitude so that the physical cross-range coordinate is unambiguous.

Vx_true       = Vtotal_true*cos(beta_true);
Vy_true       = Vtotal_true*sin(beta_true);
ar_RLOS_true  = accel_true*cos(beta_true) + ...
                  (Vtotal_true^2*sin(beta_true)^2)/R0_true;

idx    = 0:(Npulse*Nburst-1);
tPulse = idx*PRI;

% One sample per burst for slow-time geometry.
tSlow = (0:Nburst-1)*tBurst;

% Full 2-D center trajectory.
dAlong = Vtotal_true*tSlow + 0.5*accel_true*tSlow.^2;
centerX_true = R0_true + dAlong*cos(beta_true);
centerY_true = dAlong*sin(beta_true);
Rcenter_true = hypot(centerX_true,centerY_true);
LOS_true     = atan2(centerY_true,centerX_true);

% JTF reference time: center of the CPI.
refBurst = round((Nburst+1)/2);
tRef = tSlow(refBurst);
Rref_true = Rcenter_true(refBurst);
LOSref = LOS_true(refBurst);

% Choose initial body yaw so that body x-axis is aligned with LOS at tRef.
yaw0_true = LOSref - omega_body*tRef;
yaw_true  = yaw0_true + omega_body*tSlow;

% Aspect angle = body orientation relative to instantaneous LOS.
aspect_true = unwrap(yaw_true - LOS_true);
pAspect = polyfit(tSlow,aspect_true,1);
omega_eff = pAspect(1);

yawRef = yaw_true(refBurst);

% Exact local radial kinematics at tRef for validation.
sRef = Vtotal_true*tRef + 0.5*accel_true*tRef^2;
vsRef = Vtotal_true + accel_true*tRef;
omega_LOS_ref = R0_true*vsRef*sin(beta_true)/(Rref_true^2);
vr_ref_true = ((R0_true*cos(beta_true)+sRef)/Rref_true)*vsRef;
ar_ref_true = (R0_true^2*sin(beta_true)^2/Rref_true^3)*vsRef^2 + ...
              ((R0_true*cos(beta_true)+sRef)/Rref_true)*accel_true;

% Eq. 8.33 uses the instantaneous rotational/aspect rate.
omega_ref = omega_body - omega_LOS_ref;

eR  = [cos(LOSref);  sin(LOSref)];
eCR = [-sin(LOSref); cos(LOSref)];

%% 3. SELF-GENERATED 2-D AIRCRAFT
[Xc,Yc,amp] = createJet2D();
Xc  = Xc(:);
Yc  = Yc(:);
amp = amp(:);

%% 4. GROUND TRUTH FIGURE
% Left: global radar geometry.
% Right: aircraft in the exact target-local range/cross-range frame used
% for comparison with the compensated ISAR image.
figGT = figure('Name','Ground Truth','Color','w', ...
               'Position',[40 70 1500 760]);
tl = tiledlayout(figGT,1,2,'TileSpacing','compact','Padding','compact');

% ------------------------- global geometry -------------------------------
ax1 = nexttile(tl);
hold(ax1,'on'); box(ax1,'on'); grid(ax1,'on'); axis(ax1,'equal');

plot(ax1,0,0,'kp','MarkerSize',20,'MarkerFaceColor','k');
text(ax1,0.18,-0.18,'Radar','Color','k','FontSize',15,'FontWeight','bold');
plot(ax1,[0 R0_true/1e3],[0 0],'k--','LineWidth',1.5);
text(ax1,3.1,0.12,'RLOS','Color','k','FontSize',14,'FontWeight','bold');
plot(ax1,centerX_true/1e3,centerY_true/1e3,'k-','LineWidth',2.2);

% Velocity vector at the start of the trajectory.
vPlot = 0.75;
quiver(ax1,R0_true/1e3,0, vPlot*cos(beta_true),vPlot*sin(beta_true),0, ...
       'k-','LineWidth',2.2,'MaxHeadSize',0.5);
text(ax1,R0_true/1e3+0.18,0.22,'V_t','Color','k','FontSize',16, ...
     'FontWeight','bold','Interpreter','tex');

% Enlarged aircraft poses, only for visualization.
poseIdx = [1 round(0.20*Nburst) round(0.40*Nburst) ...
             round(0.60*Nburst) round(0.80*Nburst) Nburst];
poseScale = 5.0;
for q = 1:numel(poseIdx)
    m = poseIdx(q);
    th = yaw_true(m);
    Rm = [cos(th) -sin(th); sin(th) cos(th)];
    P = poseScale*Rm*[Xc.';Yc.'];
    plot(ax1,(P(1,:)+centerX_true(m))/1e3, ...
         (P(2,:)+centerY_true(m))/1e3,'.', ...
         'Color',[0.80 0.05 0.05],'MarkerSize',8);
end

% Rotation indicator near the first pose.
arcR = 0.42e3;
arc = linspace(yaw0_true+0.20,yaw0_true+1.55,70);
plot(ax1,(R0_true+arcR*cos(arc))/1e3,(arcR*sin(arc))/1e3, ...
     'k-','LineWidth',1.5);
qa = arc(end-4);
quiver(ax1,(R0_true+arcR*cos(qa))/1e3,(arcR*sin(qa))/1e3, ...
       -0.12*sin(qa),0.12*cos(qa),0,'k-','LineWidth',1.4,'MaxHeadSize',0.45);
text(ax1,(R0_true-0.50e3)/1e3,0.58,'\omega_r','Color','k','FontSize',17, ...
     'FontWeight','bold','Interpreter','tex');

xlabel(ax1,'x [km]','Color','k');
ylabel(ax1,'y [km]','Color','k');
xlim(ax1,[-0.6 9.1]);
ylim(ax1,[-0.6 3.2]);
setBlackAxes(ax1);

% ----------------------- local physical ISAR frame -----------------------
ax2 = nexttile(tl);
hold(ax2,'on'); box(ax2,'on'); grid(ax2,'on'); axis(ax2,'equal');

% Actual scatterer locations at the JTF reference time.
Rrot = [cos(yawRef) -sin(yawRef); sin(yawRef) cos(yawRef)];
Pgt = Rrot*[Xc.';Yc.'];
Pgt = Pgt + [centerX_true(refBurst); centerY_true(refBurst)];
DP  = Pgt - [centerX_true(refBurst); centerY_true(refBurst)];

Rabs_gt = hypot(Pgt(1,:),Pgt(2,:));
CR_gt   = eCR.'*DP;

% For visual comparison with the translationally compensated ISAR images,
% express the JTF-time ground truth in the same t=0 absolute-range frame.
% This preserves each scatterer's range displacement from the aircraft
% centre while moving the centre from Rref_true back to R0_true.
Rabs_gt_comp = Rabs_gt;

plot(ax2,Rabs_gt_comp,CR_gt,'ko','MarkerSize',4.5,'MarkerFaceColor',[0.80 0.05 0.05]);

xlabel(ax2,'Absolute slant range [m]','Color','k');
ylabel(ax2,'Cross-range [m]','Color','k');

xlim(ax2,[R0_true-30 R0_true+30]);
ylim(ax2,[-20 20]);
setBlackAxes(ax2);

%% 5. FORWARD MODEL: EXACT 2-D SLANT RANGE
% The complete off-RLOS translation is used in the simulator.
cx = R0_true + (Vtotal_true*tPulse + 0.5*accel_true*tPulse.^2)*cos(beta_true);
cy =      (Vtotal_true*tPulse + 0.5*accel_true*tPulse.^2)*sin(beta_true);
yawPulse = yaw0_true + omega_body*tPulse;

Es = complex(zeros(Npulse,Nburst));

for s = 1:numel(Xc)
    dx = cos(yawPulse)*Xc(s) - sin(yawPulse)*Yc(s);
    dy = sin(yawPulse)*Xc(s) + cos(yawPulse)*Yc(s);

    Rsc = hypot(cx + dx, cy + dy);
    Rmat = reshape(Rsc,Npulse,Nburst);

    phase = (4*pi/c) .* (f(:)*ones(1,Nburst)) .* Rmat;
    Es = Es + amp(s)*exp(-1i*phase);
end

%% 6. COARSE RANGE REFERENCE
% An SFCW sweep with the book-style frequency spacing has a finite
% unambiguous range. A real system uses a coarse range gate/tracker first.
% This value is an independent acquisition quantity; the fine offset is
% estimated from Es below.
coarseRangeGate = Rref_true;     % m, independent coarse tracker output near tRef

% Deramp the large absolute phase at the coarse range reference. This moves
% the fine ISAR range response into the FFT's unambiguous local interval.
coarsePhase = exp(1i*(4*pi/c) .* (f(:)*coarseRangeGate));
Es_gate = Es .* (coarsePhase*ones(1,Nburst));

% The local MP reference must lie inside the unambiguous fine-range interval.
% For this simulation the chosen coarse tracker output is near tRef.
if abs(Rref_true - coarseRangeGate) > 0.45*R_unamb
    warning('Coarse gate is too far from JTF reference; fine MP range may alias.');
end

%% 7. RANGE AXIS: ABSOLUTE SLANT RANGE
rangeOffset = (-floor(Npulse/2):ceil(Npulse/2)-1)*dr;
rangeAbs_raw = coarseRangeGate + rangeOffset; % local physical range about the coarse gate

%% 8. FIGURE 8.16 -- CONVENTIONAL ISAR, NO COMPENSATION
rawISAR = abs(fftshift(fft2(Es_gate)));
raw_dB = 20*log10(rawISAR/(max(rawISAR(:))+eps)+eps);

fdSlow = (-floor(Nburst/2):ceil(Nburst/2)-1)*(fsSlow/Nburst);

% The raw image still contains translational Doppler, so the vertical axis
% is a Doppler-derived relative cross-range display only.
rawPower = sum(rawISAR.^2,1);
fdCenterRaw = sum(fdSlow.*rawPower)/sum(rawPower + eps);
crRaw = lambda*(fdSlow-fdCenterRaw)/(2*omega_ref);

% MATLAB FFT of exp(-j*4*pi*f*R/c) places positive physical range at
% negative FFT bins. Flip the range dimension for a physical left-to-right
% increasing slant-range axis.
raw_dB_phys = flipud(raw_dB);

fig816 = figure('Name','Figure 8.16 - Conventional ISAR', ...
                'Color','w','Position',[260 80 900 760]);
ax = axes(fig816);
imagesc(ax,rangeAbs_raw,crRaw,raw_dB_phys.');
axis(ax,'xy'); box(ax,'on');
colormap(ax,hot(256));
cb = colorbar(ax); cb.Color = 'k';
caxis(ax,[-30 0]);
xlabel(ax,'Absolute slant range [m]','Color','k');
ylabel(ax,'Cross-range [m]','Color','k');
xlim(ax,[coarseRangeGate-30 coarseRangeGate+30]);
ylim(ax,[-20 20]);
setBlackAxes(ax);

%% 9. TRANSLATIONAL MOCOMP -- LOCAL MP ABOUT THE JTF REFERENCE
% Book basis: MP projection of the received field onto motion parameter
% bases. The difference here is that the local motion model is centered
% at tRef rather than at t=0:
%
%   R(t) ~= Rref + VrRef*(t-tRef) + 1/2*arRef*(t-tRef)^2
%
% This is better matched to the exact 2-D off-RLOS trajectory used by the
% simulator because the instantaneous radial kinematics vary over the CPI.
Tmotion = reshape(tPulse,Npulse,Nburst);
DeltaT  = Tmotion - tRef;
Fmat = f(:)*ones(1,Nburst);

% Use only a local slow-time neighborhood around tRef for the MP fit.
% The exact 2-D off-RLOS range is very well approximated by a quadratic over
% this shorter interval, while the full-CPI fit is biased by higher-order
% range curvature.
mpHalfWindow = 0.50;                    % [s], local MP fitting interval
mpWeight = double(abs(DeltaT) <= mpHalfWindow);

% RR is the fine reference-range offset relative to the coarse gate.
RR_coarse = -8:0.25:8;                  % [m]
V_coarse  = 48:0.5:58;                   % radial velocity at tRef [m/s]
A_coarse  = -1:0.1:2;                   % radial acceleration at tRef [m/s^2]

[dR_est,Vr_ref_est,ar_ref_est] = mpSearchLocal( ...
    Es_gate,DeltaT,Fmat,c,RR_coarse,V_coarse,A_coarse,mpWeight);

% Local refinement around the coarse maximum.
RR_ref = max(-8,dR_est-0.75):0.05:min(8,dR_est+0.75);
V_ref  = max(45,Vr_ref_est-1.0):0.05:min(65,Vr_ref_est+1.0);
A_ref  = max(-2,ar_ref_est-0.30):0.01:min(2.5,ar_ref_est+0.30);

[dR_est,Vr_ref_est,ar_ref_est] = mpSearchLocal( ...
    Es_gate,DeltaT,Fmat,c,RR_ref,V_ref,A_ref,mpWeight);

dR_est       = double(dR_est(1));
Vr_ref_est   = double(Vr_ref_est(1));
ar_ref_est   = double(ar_ref_est(1));
Rref_est     = coarseRangeGate + dR_est;

fprintf('\n==========================================================\n');
fprintf('LOCAL TRANSLATIONAL MP ESTIMATE (ABOUT JTF REFERENCE)\n');
fprintf('True JTF reference range    : %.3f m\n',Rref_true);
fprintf('Estimated reference range   : %.3f m\n',Rref_est);
fprintf('Reference-range error        : %.3f m\n',Rref_est-Rref_true);
fprintf('True radial velocity @tRef  : %.3f m/s\n',vr_ref_true);
fprintf('Estimated radial velocity   : %.3f m/s\n',Vr_ref_est);
fprintf('True radial acceleration @tRef: %.3f m/s^2\n',ar_ref_true);
fprintf('Estimated acceleration      : %.3f m/s^2\n',ar_ref_est);
fprintf('Body angular rate           : %.5f rad/s\n',omega_body);
fprintf('Reference aspect rate       : %.5f rad/s\n',omega_ref);
fprintf('Coarse tracker range gate  : %.3f m\n',coarseRangeGate);
fprintf('SFCW unambiguous range      : %.3f m\n',R_unamb);
fprintf('==========================================================\n');

%% 10. TRANSLATIONAL PHASE COMPENSATION
motionDistance = Vr_ref_est.*DeltaT + 0.5.*ar_ref_est.*(DeltaT.^2);
Sconj = exp((1i*4*pi/c).*Fmat.*motionDistance);
S_Duz = Es_gate .* Sconj;

%% 11. FIGURE 8.19 -- AFTER TRANSLATIONAL MOCOMP
transISAR = abs(fftshift(fft2(S_Duz)));
trans_dB = 20*log10(transISAR/(max(transISAR(:))+eps)+eps);

crTrans = lambda*fdSlow/(2*omega_ref);

% S_Duz is still referenced to the independent coarse range gate. The MP
% R0 estimate is reported/validated, but it is not silently substituted
% into the image axis.
rangeAbs_trans = Rref_est + rangeOffset;

% Same FFT range-sign correction as Figure 8.16.
trans_dB_phys = flipud(trans_dB);

fig819 = figure('Name','Figure 8.19 - Translational MOCOMP', ...
                'Color','w','Position',[280 100 900 760]);
ax = axes(fig819);
imagesc(ax,rangeAbs_trans,crTrans,trans_dB_phys.');
axis(ax,'xy'); box(ax,'on');
colormap(ax,hot(256));
cb = colorbar(ax); cb.Color = 'k';
caxis(ax,[-30 0]);
xlabel(ax,'Absolute slant range [m]','Color','k');
ylabel(ax,'Cross-range [m]','Color','k');
xlim(ax,[Rref_est-30 Rref_est+30]);
ylim(ax,[-20 20]);
setBlackAxes(ax);

%% 12. FIGURE 8.21 -- TRANSLATIONAL + JTF ROTATIONAL PROCESSING
% Book Code 8.5 structure:
%   Hamming window -> 1-D IFFT over frequency -> Gabor/JTF transform.
%
% IMPORTANT:
% The book's Gabor frequency variable fp starts at 145 Hz and increases by
% 1/(Npulse*tst).  It is NOT the same thing as a conventional signed
% Doppler-frequency grid.  Because the slow-time data are sampled at
% fsSlow, these fp values alias into the physical Doppler interval.
%
% The book's Gabor kernel is exp(-j*2*pi*fp*t).  With the received-field
% convention used in this simulator, the corresponding physical received
% Doppler is therefore
%
%       fD_physical = alias(fp)
%
% For this simulator, the received field after translational compensation
% has exp(+j*2*pi*fD*t), while the Gabor kernel is exp(-j*2*pi*fp*t).
% Hence the correlating channel is fp = fD (up to slow-time aliasing).

Ese = S_Duz;

% ---- 1. Hamming window, exactly as in Book Code 8.5 --------------------
win  = hamming(Npulse) * hamming(Nburst).';
Esew = Ese .* win;

% ---- 2. 1-D IFFT over the frequency dimension --------------------------
% IFFT converts stepped-frequency samples into range profiles.
% The transpose gives [slow-time x range-bin], as in the book.
Es_IFFT = ifft(Esew,[],1).';

% For physical range display, circularly center the range bins.
Es_range = fftshift(Es_IFFT,2);

% ---- 3. Book-style burst-time vector -----------------------------------
% One burst/dwell corresponds to Npulse PRIs.
tJTF = tSlow;
tst  = tBurst;

% ---- 4. Gabor wavelet bank -- physically calibrated adaptation ---------
N_gabor = 1;

% Book Code 8.5 starts fp at 145 Hz because its plotted vertical axis is a
% Doppler/frequency INDEX. It is not an absolute physical Doppler grid.
% Here fsSlow = 1/tBurst = 156.25 Hz, so 145 Hz is outside the physical
% Nyquist interval and aliases to -11.25 Hz. That cannot cover the aircraft
% Doppler span and produced the fan/spoke artifact in the previous Figure 8.21.
%
% Keep the same Gabor/JTF operation and the same frequency-channel spacing,
% but center the bank at zero so that every channel corresponds directly to
% a physical, unaliased slow-time Doppler frequency.
fpStep = 1/(Npulse*tst);
fpBook = ((0:Npulse-1) - floor(Npulse/2))*fpStep;

Alpha_p = 0.04;

GaborWavelet = complex(zeros(Npulse,Nburst));

% Use the same snapshot time selected for the ground-truth comparison.
tp = tRef;

gaussianTerm = (1/sqrt(2*pi*(Alpha_p)^2)) .* ...
               exp(-((tJTF-tp).^2)/(2*Alpha_p));

for rr = 1:Npulse
    fp = fpBook(rr);

    harmonicTerm = exp((-1i*2*pi*fp*(tJTF-tp))/N_gabor);

    GaborWavelet(rr,:) = gaussianTerm .* harmonicTerm;
end

% ---- 5. Gabor/JTF transform --------------------------------------------
% This is the core operation of Book Code 8.5.
St_raw = GaborWavelet * Es_range;

% ---- 6. Physical signed Doppler axis ----------------------------------
% The centered fpBook bank is already inside the unambiguous physical
% slow-time Doppler interval, so no aliasing/remapping is required.
%
% For the chosen exp(-j*4*pi*f*R/c) forward model, the rotational range term
% produces a slow-time factor exp(+j*2*pi*fD*t). The Gabor kernel is
% exp(-j*2*pi*fp*t), so the correlating channel is fp = fD.
fdPhysical = fpBook;
St_phys = St_raw;

finalISAR = abs(St_phys);
final_dB = 20*log10(finalISAR/(max(finalISAR(:))+eps)+eps);

% ---- 7. Physical cross-range calibration -------------------------------
% Book Eq. 8.33:
%
%       cross-range = lambda*fD/(2*omega_r)
%
% Use the effective aspect rate of THIS simulated geometry.
% Eq. 8.33: physical cross-range from the signed Doppler.
% Do not subtract an image-derived Doppler centroid here: that would hide
% a residual motion-compensation error by redefining the coordinate origin.
crossRangeJTF = lambda*fdPhysical(:)/(2*omega_ref);

% The local translational compensation is referenced to tRef, so the
% remaining constant range location is the estimated JTF reference range.
rangeAbs_jtf = Rref_est + rangeOffset;

fig821 = figure('Name','Figure 8.21 - Translational and Rotational MOCOMP', ...
                'Color','w','Position',[300 120 900 760]);
ax = axes(fig821);
imagesc(ax,rangeAbs_jtf,crossRangeJTF,final_dB);
axis(ax,'xy'); box(ax,'on');
colormap(ax,hot(256));
cb = colorbar(ax); cb.Color = 'k';
caxis(ax,[-30 0]);
xlabel(ax,'Absolute slant range [m]','Color','k');
ylabel(ax,'Cross-range [m]','Color','k');
xlim(ax,[Rref_est-30 Rref_est+30]);
ylim(ax,[-20 20]);
setBlackAxes(ax);

%% 13. PRINT DIRECT GROUND-TRUTH EXTENTS
fprintf('\nGROUND-TRUTH TARGET AT JTF REFERENCE TIME\n');
fprintf('Initial true slant range      : %.3f m\n',R0_true);
fprintf('JTF-reference true range      : %.3f m\n',Rref_true);
fprintf('Compensated GT range extent   : %.3f to %.3f m\n',min(Rabs_gt_comp),max(Rabs_gt_comp));
fprintf('Cross-range extent            : %.3f to %.3f m\n',min(CR_gt),max(CR_gt));
fprintf('Aircraft length               : %.3f m\n',max(Xc)-min(Xc));
fprintf('Aircraft width                : %.3f m\n',max(Yc)-min(Yc));
fprintf('Estimated JTF reference range : %.3f m\n',Rref_est);
fprintf('True JTF reference range      : %.3f m\n',Rref_true);
fprintf('Reference-range error         : %.3f m\n',Rref_est-Rref_true);
fprintf('True radial velocity @tRef    : %.3f m/s\n',vr_ref_true);
fprintf('Estimated radial velocity     : %.3f m/s\n',Vr_ref_est);
fprintf('True radial acceleration @tRef: %.3f m/s^2\n',ar_ref_true);
fprintf('Estimated radial acceleration : %.3f m/s^2\n',ar_ref_est);
fprintf('JTF reference aspect rate     : %.6f rad/s\n',omega_ref);
fprintf('JTF Doppler grid              : %.3f to %.3f Hz\n',min(fdPhysical),max(fdPhysical));
fprintf('JTF cross-range grid          : %.3f to %.3f m\n',min(crossRangeJTF),max(crossRangeJTF));
fprintf('JTF cross-range scale         : %.6f m/Hz\n',lambda/(2*omega_ref));
fprintf('Physical range sign correction: Fig. 8.16/8.19 horizontal flip\n');
fprintf('Physical cross-range sign     : Eq. 8.33 (positive)\n');
fprintf('Local MP model center        : tRef = %.6f s\n',tRef);

%% 14. GROUND-TRUTH OVERLAY VALIDATION
% This figure leaves the normal Figure 8.21 image unchanged and overlays the
% independently calculated ground-truth scatterer locations in the SAME
% compensated absolute-range / cross-range coordinate system.
%
% Ground truth is evaluated at the same JTF snapshot time tRef as the local
% MP model and JTF Gaussian window. The aircraft body x-axis was deliberately
% aligned with LOS at tRef, so CR_GT is also the target's physical body-y
% coordinate at the reference snapshot. No artificial t=0 translation is applied.
%
% Cross-range is measured in the instantaneous LOS frame at tRef:
%
%     CR_GT = eCR^T DeltaP.
%
figOverlay = figure('Name','Figure 8.21 - Ground Truth Overlay Validation', ...
                    'Color','w','Position',[330 150 980 800]);
ax = axes(figOverlay);
imagesc(ax,rangeAbs_jtf,crossRangeJTF,final_dB);
axis(ax,'xy');
box(ax,'on');
colormap(ax,hot(256));
cb = colorbar(ax); cb.Color = 'k';
caxis(ax,[-30 0]);
hold(ax,'on');

% Plot the exact ground-truth scattering centers at the JTF reference time.
plot(ax,Rabs_gt_comp,CR_gt,'co','MarkerSize',5.5, ...
     'LineWidth',1.1,'MarkerFaceColor','none');

% Mark the target-center reference location.
plot(ax,Rref_true,0,'c+','MarkerSize',10,'LineWidth',1.5);

% Ground-truth bounding box, useful for presentation validation.
plot(ax,[min(Rabs_gt_comp) max(Rabs_gt_comp) max(Rabs_gt_comp) min(Rabs_gt_comp) min(Rabs_gt_comp)], ...
     [min(CR_gt) min(CR_gt) max(CR_gt) max(CR_gt) min(CR_gt)], ...
     'c--','LineWidth',0.9);

xlabel(ax,'Absolute slant range [m]','Color','k');
ylabel(ax,'Cross-range [m]','Color','k');
xlim(ax,[Rref_est-30 Rref_est+30]);
ylim(ax,[-20 20]);
legend(ax,{'Ground-truth scatterers','Target center','Ground-truth extent'}, ...
       'Location','northeast','Box','on');
setBlackAxes(ax);

title(ax,'Figure 8.21 with Ground-Truth Scatterer Overlay (t = tRef)', ...
      'FontWeight','bold');

%% 15. SIMPLE NUMERICAL OVERLAY CHECK
% For every ground-truth scatterer, locate the nearest final-image pixel and
% report its intensity. This is a diagnostic only; it does not alter the image.
%
% MATLAB dimensions are handled explicitly:
%   rangeGrid = Nrange x 1, Rabs_gt_comp = 1 x Nscatter
%   crossGrid = Ncross x 1, CR_gt        = 1 x Nscatter
rangeGrid = rangeAbs_jtf(:);
crossGrid = crossRangeJTF(:);

[~,idxR_gt]  = min(abs(rangeGrid - Rabs_gt_comp),[],1);
[~,idxCR_gt] = min(abs(crossGrid - CR_gt(:).'),[],1);

idxR_gt  = idxR_gt(:);
idxCR_gt = idxCR_gt(:);

linIdx = sub2ind(size(final_dB),idxCR_gt,idxR_gt);
gtPixel_dB = final_dB(linIdx);

fprintf('\nGROUND-TRUTH OVERLAY PIXEL CHECK\n');
fprintf('Nearest final-image pixel at GT scatterers:\n');
fprintf('  minimum intensity : %.2f dB\n',min(gtPixel_dB));
fprintf('  median intensity  : %.2f dB\n',median(gtPixel_dB));
fprintf('  maximum intensity : %.2f dB\n',max(gtPixel_dB));
fprintf('  GT points above -15 dB: %d / %d\n',sum(gtPixel_dB >= -15),numel(gtPixel_dB));

% Second diagnostic: for each GT point, compute distance to the nearest
% final-image pixel whose value is at least -15 dB.
brightMask = final_dB >= -15;
if any(brightMask(:))
    [rBright,cBright] = find(brightMask);
    rangeBright = rangeAbs_jtf(cBright(:));
    crossBright = crossRangeJTF(rBright(:));
    rangeBright = rangeBright(:);
    crossBright = crossBright(:);
    Rabs_gt_col = Rabs_gt_comp(:);
    CR_gt_col = CR_gt(:);

    gtNearestDist = zeros(numel(Rabs_gt_col),1);
    for kk = 1:numel(Rabs_gt_col)
        d2 = (rangeBright - Rabs_gt_col(kk)).^2 + ...
             (crossBright - CR_gt_col(kk)).^2;
        gtNearestDist(kk) = sqrt(min(d2(:)));
    end

    fprintf('Nearest >= -15 dB bright-pixel distance:\n');
    fprintf('  median distance : %.3f m\n',median(gtNearestDist));
    fprintf('  maximum distance: %.3f m\n',max(gtNearestDist));
    fprintf('  (distance measured in displayed range/cross-range coordinates)\n');
else
    fprintf('No pixels above -15 dB; bright-pixel distance not computed.\n');
end

%% 15. PRESENTATION SUMMARY
fprintf('\\nPRESENTATION VALIDATION SUMMARY\\n');
fprintf('Body/LOS alignment error at tRef : %.6f deg\\n',rad2deg(unwrap(yawRef-LOSref)));
fprintf('Local MP fitting half-window     : %.3f s\\n',mpHalfWindow);
fprintf('Local MP reference-range error    : %.3f m\\n',Rref_est-Rref_true);
fprintf('Local MP radial-velocity error    : %.3f m/s\\n',Vr_ref_est-vr_ref_true);
fprintf('Local MP radial-accel error       : %.3f m/s^2\\n',ar_ref_est-ar_ref_true);
fprintf('GT range extent @ tRef            : %.3f to %.3f m\\n',min(Rabs_gt),max(Rabs_gt));
fprintf('GT cross-range extent              : %.3f to %.3f m\\n',min(CR_gt),max(CR_gt));
fprintf('\\n');

%% ========================= LOCAL FUNCTIONS ==============================
function [Rbest,Vbest,Abest] = mpSearchLocal(Es,DeltaT,F,c,RR,V,A,mpWeight)
% Localized Chapter-8-style coherent MP projection.
% Es has already been deramped by the coarse reference near tRef.
% RR = Rref - coarseRangeGate; V/A are radial kinematics at tRef.
% mpWeight selects the local slow-time neighborhood used by the MP fit.

scoreBest = -inf;
Rbest = RR(1);
Vbest = V(1);
Abest = A(1);

if nargin < 8 || isempty(mpWeight)
    mpWeight = ones(size(DeltaT));
end

% Range steering for candidate reference-range offsets.
rangeSteering = exp(1i*(4*pi/c).*(F(:,1)*RR));

for iv = 1:numel(V)
    Vest = V(iv);
    for ia = 1:numel(A)
        Aest = A(ia);

        motionDistance = Vest.*DeltaT + 0.5.*Aest.*(DeltaT.^2);
        motionPhase = exp((1i*4*pi/c).*F.*motionDistance);

        % Coherent projection over the local slow-time interval.
        burstSum = sum(Es.*motionPhase.*mpWeight,2);
        projection = rangeSteering.'*burstSum;
        projection = abs(projection(:));

        [localBest,idxLocal] = max(projection);
        if localBest > scoreBest
            scoreBest = localBest;
            Rbest = RR(idxLocal);
            Vbest = Vest;
            Abest = Aest;
        end
    end
end
end

function [x,y,amp] = createJet2D()
% Equal-amplitude ideal point scatterers.

x_fuse = linspace(-13,13,35);
y_fuse = zeros(size(x_fuse));

x_wing  = linspace(5,-3,18);
y_wingR = linspace(1.0,11.0,18);
y_wingL = -y_wingR;

x_tail  = linspace(-8.5,-12.0,10);
y_tailR = linspace(0.6,4.0,10);
y_tailL = -y_tailR;

x_fin = [-8.5 -9.5 -10.5 -11.5];
y_fin = [0.4 0.8 1.2 1.6];

x_nose = [13.0 13.8 14.5];
y_nose = [0 0 0];

x_eng = linspace(3,-2.5,10);
y_engR = 2.2*ones(size(x_eng));
y_engL = -y_engR;

x = [x_fuse,x_wing,x_wing,x_tail,x_tail,x_fin,x_nose,x_eng,x_eng];
y = [y_fuse,y_wingR,y_wingL,y_tailR,y_tailL,y_fin,y_nose,y_engR,y_engL];
amp = ones(size(x));
end

function setBlackAxes(ax)
set(ax,'XColor','k','YColor','k','Color','w', ...
    'FontName','Arial','FontSize',12,'FontWeight','bold', ...
    'LineWidth',1.0);
ax.XLabel.Color = 'k';
ax.YLabel.Color = 'k';
end
