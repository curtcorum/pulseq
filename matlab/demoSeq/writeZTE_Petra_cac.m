% a basic 3D ZTE/PETRA sequence
% achieves "TE" about 70 us and possibly below (100 us on GE TOPPE for now)
%
% recommend executing >> clear all;
% before running, especially for TOPPE conversion, which reads variables
% from the workspace.
%
% started from https://github.com/pulseq/pulseq/blob/master/matlab/demoSeq/writeZTE_Petra.m
%
% Curt Corum, 3/19/2024

%% high-level sequence parameters

fov     = 2*256e-3;     % 2x os FOV
dx      = 16e-3;         % Define FOV and resolution 16 for 626 blocks ( only spi center), and 4 for 9765 blocks (nr=256)
alpha   = 4;            % flip angle
Nr      = 256;          % number of readout points (with some oversampling), saved to Nx
R       = 16;           % acceleration/undersampling (for the outer shell)
R_inner = 1;            % acceleration/undersampling (angular direction of the inner area)
%xSpoil=3;%0.6;          %0.6 for 2mm % the amount of spoiling after the end of the readout (used to ramp to the next point
tSpoil  = .00088;       % s of spoiling after ro, set to make TR 5ms for now *** 
do_spi = false;          % do_spi = true; generate SPI blocks for PETRA, only ZTE if false
do_zte = true;
ge_dt = 00e-6;         % 116e-6 s (rounded to 120us) refresh time for loading new waveforms on GE

Kmax=1/2/dx;
dK=1/fov;

%% more detailed and derived params

rf_duration         = 16e-6;        % covers inner half of FOV at ~90% flat % duration of the excitation pulse
ro_duration         = 16*Nr*1e-6;   % 2x os 31.25 kHz GE % read-out time: controls RO bandwidth and T2-blurring
minRF_to_ADC_time   = 10e-6;       % the parameter wich defines TE together with ro_discard
%ro_discard=4;                       % how many ADC samples are contaminated by RF switching artifacts and alike
rfSpoilingInc       = 117;          % RF spoiling increment
nDummyShots         = 0;            % TOPPE parameter, ~= 0 not implemented yet
pislquant           = 0;            % TOPPE parameter, number of shots/ADC events used for receive gain calibration, ~= 0 not implemented yet
%Ny                  = Ns + Ns_spi; % TOPPE parameter, data index, need to calculate at bottom, Ny = Ns + Ns_spi

% TOPPE System/design parameters from v6/examples/write2DGRE.m
% Reduce gradients by 1/sqrt(3) to allow for oblique scans.
% Reduce slew a bit further to reduce PNS.
% also some previously unspecified parameters, CAC 240314 
sys = mr.opts('maxGrad',        50/sqrt(3), 'gradUnit','mT/m', ...
              'maxSlew',        120/sqrt(3), 'slewUnit', 'T/m/s', ...
              'riseTime',       [],...
              'rfDeadTime',     10e-6, ... % 100e-6
              'rfRingdownTime', 0e-6, ... % 60e-6
              'adcDeadTime',    0e-6, ... % 50e-6
              'adcRasterTime',  2e-6, ... % 2e-6
              'rfRasterTime',   1e-6,... % 1e-6
              'gradRasterTime', 20e-6,... % 20e-6, set equal to blockDurationRaster
              'blockDurationRaster', 20e-6, ... % 20e-6
              'gamma',          42576000,...
              'B0',             3.0);

seq     = mr.Sequence(sys);         % Create a new sequence object
seq_sar = mr.Sequence(sys);         % Create an auxillary sequence object for SAR testing

%% create main sequence elements

% %create alpha-degree block pulse 
rf = mr.makeBlockPulse( alpha * pi/180, 'Duration', rf_duration, 'system', sys); % **************************************** hard pulse ********************
rf_dummy = rf; rf_dummy.signal = rf.signal * 1e-6; % for debug, used declaration in function
rf_align_delay = ceil( mr.calcDuration( rf)/sys.blockDurationRaster)*sys.blockDurationRaster; % calculated again in the populate_subsequence function for now...
rf_aux = rf;
%or create alpha-degree gaussian pulse 
%rf = mr.makeGaussPulse(alpha*pi/180,'Duration',rf_duration,'timeBwProduct',3,'system',sys);

Tenc = rf_duration/2 + minRF_to_ADC_time + ro_duration; %encoding time
% Tg=sys.rfDeadTime+rf_duration/2+Tenc+sys.adcDeadTime; % constant gradient time
Tg = ceil( (sys.rfDeadTime+rf_duration/2 + Tenc + sys.adcDeadTime) / sys.gradRasterTime) * sys.gradRasterTime; % constant gradient time, additional rounding *** CAC 240414
%Tt=ceil((Tenc*(1+xSpoil)-Tg)/sys.gradRasterTime)*sys.gradRasterTime; % transition time
%Tt = ceil( ((Tenc + tSpoil) - Tg) / sys.gradRasterTime) * sys.gradRasterTime; % transition time % *** CAC 230314
Tt = ceil( ((Tenc + tSpoil) - Tg) / sys.blockDurationRaster) * sys.blockDurationRaster; % transition time % *** CAC 230314
assert( Tt > 0, 'Transition time too short, possibly add spoiling.');
rf_aux.delay = Tt - rf_align_delay + 14e-6;

Ag = Kmax/Tenc; % read gradient grdient amplitude
% % Gr=mr.makeTrapezoid('z','Amplitude',Ag,'flatTime',Tg,'riseTime',0); % this graient has no ramps, I wonder if the Matlab mr library would like it...
% Gr=mr.makeExtendedTrapezoid('z','times',[0 Tg],'amplitudes',[Ag Ag]); % this is a constant graient with no ramps

TR = Tg + Tt; % gradient defined in populate_subsequence()

adc = mr.makeAdc( Nr,'Duration', ro_duration, 'Delay', sys.rfDeadTime + rf_duration + minRF_to_ADC_time);
adc_align_delay = ceil( mr.calcDuration( adc)/sys.blockDurationRaster) *  sys.blockDurationRaster;

Gr = mr.makeExtendedTrapezoid( 'z','times', [0 TR-Tt], 'amplitudes', [Ag Ag]); % this is a constant graient with no ramps, *** defined in function, for debugging here *** CAC 240320
Gdummy = Gr; % for SPI center, needed for TOPPE

SamplesBookkeeping=[];

%rfbw=1/rf_duration;
rfbw = mr.calcRfBandwidth( rf);
fprintf( 'Pulse bandwidth %g [kHz]\n', rfbw/1e3);
fprintf( 'Read gradient amplitude %g mT/m, effective "slice thinckess" %f mm\n', 1e3*Ag/sys.gamma, 1/Ag*rfbw*1e3);

FO = 50e-3 * Ag;
FN = 0; % TR is broken now for FN~=0... % ???? slice index *** CAC 240314




%% generate the sampling set on a surface of a sphere
[phi, theta, im]    = spherical_samples( Kmax, dK, R); 
Ns                  = length( phi);
SamplesBookkeeping  = [SamplesBookkeeping Ns];

%% main ZTE loop
if do_zte
    fprintf( 'Populating ZTE loop (%g TRs * %g ms/TR = %g s)\n', Ns, TR*1000, TR*Ns);
    tic;
    trid_n = 5;
    populate_subsequence( sys, seq, rf, adc, phi, theta, im, Ns, Ag, TR, Tt, FN, FO, trid_n, ge_dt); % ZTE =============================================================================== ZTE
    toc
else
    fprintf( 'ZTE section disabled...\n');
    Ns = 0;
end

%% SPI loop
Ns_spi = 0;  % default to 0 if loop does not run
if do_spi
    KstartZTE   = Ag*(rf_duration/2 + minRF_to_ADC_time + adc.dwell);
    nKspi       = floor( KstartZTE/dK);
    dKspi       = KstartZTE/(nKspi + 1);
    Tenc_spi    = rf_duration/2 + minRF_to_ADC_time + adc.dwell/2; % name changed from Tenc so no overwrite from ZTE section, *** CAC 240314
    fprintf( 'Populating SPI loop (%d spheres)\n', nKspi);
    tic;
    trid_n = 1;
    for s = nKspi:-1:1
        % generate the sampling set on a surface of a sphere
        [phi, theta, im] = spherical_samples( dKspi*s, dKspi, R_inner); % normally no acceleration
        Ns_spi = length( phi);
        SamplesBookkeeping = [SamplesBookkeeping Ns_spi];
        % update gradients
        Ag = dKspi*s/Tenc_spi; % read gradient grdient amplitude
        % actually create the next sampling sphere
        fprintf(' Populating sphere %d (%d TRs)\n', s, Ns_spi);
        fprintf( 'Effective "slice thinckess" %f mm\n',  1/Ag*rfbw*1e3);
        populate_subsequence( sys, seq, rf, adc, phi, theta, im, Ns_spi, Ag, TR, Tt, FN, FO, trid_n, ge_dt); % SPI ======================================================================= SPI
    end

    % sample the centere of k-space
    Ns_spi = Ns_spi +1;
    fprintf(' Populating center point (1 TR)\n');
    fprintf( 'Effective "slice thinckess" INF mm\n');
    trid_n = 3; 
    %rf_cent = rf; rf_cent.delay = Tt - rf_cent.shape_dur;
    seq.addBlock( mr.makeDelay( Tt + rf_align_delay), rf_aux, mr.makeLabel( 'SET', 'TRID', trid_n)); % rf and adc cannot be in same block for toppe, 
    seq.addBlock( adc, mr.scaleGrad(Gdummy, 0)); % SPI center ===================================================================================================================== SPI center
    SamplesBookkeeping = [SamplesBookkeeping 1];
    toc
    fprintf( 'Total number of SPI samples: %d; a Cartesian patch would require %d\n', sum( SamplesBookkeeping(2:end)), ceil( nKspi^3*pi*4/3));
else
    fprintf( 'SPI section disabled...\n');
    SamplesBookkeeping = [SamplesBookkeeping 0 0];
end

%% duration and number of blocks
seq_nblocks = size( seq.blockDurations); 
fprintf('Sequence total number of blocks %d\n', seq_nblocks(2));
seq_duration = seq.duration;
fprintf('Sequence total duration %g [s]\n', seq_duration);


%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%%
%seq.plot('TimeRange',[0 1]);

seq.setDefinition( 'FOV', [fov fov fov]);
seq.setDefinition( 'Name', 'ztePetraCac');
seq.setDefinition( 'SamplesPerShell', SamplesBookkeeping);

% for TOPPE
seq.setDefinition( 'nDummyShots', nDummyShots);
seq.setDefinition( 'pislquant', pislquant);
Nx = Nr; 
seq.setDefinition( 'Nx', Nx);
Nv = sum( SamplesBookkeeping); % always does center SPI
Ny = Nv;
seq.setDefinition( 'Ny', Ny);
%seq.setDefinition( 'Nz', Nr);
Nv = sum( SamplesBookkeeping); % always does center SPI
seq.setDefinition( 'Nv', Nv);

seq.write( 'ztePetraCac.seq');

% return % to suppress SAR calc and plots

%% create an RF-only version of the sequence (e.g. for the SAR or signal evolution testing)

tic;
[total_duration, total_numBlocks]=seq.duration();
for iB=1:total_numBlocks
    b=seq.getBlock(iB);
    bd=seq.blockDurations(iB);
    bs={mr.makeDelay(bd)};
    if ~isempty(b.rf)
        bs{end+1}=b.rf;
    end
    if ~isempty(b.adc)
        bs{end+1}=b.adc;
    end
    seq_sar.addBlock(bs);
end
toc
seq_sar.write('ztePetraCacSar.seq');

%return %uncomment to supress plots, etc.

% %% test binary storing
% 
% seq.writeBinary('zte_petra.bin');
% seq_bin=mr.Sequence();          
% seq_bin.readBinary('zte_petra.bin');
% seq_bin.write('zte_petra_bin.seq');
% return

%% visualize the 3D k-space 
tic;
[kfa,~,kf]=seq.calculateKspacePP();
toc

% K-Space Sampling Trajectory
figure; plot3(kf(1,:),kf(2,:),kf(3,:));
hold on; plot3(kfa(1,:),kfa(2,:),kfa(3,:),'r.');
%figure; plot3(kfa(1,:),kfa(2,:),kfa(3,:));

% Sequence plots
seq.plot

% Nice figure plots
seq.paperPlot

%% local functions

%**********************************************************************************************************************
% function [phi, theta, im]=spherical_samples( Kr, dK, R)
function [phi, theta, im]=spherical_samples( Kr, dK, R)
% the number of samples equals the ceil of the area of the sphere divided
% by the area around every sample
Ns=ceil(4*pi*((Kr/dK)^2)/R); 
np=0:(Ns-1);
alpha_gold=pi*(3-sqrt(5));
phi=np*alpha_gold;
%theta=0.5*pi*sqrt(np/Ns); % from Davide Piccini  https://doi.org/10.1002/mrm.22898
theta=acos(1-2*np/(Ns-1));  % from  Anton Semechko (2020). Suite of functions to perform uniform sampling of a sphere (https://github.com/AntonSemechko/S2-Sampling-Toolbox), GitHub. Retrieved October 3, 2020. 

xp=sin(theta).*cos(phi);
yp=sin(theta).*sin(phi);
zp=cos(theta);

figure; sphere; colormap gray;
hold on; plot3(xp,yp,zp,'.');

% looking for the optimal interleaving factor
nm=round(Ns/2); % middle of the trajetory
sr=round(sqrt(Ns));
v0=[xp(nm); yp(nm); zp(nm)];
v=[xp(nm+(1:sr)); yp(nm+(1:sr)); zp(nm+(1:sr))];
%figure;plot(vecnorm((v-v0(:,ones(size(v,2),1)))));
[dKm, im] = min( vecnorm( (v-v0(:, ones( size( v, 2), 1)))));

% fprintf('requested dK=%g achieved minimal dK=%g, min acceleration: %g, Ns=%d\n', dK/sqrt(R), dKm*Kr, (dKm*Kr/dK)^2,Ns);

% v=[xp; yp; zp]*Kr/dK/sqrt(R)*100;
% md=zeros(1,Ns);
% d=zeros(1,Ns);
% for i=1:Ns
%     d=vecnorm(v(:,i*ones(1,Ns))-v);
%     d(i)=NaN;
%     md(i)=min(d);
% end
% fprintf('dKmin=%g%%, dKmed=%g%% dKmax=%g%%\n', min(md),median(md),max(md));

end % function [phi, theta, im]=spherical_samples( Kr, dK, R) *********************************************************


%**********************************************************************************************************************
% function populate_subsequence( sys, seq, rf, adc, phi, theta, im, Ns, Ag, TR, Tt, FN, FO, trid_n, ge_dt)
function populate_subsequence( sys, seq, rf, adc, phi, theta, im, Ns, Ag, TR, Tt, FN, FO, trid_n, ge_dt)
Azc=Ag*(TR-Tt)/(TR+Tt);  %*0.35;
%Azc=0; % for MoCo we need "no-gradient" event blocks to be able to apply updates

rf_align_delay = ceil( mr.calcDuration( rf)/sys.blockDurationRaster)*sys.blockDurationRaster;
rf.delay = rf_align_delay - mr.calcDuration( rf);
%rf_dummy = rf; rf_dummy.signal = rf.signal * 1e-6; % *** dummy pulse for toppe ***
rf_aux = rf; rf_aux.delay = Tt - rf_align_delay + 14e-6;

Gr = mr.makeExtendedTrapezoid( 'z','times', [0 TR-Tt], 'amplitudes', [Ag Ag]); % this is a constant graient with no ramps

% pre-ramp the gradient to Azc
% if abs(Azc) > eps
%     %Tpr=max(2,ceil(Azc/sys.maxSlew/sys.gradRasterTime))*sys.gradRasterTime;
%     Tpr = max( 2, ceil( Azc/sys.maxSlew/sys.blockDurationRaster)) * sys.blockDurationRaster;
%     assert( Tpr <= TR);
%     % this "dummy TR" does not have an actual RF pulse, which is not good
%     % when it is called for the inner shells... but otherwise there would be no enough spoiling...
%     seq.addBlock( rf_dummy, mr.makeDelay( rf_align_delay), mr.makeLabel( 'SET', 'TRID', trid_n));  % label needed to mark TOPPE segment, use actual rf for steady state between zte and spi
%     seq.addBlock( mr.align( 'right', mr.makeDelay( TR), mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, Tpr], 'amplitudes', [0, Azc])));
% end

% the loop itself

for j = 1:im
    Glast = struct( 'x', 0, 'y', 0, 'z', 0);
    Gcr = mr.rotate( 'z', phi(j), mr.rotate( 'y', theta(j), Gr));
    Gcurr = struct( 'x', 0, 'y', 0, 'z', 0);
    
    for g = 1:length( Gcr)
        Gcurr.(Gcr{g}.channel) = Gcr{g}.waveform(1);
    end

    % ramp from zero
    seq.addBlock( ...
        mr.makeExtendedTrapezoid( 'x', 'system', sys, 'times', [0, Tt + rf_align_delay], 'amplitudes', [0, 0]), ...
        mr.makeExtendedTrapezoid( 'y', 'system', sys, 'times', [0, Tt + rf_align_delay], 'amplitudes', [0, 0]), ...
        mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, Tt + rf_align_delay], 'amplitudes', [0, 0]), ...
        rf_aux, mr.makeLabel( 'SET', 'TRID', trid_n + j)); % label start of TRID segment
    seq.addBlock( ...
        mr.makeExtendedTrapezoid( 'x', 'system', sys, 'times', [0, TR - Tt - rf_align_delay], 'amplitudes', [0, Gcurr.x]), ...
        mr.makeExtendedTrapezoid( 'y', 'system', sys, 'times', [0, TR - Tt - rf_align_delay], 'amplitudes', [0, Gcurr.y]), ...
        mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, TR - Tt - rf_align_delay], 'amplitudes', [0, Gcurr.z]));

    Glast = Gcurr;

    % acquisition segment
    for i = j:im:Ns
        Gcr = mr.rotate( 'z', phi(i), mr.rotate( 'y', theta(i), Gr));
        Gcurr = struct( 'x', 0, 'y', 0, 'z', 0);
        
        for g = 1:length( Gcr)
            Gcurr.(Gcr{g}.channel) = Gcr{g}.waveform(1);
        end

        seq.addBlock( ...
            mr.makeExtendedTrapezoid( 'x', 'system', sys, 'times', [0, Tt - rf_align_delay], 'amplitudes', [Glast.x, Gcurr.x]), ...
            mr.makeExtendedTrapezoid( 'y', 'system', sys, 'times', [0, Tt - rf_align_delay], 'amplitudes', [Glast.y, Gcurr.y]), ...
            mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, Tt - rf_align_delay], 'amplitudes', [Glast.z, Gcurr.z]));
        seq.addBlock( ...
            mr.makeExtendedTrapezoid( 'x', 'system', sys, 'times', [0, rf_align_delay], 'amplitudes', [Gcurr.x, Gcurr.x]), ...
            mr.makeExtendedTrapezoid( 'y', 'system', sys, 'times', [0, rf_align_delay], 'amplitudes', [Gcurr.y, Gcurr.y]), ...
            mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, rf_align_delay], 'amplitudes', [Gcurr.z, Gcurr.z]), ...
            rf); % no label, first only
        seq.addBlock( [{adc}, Gcr]); % toppe does not allow rf and adc in same block

        Glast = Gcurr;
        shot = 2;
    end

    % ramp to zero
    seq.addBlock( ...
        mr.makeExtendedTrapezoid( 'x', 'system', sys, 'times', [0, Tt - rf_align_delay], 'amplitudes', [Glast.x, Glast.x]), ...
        mr.makeExtendedTrapezoid( 'y', 'system', sys, 'times', [0, Tt - rf_align_delay], 'amplitudes', [Glast.y, Glast.y]), ...
        mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, Tt - rf_align_delay], 'amplitudes', [Glast.z, Glast.z]));
    seq.addBlock( ...
        mr.makeExtendedTrapezoid( 'x', 'system', sys, 'times', [0, TR - Tt + rf_align_delay - ge_dt], 'amplitudes', [Glast.x, 0]), ...
        mr.makeExtendedTrapezoid( 'y', 'system', sys, 'times', [0, TR - Tt + rf_align_delay - ge_dt], 'amplitudes', [Glast.y, 0]), ...
        mr.makeExtendedTrapezoid( 'z', 'system', sys, 'times', [0, TR - Tt + rf_align_delay - ge_dt], 'amplitudes', [Glast.z, 0]), ...
        rf); % to maintain ss
end
end % function populate_subsequence( sys, seq, rf, adc, phi, theta, im, Ns, Ag, TR, Tt, FN, FO) ***********************

