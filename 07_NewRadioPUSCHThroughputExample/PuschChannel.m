classdef PuschChannel < handle
    %PUSCHCHANNEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        simParameters
        puschIndices
        puschIndicesInfo
        channel
        maxChDelay
        encodeULSCH
        decodeULSCH
        waveformInfo
        constellationDiagram
        puschmcsTables = nrPUSCHMCSTables;
        MCSIndex;
        colorLimits = [0 2];
        grids;
    end

    methods
        function obj = PuschChannel(plot)
            % obj.constellationDiagram = comm.ConstellationDiagram;
            % obj.constellationDiagram.EnableMeasurements = true;

            obj.simParameters = struct();

            if nargin > 0
                obj.simParameters.Plot = plot;
            else
                obj.simParameters.Plot = false;
            end

            obj.simParameters.SNR = 45;
            obj.simParameters.MCSIndex = 28; % MCS index (0-28)

            obj.simParameters.PerfectChannelEstimator = false;
            obj.simParameters.DisplaySimulationInformation = true;

            if obj.simParameters.Plot
                % Create rxwaveform directory if it doesn't exist
                if ~exist('00_rxwaveform', 'dir')
                    mkdir('00_rxwaveform');
                end
                % Delete existing files in rxwaveform directory
                if exist('00_rxwaveform', 'dir')
                    delete('00_rxwaveform/*');
                end

                % Create rxgrid directory if it doesn't exist
                if ~exist('01_rxgrid', 'dir')
                    mkdir('01_rxgrid');
                end            
                % Delete existing files in rxgrid directory
                if exist('01_rxgrid', 'dir')
                    delete('01_rxgrid/*');
                end

                % Create constellation directory if it doesn't exist
                if ~exist('02_constellation', 'dir')
                    mkdir('02_constellation');
                end
                % Delete existing files in constellation directory
                if exist('02_constellation', 'dir')
                    delete('02_constellation/*');
                end

                % Create channel directory if it doesn't exist
                if ~exist('03_ulchannel', 'dir')
                    mkdir('03_ulchannel');
                end
                % Delete existing files in channel directory
                if exist('03_ulchannel', 'dir')
                    delete('03_ulchannel/*');
                end

                % Create tx directory if it doesn't exist
                if ~exist('04_txgrid', 'dir')
                    mkdir('04_txgrid');
                end
                % Delete existing files in tx directory
                if exist('04_txgrid', 'dir')
                    delete('04_txgrid/*');
                end

                % Create interference directory if it doesn't exist
                if ~exist('05_interference', 'dir')
                    mkdir('05_interference');
                end
                % Delete existing files in interference directory
                if exist('05_interference', 'dir')
                    delete('05_interference/*');
                end
            end

            % Set waveform type and PUSCH numerology (SCS and CP type)
            obj.simParameters.Carrier = nrCarrierConfig;        % Carrier resource grid configuration
            obj.simParameters.Carrier.NSizeGrid = 25;           % Bandwidth in number of resource blocks (25 RBs at 15 kHz SCS for 5 MHz BW)
            obj.simParameters.Carrier.SubcarrierSpacing = 15;   % 15, 30, 60, 120 (kHz)
            obj.simParameters.Carrier.CyclicPrefix = 'Normal';  % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
            obj.simParameters.Carrier.NCellID = 0;              % Cell identity

            obj.simParameters.slotDurationS = 1e-3 / (obj.simParameters.Carrier.SubcarrierSpacing/15); % Slot duration in seconds based on SCS

            % PUSCH/UL-SCH parameters
            obj.simParameters.PUSCH = nrPUSCHConfig;      % This PUSCH definition is the basis for all PUSCH transmissions in the BLER simulation
            
            % Define PUSCH time-frequency resource allocation per slot to be full grid (single full grid BWP)
            obj.simParameters.PUSCH.PRBSet =  0:obj.simParameters.Carrier.NSizeGrid-1; % PUSCH PRB allocation
            obj.simParameters.PUSCH.SymbolAllocation = [0,obj.simParameters.Carrier.SymbolsPerSlot]; % PUSCH symbol allocation in each slot
            obj.simParameters.PUSCH.MappingType = 'A'; % PUSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
            
            % Scrambling identifiers
            obj.simParameters.PUSCH.NID = obj.simParameters.Carrier.NCellID;
            obj.simParameters.PUSCH.RNTI = 1;
            
            % Define the transform precoding enabling, layering and transmission scheme
            obj.simParameters.PUSCH.TransformPrecoding = false; % Enable/disable transform precoding
            obj.simParameters.PUSCH.NumLayers = 1;              % Number of PUSCH transmission layers
            obj.simParameters.PUSCH.TransmissionScheme = 'nonCodebook'; % Transmission scheme ('nonCodebook','codebook')
            obj.simParameters.PUSCH.NumAntennaPorts = 1;        % Number of antenna ports for codebook based precoding
            obj.simParameters.PUSCH.TPMI = 0;                   % Precoding matrix indicator for codebook based precoding
            
            % PUSCH DM-RS configuration
            obj.simParameters.PUSCH.DMRS.DMRSTypeAPosition = 2;       % Mapping type A only. First DM-RS symbol position (2,3)
            obj.simParameters.PUSCH.DMRS.DMRSLength = 1;              % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
            obj.simParameters.PUSCH.DMRS.DMRSAdditionalPosition = 1;  % Additional DM-RS symbol positions (max range 0...3)
            obj.simParameters.PUSCH.DMRS.DMRSConfigurationType = 1;   % DM-RS configuration type (1,2)
            obj.simParameters.PUSCH.DMRS.NumCDMGroupsWithoutData = 2; % Number of CDM groups without data
            obj.simParameters.PUSCH.DMRS.NIDNSCID = 0;                % Scrambling identity (0...65535)
            obj.simParameters.PUSCH.DMRS.NSCID = 0;                   % Scrambling initialization (0,1)
            obj.simParameters.PUSCH.DMRS.NRSID = 0;                   % Scrambling ID for low-PAPR sequences (0...1007)
            obj.simParameters.PUSCH.DMRS.GroupHopping = 0;            % Group hopping (0,1)
            obj.simParameters.PUSCH.DMRS.SequenceHopping = 0;         % Sequence hopping (0,1)
            
            % Additional simulation and UL-SCH related parameters
            obj.simParameters.PUSCHExtension = struct();  % This structure is to hold additional simulation parameters for the UL-SCH and PUSCH

            % HARQ process and rate matching/TBS parameters
            obj.simParameters.PUSCHExtension.XOverhead = 0;       % Set PUSCH rate matching overhead for TBS (Xoh)
            obj.simParameters.PUSCHExtension.NHARQProcesses = 1;  % Number of parallel HARQ processes to use
            obj.simParameters.PUSCHExtension.EnableHARQ = true;   % Enable retransmissions for each process, using RV sequence [0,2,3,1]
            
            % LDPC decoder parameters
            % Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
            obj.simParameters.PUSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
            obj.simParameters.PUSCHExtension.MaximumLDPCIterationCount = 6;
            
            % Define the overall transmission antenna geometry at end-points
            % If using a CDL propagation channel then the integer number of antenna elements is
            % turned into an antenna panel configured when the channel model object is created
            obj.simParameters.NTxAnts = 1; % Number of transmit antennas
            obj.simParameters.NRxAnts = 1; % Number of receive antennas
            
            % Define the general CDL/TDL propagation channel parameters
            obj.simParameters.DelayProfile = 'TDL-A'; % Use TDL-A model (Indoor hotspot model)
            obj.simParameters.DelaySpread = 30e-9;
            obj.simParameters.MaximumDopplerShift = 10;

            init_channel(obj);
            init_encodeULSCH(obj);
            init_decodeULSCH(obj);

            obj.grids =  struct('Slot', {}, 'Rx', {}, 'NoiseInterference', {}, 'Interference', {});
        end

        function init_channel(obj)
            obj.waveformInfo = nrOFDMInfo(obj.simParameters.Carrier); % Get information about the baseband waveform after OFDM modulation step

            obj.channel = nrTDLChannel; % TDL channel object

            % Swap transmit and receive sides as the default TDL channel is
            % configured for downlink transmissions
            swapTransmitAndReceive(obj.channel);

            % Set the channel geometry
            obj.channel.NumTransmitAntennas = obj.simParameters.NTxAnts;
            obj.channel.NumReceiveAntennas = obj.simParameters.NRxAnts;

            % Assign simulation channel parameters and waveform sample rate to the object
            obj.channel.DelayProfile = obj.simParameters.DelayProfile;
            obj.channel.DelaySpread = obj.simParameters.DelaySpread;
            obj.channel.MaximumDopplerShift = obj.simParameters.MaximumDopplerShift;
            obj.channel.SampleRate = obj.waveformInfo.SampleRate;

            % Get the maximum channel delay.
            chInfo = info(obj.channel);
            obj.maxChDelay = chInfo.MaximumChannelDelay;
        end

        function init_encodeULSCH(obj)
            obj.encodeULSCH = nrULSCH;
            obj.encodeULSCH.MultipleHARQProcesses = false;
        end

        function init_decodeULSCH(obj)
            obj.decodeULSCH = nrULSCHDecoder;
            obj.decodeULSCH.MultipleHARQProcesses = false;
            obj.decodeULSCH.LDPCDecodingAlgorithm = obj.simParameters.PUSCHExtension.LDPCDecodingAlgorithm;
            obj.decodeULSCH.MaximumLDPCIterationCount = obj.simParameters.PUSCHExtension.MaximumLDPCIterationCount;
        end
        
        function ret = tranceiver(obj, crc_count, interference)
            % Transceiver function that simulates PUSCH transmission and reception
            % Inputs:
            %   tx - Transmitted bits
            % Outputs:
            %   rx - Received bits
            %   time_s - Simulation time in seconds

            % TODO: Implement actual PUSCH transmission and reception chain
            % For now, just pass through the bits and return dummy timing

            carrier = obj.simParameters.Carrier;
            pusch = obj.simParameters.PUSCH;

            pusch.Modulation = obj.puschmcsTables.TransformPrecodingQAM64Table.Modulation{obj.simParameters.MCSIndex};
            obj.simParameters.PUSCHExtension.TargetCodeRate = obj.puschmcsTables.TransformPrecodingQAM64Table.TargetCodeRate(obj.simParameters.MCSIndex);

            obj.encodeULSCH.TargetCodeRate = obj.simParameters.PUSCHExtension.TargetCodeRate;
            obj.decodeULSCH.TargetCodeRate = obj.simParameters.PUSCHExtension.TargetCodeRate;

            decodeULSCHLocal = obj.decodeULSCH;  % Copy of the decoder handle to help PCT classification of variable
            decodeULSCHLocal.reset();        % Reset decoder at the start of each SNR point

            puschNonCodebook = pusch;
            puschNonCodebook.TransmissionScheme = 'nonCodebook';

            % Specify the fixed order in which we cycle through the HARQ process IDs
            harqSequence = 0:obj.simParameters.PUSCHExtension.NHARQProcesses-1; 
            rvSeq = [0 2 3 1];
            harqEntity = HARQEntity(harqSequence,rvSeq);

            offset = 0;

            time_s = 0;
            nslot = 0;

            timeout_counter = 0;

            while crc_count > 0 
                % Set the slot number
                carrier.NSlot = nslot;
                nslot = nslot + 1;

                % Calculate the PUSCH indices
                [obj.puschIndices, obj.puschIndicesInfo] = nrPUSCHIndices(carrier, pusch);
                % Calculate the transport block size for the transmission in the slot
                MRB = numel(obj.puschIndicesInfo.PRBSet);
                trBlkSize = nrTBS(pusch.Modulation, pusch.NumLayers, MRB, obj.puschIndicesInfo.NREPerPRB, obj.simParameters.PUSCHExtension.TargetCodeRate, obj.simParameters.PUSCHExtension.XOverhead);

                % HARQ processing
                % If new data for current process then create a new UL-SCH transport block
                if harqEntity.NewData
                    fprintf("New data for HARQ process %d\n", harqEntity.HARQProcessID);
                    trBlk = randi([0 1], trBlkSize, 1);

                    setTransportBlock(obj.encodeULSCH, trBlk, harqEntity.HARQProcessID);
                    if harqEntity.SequenceTimeout
                        fprintf("Resetting soft buffer for HARQ process %d\n", harqEntity.HARQProcessID);
                        resetSoftBuffer(decodeULSCHLocal, harqEntity.HARQProcessID);
                        timeout_counter = timeout_counter + 1;
                        if timeout_counter > 10
                            fprintf("HARQ process %d timeout counter exceeded 10\n", harqEntity.HARQProcessID);
                            ret = false;
                            return;
                        end
                    else
                        timeout_counter = 0;
                    end
                end

                % Encode the UL-SCH transport block
                codedTrBlock = obj.encodeULSCH(pusch.Modulation, pusch.NumLayers, ...
                    obj.puschIndicesInfo.G, harqEntity.RedundancyVersion);

                % Create resource grid for a slot
                puschGrid = nrResourceGrid(carrier,obj.simParameters.NTxAnts);

                % PUSCH modulation, including codebook based MIMO precoding if TxScheme = 'codebook'
                puschSymbols = nrPUSCH(carrier,pusch,codedTrBlock);

                % Implementation-specific PUSCH MIMO precoding and mapping. This 
                % MIMO precoding step is in addition to any codebook based 
                % MIMO precoding done during PUSCH modulation above
                if (strcmpi(pusch.TransmissionScheme,'codebook'))
                    % Codebook based MIMO precoding, F precodes between PUSCH
                    % transmit antenna ports and transmit antennas
                    F = eye(pusch.NumAntennaPorts,obj.simParameters.NTxAnts);
                else
                    % Non-codebook based MIMO precoding, F precodes between PUSCH 
                    % layers and transmit antennas
                    F = eye(pusch.NumLayers,obj.simParameters.NTxAnts);
                end

                [~,puschAntIndices] = nrExtractResources(obj.puschIndices,puschGrid);
                puschGrid(puschAntIndices) = puschSymbols * F;

                % Implementation-specific PUSCH DM-RS MIMO precoding and mapping.
                % The first DM-RS creation includes codebook based MIMO precoding if applicable
                dmrsSymbols = nrPUSCHDMRS(carrier,pusch);
                dmrsIndices = nrPUSCHDMRSIndices(carrier,pusch);

                for p = 1:size(dmrsSymbols,2)
                    [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),puschGrid);
                    puschGrid(dmrsAntIndices) = puschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * F(p,:);
                end
                if obj.simParameters.Plot
                    h = figure('Visible', 'off');
                    imagesc(abs(puschGrid(:,:,1)));
                    clim(obj.colorLimits);
                    colorbar;
                    title('PUSCH Grid Magnitude');
                    saveas(h, sprintf('04_txgrid/puschgrid_slot_%d.png', nslot));
                    close(h);
                end
                if obj.simParameters.Plot
                    h = figure('Visible', 'off');
                    scatter(real(puschGrid(puschAntIndices)), imag(puschGrid(puschAntIndices)), '.');
                    grid on;
                    title(sprintf('PUSCH Tx Constellation Diagram - Slot %d', nslot));
                    xlabel('In-Phase');
                    ylabel('Quadrature');
                    axis equal;
                    saveas(h, sprintf('04_txgrid/pusch_tx_constellation_slot_%d.png', nslot));
                    close(h);
                end

                % OFDM modulation
                txWaveform = nrOFDMModulate(carrier, puschGrid);

                % Pass data through channel model. Append zeros at the end of the
                % transmitted waveform to flush channel content. These zeros take
                % into account any delay introduced in the channel. This is a mix
                % of multipath delay and implementation delay. This value may 
                % change depending on the sampling rate, delay profile and delay
                % spread
                txWaveform = [txWaveform; zeros(obj.maxChDelay,size(txWaveform,2))]; %#ok<AGROW>
                [rxWaveform,pathGains,sampleTimes] = obj.channel(txWaveform);

                % Add AWGN to the received time domain waveform 
                % Normalize noise power by the IFFT size used in OFDM modulation,
                % as the OFDM modulator applies this normalization to the
                % transmitted waveform. Also normalize by the number of receive
                % antennas, as the channel model applies this normalization to the
                % received waveform, by default
                fprintf('Required SNR: %.2f dB\n', obj.simParameters.SNR);
                SNRdB = obj.simParameters.SNR;
                SNR = 10^(SNRdB/10);
                N0 = 1/sqrt(obj.simParameters.NRxAnts*double(obj.waveformInfo.Nfft)*SNR);
                noise = N0*randn(size(rxWaveform),"like",rxWaveform);

                % Plot the received waveform magnitude and phase for each receive antenna
                if obj.simParameters.Plot
                    h = figure('Visible', 'off');
                    numRxAnts = size(rxWaveform, 2);
                    for ant = 1:numRxAnts
                        subplot(2*numRxAnts,1,2*ant-1);
                        plot(abs(rxWaveform(:,ant)));
                        title(sprintf('Received Waveform Magnitude - Slot %d, Rx Ant %d', nslot, ant));
                        xlabel('Time');
                        ylabel('Magnitude');

                        subplot(2*numRxAnts,1,2*ant);
                        spectrogram(rxWaveform(:,ant), 256, 240, 256, obj.waveformInfo.SampleRate, 'xaxis');
                        title(sprintf('Received Waveform Spectrogram - Slot %d, Rx Ant %d', nslot, ant));
                        colorbar;
                    end
                    saveas(h, sprintf('00_rxwaveform/rxwaveform_slot_%d.png', nslot));
                    close(h);
                end

                % Add interference
                if nargin > 2 && ~isempty(interference)
                    interferenceWaveform = interference.getInterference(rxWaveform, obj.waveformInfo.SampleRate);
                    if obj.simParameters.Plot
                        h = figure('Visible', 'off');
                        subplot(2,1,1);
                        plot(abs(fftshift(fft(interferenceWaveform))));
                        title(sprintf('Interference FFT Magnitude - Slot %d', nslot));
                        xlabel('Frequency');
                        ylabel('Magnitude');
                        
                        subplot(2,1,2);
                        spectrogram(interferenceWaveform, 256, 240, 256, obj.waveformInfo.SampleRate, 'xaxis');
                        title(sprintf('Interference Spectrogram - Slot %d', nslot));
                        colorbar;
                        saveas(h, sprintf('05_interference/interference_spectrogram_slot_%d.png', nslot));
                        close(h);
                    end
                else
                    interferenceWaveform = zeros(size(rxWaveform));
                end

                rxWaveform = rxWaveform + noise + interferenceWaveform;

                % Plot the received waveform magnitude and phase for each receive antenna
                if obj.simParameters.Plot
                    h = figure('Visible', 'off');
                    numRxAnts = size(rxWaveform, 2);
                    for ant = 1:numRxAnts
                        subplot(2*numRxAnts,1,2*ant-1);
                        plot(abs(rxWaveform(:,ant)));
                        title(sprintf('Received Waveform Magnitude - Slot %d, Rx Ant %d', nslot, ant));
                        xlabel('Time');
                        ylabel('Magnitude');

                        subplot(2*numRxAnts,1,2*ant);
                        spectrogram(rxWaveform(:,ant), 256, 240, 256, obj.waveformInfo.SampleRate, 'xaxis');
                        title(sprintf('Received Waveform Spectrogram - Slot %d, Rx Ant %d', nslot, ant));
                        colorbar;
                    end
                    saveas(h, sprintf('00_rxwaveform/rxwaveform_slot_%d.png', nslot));
                    close(h);
                end

                % % Plot the tone noise in time domain
                % h = figure('Visible', 'on');
                % t = (0:length(toneNoise)-1)' / obj.waveformInfo.SampleRate;  % Time vector in seconds
                
                % subplot(2,1,1);
                % plot(t, real(toneNoise));
                % title(sprintf('Tone Noise Real Component - Slot %d', nslot));
                % xlabel('Time (s)');
                % ylabel('Amplitude');
                % grid on;
                
                % subplot(2,1,2);
                % plot(t, imag(toneNoise));
                % title(sprintf('Tone Noise Imaginary Component - Slot %d', nslot));
                % xlabel('Time (s)');
                % ylabel('Amplitude');
                % grid on;
                % % saveas(h, sprintf('tonenoise/tonenoise_slot_%d.png', nslot));
                % % close(h);

                if (obj.simParameters.PerfectChannelEstimator)
                    % Perfect synchronization. Use information provided by the
                    % channel to find the strongest multipath component
                    pathFilters = getPathFilters(obj.channel);
                    [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
                else
                    % Practical synchronization. Correlate the received waveform 
                    % with the PUSCH DM-RS to give timing offset estimate 't' and
                    % correlation magnitude 'mag'. The function
                    % hSkipWeakTimingOffset is used to update the receiver timing
                    % offset. If the correlation peak in 'mag' is weak, the current
                    % timing estimate 't' is ignored and the previous estimate
                    % 'offset' is used
                    [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
                    offset = hSkipWeakTimingOffset(offset,t,mag);
                    % Display a warning if the estimated timing offset exceeds the
                    % maximum channel delay
                    if offset > obj.maxChDelay
                        warning(['Estimated timing offset (%d) is greater than the maximum channel delay (%d).' ...
                            ' This will result in a decoding failure. This may be caused by low SNR,' ...
                            ' or not enough DM-RS symbols to synchronize successfully.'],offset,obj.maxChDelay);
                    end
                end

                rxWaveform = rxWaveform(1+offset:end,:);

                % Perform OFDM demodulation on the received data to recreate the
                % resource grid, including padding in the event that practical
                % synchronization results in an incomplete slot being demodulated   
                rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
                [K,L,R] = size(rxGrid);
                if (L < carrier.SymbolsPerSlot)
                    rxGrid = cat(2, rxGrid, zeros(K,carrier.SymbolsPerSlot-L, R));
                end

                time_s = time_s + obj.simParameters.slotDurationS;
                % 7D1S2U
                if mod(nslot, 2) == 0
                    time_s = time_s + 8 * obj.simParameters.slotDurationS;
                end

                % h = figure('Visible', 'off');
                % imagesc(abs(rxGrid(:,:,1)));
                % colorbar;
                % title(sprintf('Received Resource Grid Magnitude - Slot %d', nslot));
                % xlabel('OFDM Symbols');
                % ylabel('Subcarriers');
                % saveas(h, sprintf('rxgrid/rxgrid_slot_%d.png', nslot));
                % close(h);

                if obj.simParameters.Plot
                    h = figure('Visible', 'off');
                    imagesc(abs(rxGrid(:,:,1)));
                    clim(obj.colorLimits);
                    set(gca, 'YDir', 'normal');
                    colorbar;
                    title(sprintf('Received Resource Grid Magnitude - Slot %d', nslot));
                    xlabel('OFDM Symbols');
                    ylabel('Subcarriers');

                    % Highlight DMRS symbols with border
                    hold on;
                    dmrsIndices = nrPUSCHDMRSIndices(carrier, pusch);
                    [dmrsRows, dmrsCols] = ind2sub(size(rxGrid(:,:,1)), dmrsIndices);
                    plot(dmrsCols, dmrsRows, 'r.', 'MarkerSize', 1);
                    % Add border around DMRS symbols
                    for i = 1:length(dmrsRows)
                        rectangle('Position', [dmrsCols(i)-0.5, dmrsRows(i)-0.5, 1, 1], ...
                                'EdgeColor', 'r', 'LineWidth', 0.5);
                    end
                    hold off;
                    
                    % Zoom the whole image
                    set(gca, 'FontSize', 20); % Increase font size
                    set(gcf, 'Position', [100, 100, 1200, 800]); % Make figure larger
                    
                    saveas(h, sprintf('01_rxgrid/rxgrid_slot_%d.png', nslot));
                    close(h);
                end

                if (obj.simParameters.PerfectChannelEstimator)

                    % Perfect channel estimation, use the value of the path gains
                    % provided by the channel
                    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);

                    % Get perfect noise estimate (from the noise realization)
                    noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end,:));
                    noiseEst = var(noiseGrid(:));

                    % Get perfect noise estimate (from the noise realization)
                    noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end,:));
                    noiseEst = var(noiseGrid(:));

                    % Apply MIMO deprecoding to estChannelGrid to give an estimate
                    % per transmission layer
                    K = size(estChannelGrid,1);
                    estChannelGrid = reshape(estChannelGrid,K*carrier.SymbolsPerSlot*obj.simParameters.NTxAnts,obj.simParameters.NRxAnts);
                    estChannelGrid = estChannelGrid * F.';
                    if (strcmpi(pusch.TransmissionScheme,'codebook'))
                        W = nrPUSCHCodebook(pusch.NumLayers,pusch.NumAntennaPorts,pusch.TPMI,pusch.TransformPrecoding);
                        estChannelGrid = estChannelGrid * W.';
                    end
                    estChannelGrid = reshape(estChannelGrid,K,carrier.SymbolsPerSlot,obj.simParameters.NTxAnts,[]);
                else
                    % Practical channel estimation between the received grid and
                    % each transmission layer, using the PUSCH DM-RS for each layer
                    % which are created by specifying the non-codebook transmission
                    % scheme
                    dmrsLayerSymbols = nrPUSCHDMRS(carrier,puschNonCodebook);
                    dmrsLayerIndices = nrPUSCHDMRSIndices(carrier,puschNonCodebook);
                    [estChannelGrid, noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsLayerIndices,dmrsLayerSymbols,'CDMLengths',pusch.DMRS.CDMLengths);
                    fprintf('Noise estimate dB: %.2f\n', 10*log10(noiseEst));

                    tx_dmrs = puschGrid(dmrsLayerIndices);
                    rx_dmrs = rxGrid(dmrsLayerIndices);
                    tx_dmrs_power = mean(abs(tx_dmrs).^2);
                    rx_dmrs_power = mean(abs(rx_dmrs).^2);
                    noise_power = mean(abs(rx_dmrs - estChannelGrid(dmrsLayerIndices).*tx_dmrs).^2);
                    snr = rx_dmrs_power / noise_power;
                    fprintf('Calculated SNRdB: %.2f\n', 10*log10(snr));

                    if obj.simParameters.Plot
                        receivedDMRS = rxGrid(dmrsLayerIndices);
                        h = figure('Visible', 'off');
                        scatter(real(receivedDMRS), imag(receivedDMRS), '.');
                        grid on;
                        title('DMRS Symbols Constellation');
                        xlabel('In-Phase (I)');
                        ylabel('Quadrature (Q)');
                        axis equal;
                        saveas(h, sprintf('02_constellation/rx_dmrs_constellation_slot_%d.png', nslot));
                        close(h);
                    end
                    if obj.simParameters.Plot
                        h = figure('Visible', 'off');
                        imagesc(abs(estChannelGrid(:,:,1)));
                        clim(obj.colorLimits);
                        colorbar;
                        title('Estimated Channel Grid Magnitude');
                        xlabel('OFDM Symbols');
                        ylabel('Subcarriers');
                        saveas(h, sprintf('03_ulchannel/channel_slot_%d.png', nslot));
                        close(h);
                    end
                end

                % Extract the PUSCH symbols from the received resource grid
                [~, puschAntIndices] = nrExtractResources(obj.puschIndices, rxGrid);

                % Get PUSCH resource elements from the received grid
                [puschRx, puschHest] = nrExtractResources(obj.puschIndices, rxGrid, estChannelGrid);

                if obj.simParameters.Plot
                    % Plot constellation diagram
                    h = figure('Visible', 'off');
                    scatter(real(puschRx), imag(puschRx), '.');
                    grid on;
                    title(sprintf('PUSCH Rx Constellation Diagram - Slot %d', nslot));
                    xlabel('In-Phase');
                    ylabel('Quadrature');
                    axis equal;
                    saveas(h, sprintf('02_constellation/pusch_rx_constellation_slot_%d.png', nslot));
                    close(h);
                end

                % Equalization
                [puschEq, csi] = nrEqualizeMMSE(puschRx, puschHest, noiseEst);

                % Plot constellation diagram using ConstellationDiagram                
                % obj.constellationDiagram.SamplesPerSymbol = 1;
                % obj.constellationDiagram.ShowTrajectory = false;
                % obj.constellationDiagram.ShowReferenceConstellation = true;
                % obj.constellationDiagram.ReferenceConstellation = obj.puschEq;
                % obj.constellationDiagram.Title = 'PUSCH Constellation Diagram';
                % obj.constellationDiagram.XLabel = 'In-Phase';
                % obj.constellationDiagram.YLabel = 'Quadrature';
                % obj.constellationDiagram(obj.puschEq);

                if obj.simParameters.Plot
                    % Plot constellation diagram
                    h = figure('Visible', 'off');
                    scatter(real(puschEq), imag(puschEq), '.');
                    grid on;
                    title(sprintf('PUSCH Constellation Diagram - Slot %d', nslot));
                    xlabel('In-Phase');
                    ylabel('Quadrature');
                    axis equal;
                    saveas(h, sprintf('02_constellation/pusch_constellation_slot_%d.png', nslot));
                    close(h);
                end

                % Decode PUSCH physical channel
                [ulschLLRs,rxSymbols] = nrPUSCHDecode(carrier, puschNonCodebook, puschEq, noiseEst);

                % Apply channel state information (CSI) produced by the equalizer,
                % including the effect of transform precoding if enabled
                if (pusch.TransformPrecoding)
                    MSC = MRB * 12;
                    csi = nrTransformDeprecode(csi,MRB) / sqrt(MSC);
                    csi = repmat(csi((1:MSC:end).'),1,MSC).';
                    csi = reshape(csi,size(rxSymbols));
                end
                csi = nrLayerDemap(csi);
                Qm = length(ulschLLRs) / length(rxSymbols);
                csi = reshape(repmat(csi{1}.',Qm,1),[],1);
                ulschLLRs = ulschLLRs .* csi;

                % Decode the UL-SCH transport channel
                decodeULSCHLocal.TransportBlockLength = trBlkSize;
                [decbits, blkerr] = decodeULSCHLocal(ulschLLRs, pusch.Modulation, pusch.NumLayers, harqEntity.RedundancyVersion);

                if (~blkerr)
                    NoiseInterferenceGrid = getNoiseInterference(decbits, obj.encodeULSCH, harqEntity, pusch, carrier, obj.simParameters.NTxAnts, rxGrid, estChannelGrid);

                    if (obj.simParameters.Plot)
                        % Split H into I and Q components
                        NoiseInterferenceGrid_I = real(NoiseInterferenceGrid);
                        NoiseInterferenceGrid_Q = imag(NoiseInterferenceGrid);
                        NoiseInterferenceGrid_zero = zeros(size(NoiseInterferenceGrid_I));
                        NoiseInterferenceGrid_IQ = cat(3, NoiseInterferenceGrid_I, NoiseInterferenceGrid_Q, NoiseInterferenceGrid_zero);
                        % Save H_IQ matrix as image without plotting
                        imwrite(NoiseInterferenceGrid_IQ, sprintf('03_ulchannel/noise_interference_slot_%d.png', nslot));
                        save(sprintf('03_ulchannel/noise_interference_slot_%d.mat', nslot), 'NoiseInterferenceGrid_IQ');
                    end
                    if (obj.simParameters.Plot)
                        h = figure('Visible', 'off');
                        imagesc(abs(NoiseInterferenceGrid(:,:,1)));
                        clim(obj.colorLimits);
                        colorbar;
                        title('Noise Interference Grid Magnitude');
                        xlabel('OFDM Symbols');
                        ylabel('Subcarriers');
                        saveas(h, sprintf('03_ulchannel/noise_interference_rgb_slot_%d.png', nslot));
                        close(h);
                    end

                    if nargin > 2 && ~isempty(interference)
                        InterferenceGrid = nrOFDMDemodulate(carrier, interferenceWaveform(1+offset:end,:));

                        if obj.simParameters.Plot
                            InterferenceGrid_abs = abs(InterferenceGrid(:,:,1));
                            imwrite(mat2gray(InterferenceGrid_abs), sprintf('03_ulchannel/interference_slot_%d.png', nslot));
                            save(sprintf('03_ulchannel/interference_slot_%d.mat', nslot), 'InterferenceGrid_abs');
                        end
                        if obj.simParameters.Plot
                            h = figure('Visible', 'off');
                            imagesc(abs(InterferenceGrid(:,:,1)));
                            clim(obj.colorLimits);
                            colorbar;
                            title('Interference Grid Magnitude');
                            xlabel('OFDM Symbols');
                            ylabel('Subcarriers');
                            saveas(h, sprintf('03_ulchannel/interference_rgb_slot_%d.png', nslot));
                            close(h);
                        end
                    else
                        InterferenceGrid = [];
                    end

                    row = struct('Slot', nslot, 'Rx', rxGrid, 'NoiseInterference', NoiseInterferenceGrid, 'Interference', InterferenceGrid);
                    obj.grids = [obj.grids; row];

                    crc_count = crc_count - 1;
                end

                % Update current process with CRC error and advance to next process
                procstatus = updateAndAdvance(harqEntity,blkerr,trBlkSize,obj.puschIndicesInfo.G);
                if (obj.simParameters.DisplaySimulationInformation)
                    fprintf("NSlot: %d, Procstatus: %s\n", nslot, procstatus);
                end

                % Add the transport block to the received bits
                if (blkerr)
                    fprintf("Block error\n");
                    continue
                end
            end

            ret = true;
        end
    end
end

