

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
        MCS_Table
    end

    methods
        function obj = PuschChannel()

            % Table 5.1.3.1-1: MCS index table 1 for PDSCH
            % MCS Table 1 (QPSK, 16QAM, 64QAM) for PUSCH
            obj.MCS_Table.Modulation = {...
                'QPSK',  'QPSK',  'QPSK',  'QPSK',  'QPSK',  'QPSK',  'QPSK',  'QPSK',  'QPSK',  'QPSK', ...
                '16QAM', '16QAM', '16QAM', '16QAM', '16QAM', '16QAM', '16QAM', ...
                '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', '64QAM', ...
                'reserved', 'reserved', 'reserved'};
            obj.MCS_Table.TargetCodeRate = [...
                120, 157, 193, 251, 308, 379, 449, 526, 602, 679, ...
                340, 378, 434, 490, 553, 616, 658, ...
                438, 466, 517, 567, 616, 666, 719, 772, 822, 873, 910, 948, ...
                0, 0, 0] / 1024;

            obj.simParameters = struct();
            obj.simParameters.SNR = 40;
            obj.simParameters.PerfectChannelEstimator = true;
            obj.simParameters.DisplaySimulationInformation = true;

            % Set waveform type and PUSCH numerology (SCS and CP type)
            obj.simParameters.Carrier = nrCarrierConfig;        % Carrier resource grid configuration
            obj.simParameters.Carrier.NSizeGrid = 273;           % Bandwidth in number of resource blocks (52 RBs at 15 kHz SCS for 10 MHz BW)
            obj.simParameters.Carrier.SubcarrierSpacing = 30;   % 15, 30, 60, 120 (kHz)
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

            % Define codeword modulation and target code rate
            mcsIndex = 28; % MCS index (0-28)
            obj.simParameters.PUSCH.Modulation = obj.MCS_Table.Modulation{mcsIndex + 1}; % +1 for 1-based indexing
            obj.simParameters.PUSCHExtension.TargetCodeRate = obj.MCS_Table.TargetCodeRate(mcsIndex + 1);

            init_channel(obj);
            init_encodeULSCH(obj);
            init_decodeULSCH(obj);
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
            obj.encodeULSCH.TargetCodeRate = obj.simParameters.PUSCHExtension.TargetCodeRate;
        end

        function init_decodeULSCH(obj)
            obj.decodeULSCH = nrULSCHDecoder;
            obj.decodeULSCH.MultipleHARQProcesses = false;
            obj.decodeULSCH.TargetCodeRate = obj.simParameters.PUSCHExtension.TargetCodeRate;
            obj.decodeULSCH.LDPCDecodingAlgorithm = obj.simParameters.PUSCHExtension.LDPCDecodingAlgorithm;
            obj.decodeULSCH.MaximumLDPCIterationCount = obj.simParameters.PUSCHExtension.MaximumLDPCIterationCount;
        end
        
        function [rx, time_s] = tranceiver(obj, tx)
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

            decodeULSCHLocal = obj.decodeULSCH;  % Copy of the decoder handle to help PCT classification of variable
            decodeULSCHLocal.reset();        % Reset decoder at the start of each SNR point

            puschNonCodebook = pusch;
            puschNonCodebook.TransmissionScheme = 'nonCodebook';

            % Specify the fixed order in which we cycle through the HARQ process IDs
            harqSequence = 0:obj.simParameters.PUSCHExtension.NHARQProcesses-1; 
            rvSeq = [0 2 3 1];
            harqEntity = HARQEntity(harqSequence,rvSeq);

            rx = [];
            startIdx = 1;
            time_s = 0;
            nslot = 0;

            while startIdx < length(tx)
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
                    % Extract the transport block from the transmitted bits
                    endIdx = min(startIdx + trBlkSize - 1, length(tx));
                    trx_size = endIdx - startIdx + 1;
                    trBlk = tx(startIdx:endIdx);
                    fprintf("Transmitting %d bits\n", trx_size);

                    % Pad the transport block with zeros if needed
                    if trx_size < trBlkSize
                        trBlk = [trBlk; zeros(trBlkSize - trx_size, 1)];
                    end

                    setTransportBlock(obj.encodeULSCH, trBlk, harqEntity.HARQProcessID);
                    if harqEntity.SequenceTimeout
                        resetSoftBuffer(obj.decodeULSCHLocal, harqEntity.HARQProcessID);
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
                SNRdB = obj.simParameters.SNR;
                SNR = 10^(SNRdB/10);
                N0 = 1/sqrt(obj.simParameters.NRxAnts*double(obj.waveformInfo.Nfft)*SNR);
                noise = N0*randn(size(rxWaveform),"like",rxWaveform);
                rxWaveform = rxWaveform + noise;


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
                    [estChannelGrid,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsLayerIndices,dmrsLayerSymbols,'CDMLengths',pusch.DMRS.CDMLengths);
                end

                % Extract the PUSCH symbols from the received resource grid
                [~, puschAntIndices] = nrExtractResources(obj.puschIndices, rxGrid);

                % Get PUSCH resource elements from the received grid
                [puschRx, puschHest] = nrExtractResources(obj.puschIndices, rxGrid, estChannelGrid);

                % Equalization
                [puschEq, csi] = nrEqualizeMMSE(puschRx, puschHest, noiseEst);

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

                % Update current process with CRC error and advance to next process
                procstatus = updateAndAdvance(harqEntity,blkerr,trBlkSize,obj.puschIndicesInfo.G);
                if (obj.simParameters.DisplaySimulationInformation)
                    fprintf("NSlot: %d, Procstatus: %s\n", nslot, procstatus);
                end

                % Add the transport block to the received bits
                if (~blkerr)
                    % Update the start index for the next transmission
                    startIdx = startIdx + trBlkSize;
                    
                    if (trx_size < trBlkSize)
                        rx = [rx; decbits(1:trx_size)];
                    else
                        rx = [rx; decbits];
                    end
                else
                    fprintf("Block error\n");
                end
                time_s = time_s + obj.simParameters.slotDurationS;
                % 7D1S2U
                if mod(nslot, 2) == 0
                    time_s = time_s + 8 * obj.simParameters.slotDurationS;
                end
            end
        end
    end
end

