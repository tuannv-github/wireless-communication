classdef InterferenceModulatedSignal
    %INTERFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SIR;
        carrier;
        pusch;
        puschIndices;
        puschIndicesInfo;

        puschmcsTables = nrPUSCHMCSTables;
        harqEntity;
        encodeULSCH;
        channel;
    end
    
    methods
        function obj = InterferenceModulatedSignal(SIR)
            if nargin < 1
                SIR = 15;
            end
            obj.SIR = SIR;

            obj.carrier = nrCarrierConfig;        % Carrier resource grid configuration
            obj.carrier.NSizeGrid = 25;           % Bandwidth in number of resource blocks (25 RBs at 15 kHz SCS for 5 MHz BW)
            obj.carrier.SubcarrierSpacing = 15;   % 15, 30, 60, 120 (kHz)
            obj.carrier.CyclicPrefix = 'Normal';  % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
            obj.carrier.NCellID = 0;              % Cell identity

            % PUSCH/UL-SCH parameters
            obj.pusch = nrPUSCHConfig;      % This PUSCH definition is the basis for all PUSCH transmissions in the BLER simulation
            % Define PUSCH time-frequency resource allocation per slot to be full grid (single full grid BWP)
            obj.pusch.PRBSet =  0:obj.carrier.NSizeGrid-1; % PUSCH PRB allocation
            obj.pusch.SymbolAllocation = [0,obj.carrier.SymbolsPerSlot]; % PUSCH symbol allocation in each slot
            obj.pusch.MappingType = 'A'; % PUSCH mapping type ('A'(slot-wise),'B'(non slot-wise))
            % Scrambling identifiers
            obj.pusch.NID = obj.carrier.NCellID;
            obj.pusch.RNTI = 1;
            % Define the transform precoding enabling, layering and transmission scheme
            obj.pusch.TransformPrecoding = false; % Enable/disable transform precoding
            obj.pusch.NumLayers = 1;              % Number of PUSCH transmission layers
            obj.pusch.TransmissionScheme = 'nonCodebook'; % Transmission scheme ('nonCodebook','codebook')
            obj.pusch.NumAntennaPorts = 1;        % Number of antenna ports for codebook based precoding
            obj.pusch.TPMI = 0;                   % Precoding matrix indicator for codebook based precoding
            % PUSCH DM-RS configuration
            obj.pusch.DMRS.DMRSTypeAPosition = 2;       % Mapping type A only. First DM-RS symbol position (2,3)
            obj.pusch.DMRS.DMRSLength = 1;              % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
            obj.pusch.DMRS.DMRSAdditionalPosition = 1;  % Additional DM-RS symbol positions (max range 0...3)
            obj.pusch.DMRS.DMRSConfigurationType = 1;   % DM-RS configuration type (1,2)
            obj.pusch.DMRS.NumCDMGroupsWithoutData = 2; % Number of CDM groups without data
            obj.pusch.DMRS.NIDNSCID = 0;                % Scrambling identity (0...65535)
            obj.pusch.DMRS.NSCID = 0;                   % Scrambling initialization (0,1)
            obj.pusch.DMRS.NRSID = 0;                   % Scrambling ID for low-PAPR sequences (0...1007)
            obj.pusch.DMRS.GroupHopping = 0;            % Group hopping (0,1)
            obj.pusch.DMRS.SequenceHopping = 0;         % Sequence hopping (0,1)

            obj.encodeULSCH = nrULSCH;
            obj.encodeULSCH.MultipleHARQProcesses = false;

            obj.channel = nrTDLChannel; % TDL channel object
            % Swap transmit and receive sides as the default TDL channel is
            % configured for downlink transmissions
            swapTransmitAndReceive(obj.channel);
            % Set the channel geometry
            obj.channel.NumTransmitAntennas = 1;
            obj.channel.NumReceiveAntennas = 1;
            % Assign simulation channel parameters and waveform sample rate to the object
            obj.channel.DelayProfile = 'TDL-A';
            obj.channel.DelaySpread = 30e-9;
            obj.channel.MaximumDopplerShift = 10;
            obj.channel.SampleRate = nrOFDMInfo(obj.carrier).SampleRate;
        end

        function rxWaveform = transmit(obj, modulation)
            XOverhead = 0;
            obj.pusch.Modulation = modulation;
            targetCodeRate = 0.8887;

            [obj.puschIndices, obj.puschIndicesInfo] = nrPUSCHIndices(obj.carrier, obj.pusch);
            trBlkSize = nrTBS(obj.pusch.Modulation, obj.pusch.NumLayers, numel(obj.puschIndicesInfo.PRBSet), ...
                                obj.puschIndicesInfo.NREPerPRB, targetCodeRate, XOverhead);
            
            trBlk = randi([0 1], trBlkSize, 1);

            setTransportBlock(obj.encodeULSCH, trBlk, 1);
            codedTrBlock = obj.encodeULSCH(obj.pusch.Modulation, obj.pusch.NumLayers, obj.puschIndicesInfo.G, 0);

            % Create resource grid for a slot
            puschGrid = nrResourceGrid(obj.carrier, 1);

            % PUSCH modulation, including codebook based MIMO precoding if TxScheme = 'codebook'
            puschSymbols = nrPUSCH(obj.carrier, obj.pusch, codedTrBlock);

            [~,puschAntIndices] = nrExtractResources(obj.puschIndices,puschGrid);
            puschGrid(puschAntIndices) = puschSymbols;

            dmrsSymbols = nrPUSCHDMRS(obj.carrier, obj.pusch);
            dmrsIndices = nrPUSCHDMRSIndices(obj.carrier, obj.pusch);
            
            F = eye(obj.pusch.NumLayers, 1);
            for p = 1:size(dmrsSymbols,2)
                [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),puschGrid);
                puschGrid(dmrsAntIndices) = puschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * F(p,:);
            end
            
            txWaveform = nrOFDMModulate(obj.carrier, puschGrid);
            txWaveform = [txWaveform; zeros(info(obj.channel).MaximumChannelDelay, size(txWaveform,2))];
            [rxWaveform, pathGains, sampleTimes] = obj.channel(txWaveform);
        end

        function interference = getInterference(obj, rxWaveform, SampleRate)
            P_signal = mean(abs(rxWaveform).^2);
            P_noise = P_signal / 10^(obj.SIR/10);
            A_noise = sqrt(P_noise);   % Amplitude of the tone

            interference_waveform = transmit(obj, 'QPSK');
            P_current_interference = mean(abs(interference_waveform).^2);
            A_current_interference = sqrt(P_current_interference);
            interference = interference_waveform / A_current_interference * A_noise;
        end
    end
end

