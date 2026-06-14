classdef InterferenceFilterWGN
    %INTERFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SIR;
    end
    
    methods
        function obj = InterferenceFilterWGN(SIR)
            obj.SIR = SIR;
        end

        function interference = getInterference(obj, rxWaveform, SampleRate)
            P_signal = mean(abs(rxWaveform).^2);
            P_interference = P_signal / 10^(obj.SIR/10);
            A_interference = sqrt(P_interference);
            t = (0:length(rxWaveform)-1)' / SampleRate;
            interference = A_interference * randn(size(t));

            % Design a bandpass filter for 1e6 to 2e6 Hz
            filterOrder = 8;
            cutoffFreqs = [0.5e6, 1.5e6] / (SampleRate/2);  % Normalize to Nyquist frequency
            [b, a] = butter(filterOrder, cutoffFreqs, 'bandpass');

            % Apply the filter to the interference
            interference = filter(b, a, interference);
        end
    end
end

