classdef InterferenceChirp
    %INTERFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SIR;
    end
    
    methods
        function obj = InterferenceChirp(SIR)
            obj.SIR = SIR;
        end

        function interference = getInterference(obj, rxWaveform, SampleRate)
            P_signal = mean(abs(rxWaveform).^2);
            P_interference = P_signal / 10^(obj.SIR/10);
            A_interference = sqrt(P_interference);

            t = (0:length(rxWaveform)-1)' / SampleRate;
            interference = A_interference * chirp(t, 1e6, 1e-3, 2e6);
        end
    end
end
