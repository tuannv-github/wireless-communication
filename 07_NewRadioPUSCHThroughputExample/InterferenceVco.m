classdef InterferenceVco
    %INTERFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SIR;
    end
    
    methods
        function obj = InterferenceVco(SIR)
            obj.SIR = SIR;
        end

        function interference = getInterference(obj, rxWaveform, SampleRate)

            t = (0:length(rxWaveform)-1)' / SampleRate;

            P_signal = mean(abs(rxWaveform).^2);
            P_noise = P_signal / 10^(obj.SIR/10);
            amp = sqrt(P_noise);   % Amplitude of the tone

            interference = amp*vco(chirp(t, 0, t(end), 10e3), [1e6 2e6], SampleRate);
            return;
        end
    end
end

