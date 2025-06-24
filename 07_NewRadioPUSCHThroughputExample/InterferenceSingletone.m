classdef InterferenceSingletone
    %INTERFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SIR;
    end
    
    methods
        function obj = InterferenceSingletone(SIR)
            obj.SIR = SIR;
        end

        function toneInterference = getInterference(obj, rxWaveform, SampleRate)
            P_signal = mean(abs(rxWaveform).^2);
            P_noise = P_signal / 10^(obj.SIR/10);
            toneAmp = sqrt(P_noise);   % Amplitude of the tone
            t = (0:length(rxWaveform)-1)' / SampleRate;
            toneFreq = 1e6;
            toneInterference = toneAmp * exp(1j*2*pi*toneFreq*t);
            return;
        end
    end
end

