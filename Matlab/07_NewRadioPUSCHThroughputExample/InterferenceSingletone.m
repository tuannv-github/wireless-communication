classdef InterferenceSingletone
    %INTERFERENCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        SIR;
        toneFreq;
    end
    
    methods
        function obj = InterferenceSingletone(SIR)
            obj.SIR = SIR;
            obj.toneFreq = 1e6;
        end

        function toneInterference = getInterference(obj, rxWaveform, SampleRate)
            P_signal = mean(abs(rxWaveform).^2);
            P_noise = P_signal / 10^(obj.SIR/10);
            toneAmp = sqrt(P_noise);   % Amplitude of the tone
            t = (0:length(rxWaveform)-1)' / SampleRate;
            toneInterference = toneAmp * exp(1j*2*pi*obj.toneFreq*t);
            return;
        end
    end
end

