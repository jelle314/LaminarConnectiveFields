function tSeries = ccmLowPass(tSeries, params)

TR = params.stim(1).framePeriod;

for n = 1:numel(params.stim)
    if params.stim(n).framePeriod ~= TR;
        fprintf(1, ...
            '[%s]: Warning: Scans have different TRs. Assuming TR = %d sec for all scans.\n', ...
            mfilename,TR);
    end
end

fc = params.analysis.fc; % cutoff frequency (Hz) 
fs = 1/TR; % sampling frequency (Hz)

[b, a] = butter(10,fc/(fs/2),'low'); 
tSeries = filtfilt(b,a,tSeries);

end