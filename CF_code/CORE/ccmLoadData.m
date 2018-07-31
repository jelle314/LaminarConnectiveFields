function [data, params, coords] = ccmLoadData(view, params, roi)
% ccmLoadData - load and preprocess time series
%
% 2009: KVH adapted from rmLoadData.

data = [];

% place datasets behind each other. 
for ds = 1:numel(params.stim),
    [tSeries coords params] = ccmLoadDataROI(view, params, ds, roi);
	if isempty(data),
		dii.end = cumsum([params.stim(:).nFrames]./[params.stim(:).nUniqueRep]);
		dii.start = [1 dii.end(1:end-1)+1];
		data = zeros(dii.end(end), size(tSeries ,2));
	end;
	data(dii.start(ds):dii.end(ds),:) = tSeries;
end;

return;


function data=raw2pc(data)
    % convert to percent signal change
    dc = ones(size(data,1),1)*mean(data);
    data = ((data./dc) - 1) .*100;
return;


function [tSeries coords params] = ccmLoadDataROI(view, params, ds, r)
    % get ROI coords
    coords = view.ROIs(r).coords;

    % index into view's data
    [coordsIndex coords] = roiIndices(view, coords);

    % store roi info
    if r == view.selectedROI,
        params = ccmSet(params,'roiName',view.ROIs(r).name);
        params = ccmSet(params,'roiCoords',coords);
        params = ccmSet(params,'roiIndex',coordsIndex);
    end

    % process everything
    tSeries  = loadtSeries(view, ds);
    roiIndex = coordsIndex;
    tSeries = tSeries(:,roiIndex(:));  
    coords = roiIndex(:);

    % only convert to percent change if the flag is set
    if params.analysis.calcPC
        tSeries = raw2pc(tSeries);    
    end

    % apply low-pass filtering 
    if params.analysis.fc  
        tSeries = ccmLowPass(tSeries,params);
    end

return

