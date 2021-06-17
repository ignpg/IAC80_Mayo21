function VarAnalisys(fitspath, RAref, Decref)
	
	% Images read and extraction
	s = FITS.read2sim(fitspath);
	for idx = 1:numel(s)
	    s(idx) = mextractor(s(idx),'Thresh',20);
	end

	% Set mid-range image as reference while getting sources to correct by atmospheric effects
	sInit = s(round(end/2));

	% Sigmoid function to bestow contraste when showing image
	PixelData = reshape(sInit.Im,[1,numel(sInit.Im)]);
	BackMean = mean(PixelData(~isoutlier(PixelData)));
	BackStd = std(PixelData(~isoutlier(PixelData)));

	NormPixelData = 1./(1 + exp((-PixelData + (BackMean + 1*BackStd))./(2*BackStd)));
	imNorm = reshape(NormPixelData,size(sInit.Im));

	% Normalized-CrossCorrelation for each detection using PSF as kernel
	Surface = normxcorr2(sInit.PSF,sInit.Im);
	Surface(size(Surface,1)-9:size(Surface,1),:) = [];
	Surface(1:10,:) = [];
	Surface(:,size(Surface,2)-9:size(Surface,2)) = [];
	Surface(:,1:10) = [];

	Idx = sub2ind(size(Surface),round(sInit.Cat(:,sInit.Col.YPEAK_IMAGE)),round(sInit.Cat(:,sInit.Col.XPEAK_IMAGE)));
	CorrValue = Surface(Idx);
	% Delete those above correlation threshold
	CorrThresh = 0.7;
	UnderThresh = CorrValue < CorrThresh;
	sInit.Cat(UnderThresh,:) = [];

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WARNING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% UNCOMMENT IF IMAGE IS NOT CUTTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Calculate distance from image center                                                       %
	CenterDist = sqrt((size(sInit.Im,1)/2 - sInit.Cat(:,sInit.Col.XWIN_IMAGE)).^2 + (size(sInit.Im,2)/2 - sInit.Cat(:,sInit.Col.YWIN_IMAGE)).^2);
	DistThresh = 0.5*sqrt((size(sInit.Im,1)/2).^2 + (size(sInit.Im,2)/2).^2);                    %
	% And delete those woo far
	OverThresh = CenterDist > DistThresh;                                                        %
	sInit.Cat(OverThresh,:) = [];                                                                %
	%                                                                                            %
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	
	% Delete those above SNR threshold
	SNRThresh = 500;
	sInit.Cat(sInit.Cat(:,sInit.Col.SN_PSF) < SNRThresh,:) = [];

	% Plot reference detections over image
	figure(1)
	imshow(imNorm)
	hold on
	plot(sInit.Cat(:,sInit.Col.XWIN_IMAGE),sInit.Cat(:,sInit.Col.YWIN_IMAGE),'ro')

	% Get reference coordinates
	CoorRef = [sInit.Cat(:,sInit.Col.ALPHAWIN_J2000),sInit.Cat(:,sInit.Col.DELTAWIN_J2000)];

	% And delete target source from them
	[minDist,minIdx] = min(AngDistance([RAref Decref],CoorRef));
	CoorRef(minIdx,:) = [];

	% Plot target source
	plot(sInit.Cat(minIdx,sInit.Col.XWIN_IMAGE),sInit.Cat(minIdx,sInit.Col.YWIN_IMAGE),'go')

	% Initialize missing-meassures array
	NaNArray = zeros([size(CoorRef,1) 1]);

	% For each image...
	for idx = 1:numel(s)
	    % ...get datetime...
	    dateTime(idx) = datetime(s(idx).Header{8,2},'InputFormat','uuuu-MM-dd''T''HH:mm:ss');
	    % ...and for each reference star...
	    for refidx = 1:size(CoorRef,1)
	        % ...find nearest...
		[minDist,minIdx] = min(AngDistance(CoorRef(refidx,:),[s(idx).Cat(:,s(idx).Col.ALPHAWIN_J2000) s(idx).Cat(:,s(idx).Col.DELTAWIN_J2000)]));
		% ...and check if distance is below 10 arcsec...
		if minDist < 10
		    % ...get flux and flux error
		    RefFlux(idx,refidx) = s(idx).Cat(minIdx,s(idx).Col.FLUX_PSF);
		    RefFluxErr(idx,refidx) = s(idx).Cat(minIdx,s(idx).Col.FLUXERR_PSF);
		else
		    NaNArray(refidx) = 1;
		end
	    end
	end
	% Delete those reference detecions with missing meassures
	RefFlux(:,find(NaNArray))=[];

	% Reference flux normalization
	MeanFlux = nanmean(RefFlux);
	NormRefFlux = RefFlux./MeanFlux;

	AtmosphericCalibration = mean(NormRefFlux,2);


	% Target source analysis
	for idx = 1:numel(s)
[minDist,minIdx] = min(AngDistance([RAref Decref],[s(idx).Cat(:,s(idx).Col.ALPHAWIN_J2000) s(idx).Cat(:,s(idx).Col.DELTAWIN_J2000)]));
	    if minDist < 10
		VarFlux(idx) = s(idx).Cat(minIdx,s(idx).Col.FLUX_PSF);
		VarFluxErr(idx) = s(idx).Cat(minIdx,s(idx).Col.FLUXERR_PSF);
	    else
		VarFlux(idx) = NaN;
		VarFluxErr(idx) = NaN;
	    end
	end

	% Flux variation plot
	figure(2)
	errorbar(juliandate(dateTime),VarFlux'./AtmosphericCalibration,VarFluxErr'./AtmosphericCalibration)
	grid on
	xlabel('Julian day')
	ylabel('Flux')
end
