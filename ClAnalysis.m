function ClAnalysis(Bfitspath,Vfitspath,Ifitspath,arcsecThresh,showBool)
    % Image analysis
    sB = mextractor(FITS.read2sim(Bfitspath),'Thresh',5,'GAIN',3.78);
    sB.Cat(isnan(sB.Cat(:,sB.Col.MAG_PSF)),:)=[];
    sV = mextractor(FITS.read2sim(Vfitspath),'Thresh',5,'GAIN',3.78);
    sV.Cat(isnan(sV.Cat(:,sV.Col.MAG_PSF)),:)=[];
    sI = mextractor(FITS.read2sim(Ifitspath),'Thresh',5,'GAIN',3.78);
    sI.Cat(isnan(sI.Cat(:,sI.Col.MAG_PSF)),:)=[];
    
    % Load catalog from simbad website
    CatData = SimbadLoad(sV);
    CatData(CatData(:,4)==999 | CatData(:,5)==999 | CatData(:,7)==999,:) = [];
    
    % Initialize arrays used for atmospheric and standar calibration
    bIdx = []; vIdx = []; iIdx = []; bM_ins = []; bM_cat = []; vM_ins = []; vM_cat = []; iM_ins = [];  iM_cat = [];
    % For each detection in V image...
    for i=1:size(sV.Cat,1)
        % ... find same detection in B, I images and catalog reference
        [minDistB,minIdxB] = min(AngDistance([sV.Cat(i,sV.Col.ALPHAWIN_J2000),sV.Cat(i,sV.Col.DELTAWIN_J2000)],...
                                             [sB.Cat(:,sB.Col.ALPHAWIN_J2000),sB.Cat(:,sB.Col.DELTAWIN_J2000)]));
        [minDistI,minIdxI] = min(AngDistance([sV.Cat(i,sV.Col.ALPHAWIN_J2000),sV.Cat(i,sV.Col.DELTAWIN_J2000)],...
                                             [sI.Cat(:,sI.Col.ALPHAWIN_J2000),sI.Cat(:,sI.Col.DELTAWIN_J2000)]));
        [minDistC,minIdxC] = min(AngDistance([sV.Cat(i,sV.Col.ALPHAWIN_J2000),sV.Cat(i,sV.Col.DELTAWIN_J2000)],...
                                             [CatData(:,1),                   CatData(:,2)]));
        if minDistB < arcsecThresh && minDistI < arcsecThresh
            bIdx = [bIdx ; minIdxB];
            vIdx = [vIdx ; i];
            iIdx = [iIdx ; minIdxI];
        end
        if minDistB < arcsecThresh && minDistI < arcsecThresh && minDistC < arcsecThresh
            bM_ins = [bM_ins ; sB.Cat(minIdxB,sB.Col.MAG_PSF)];
            bM_cat = [bM_cat ; CatData(minIdxC,4)];
            vM_ins = [vM_ins ; sV.Cat(i,sV.Col.MAG_PSF)];
            vM_cat = [vM_cat ; CatData(minIdxC,5)];
            iM_ins = [iM_ins ; sI.Cat(minIdxI,sI.Col.MAG_PSF)];
            iM_cat = [iM_cat ; CatData(minIdxC,7)];
        end
    end
    
    % Discard missmatched detections
    outliersMatched = isoutlier(bM_ins - bM_cat) | isoutlier(vM_ins - vM_cat) | isoutlier(iM_ins - iM_cat);
    bM_ins(outliersMatched) = []; bM_cat(outliersMatched) = []; 
    vM_ins(outliersMatched) = []; vM_cat(outliersMatched) = []; 
    iM_ins(outliersMatched) = []; iM_cat(outliersMatched) = [];
    
    % Calculate zero point for each filter
    fb = fit(bM_ins,bM_cat,'poly1');
    if showBool
        figure(1)
        plot(bM_ins,bM_cat,'k.')
        hold on
        plot(fb,'b')
        grid on
        xlabel('M_i_n_s')
        ylabel('M_c_a_t')
        title('Filter B')
    end
    bZP = fb.p2;
    fv = fit(vM_ins,vM_cat,'poly1');
    if showBool
        figure(2)
        plot(vM_ins,vM_cat,'k.')
        hold on
        plot(fv,'g')
        grid on
        xlabel('M_i_n_s')
        ylabel('M_c_a_t')
        title('Filter V')
    end
    vZP = fv.p2;
    fi = fit(iM_ins,iM_cat,'poly1');
    if showBool
        figure(3)
        plot(iM_ins,iM_cat,'k.')
        hold on
        plot(fi,'r')
        grid on
        xlabel('M_i_n_s')
        ylabel('M_c_a_t')
        title('Filter I')
    end
    iZP = fi.p2;
    
    % And get C1, C2, C3 and C4 by color calibration
    cal_Bobs = bM_ins + bZP;
    cal_Vobs = vM_ins + vZP;
    cal_BmVobs = cal_Bobs - cal_Vobs;
    cal_BmV = bM_cat - vM_cat;
    cal_Iobs = iM_ins + iZP;
    cal_VmIobs = cal_Vobs - cal_Iobs;
    cal_VmI = vM_cat - iM_cat;
    
    bvfstd12 = fit(cal_BmV,vM_cat-cal_Vobs,'poly1');
    bvC1 = bvfstd12.p1;
    bvC2 = bvfstd12.p2;
    bvfstd34 = fit(cal_BmVobs,cal_BmV,'poly1');
    bvC3 = bvfstd34.p1;
    bvC4 = bvfstd34.p2;
    
    vifstd12 = fit(cal_VmI,iM_cat-cal_Iobs,'poly1');
    viC1 = vifstd12.p1;
    viC2 = vifstd12.p2;
    vifstd34 = fit(cal_VmIobs,cal_VmI,'poly1');
    viC3 = vifstd34.p1;
    viC4 = vifstd34.p2;
    
    
    % Correct all detections by zero point
    sB.Cat(bIdx,sB.Col.MAG_PSF) = sB.Cat(bIdx,sB.Col.MAG_PSF) + bZP;
    sV.Cat(vIdx,sV.Col.MAG_PSF) = sV.Cat(vIdx,sV.Col.MAG_PSF) + vZP;
    sI.Cat(iIdx,sI.Col.MAG_PSF) = sI.Cat(iIdx,sI.Col.MAG_PSF) + iZP;
    
    % And apply standar calibration equations
    BmV = bvC3*(sB.Cat(bIdx,sB.Col.MAG_PSF) - sV.Cat(vIdx,sV.Col.MAG_PSF)) + bvC4;
    V = sV.Cat(vIdx,sV.Col.MAG_PSF) + bvC1*BmV + bvC2;
    
    VmI = viC3*(sV.Cat(vIdx,sV.Col.MAG_PSF) - sI.Cat(iIdx,sI.Col.MAG_PSF)) + viC4;
    I = sI.Cat(iIdx,sI.Col.MAG_PSF) + viC1*VmI + viC2;
    
    % Plot the result as HR diagram
    figure(3)
    subplot(121)
    plot(BmV,V,'k.','MarkerSize',3)
    grid on
    axis([-0.3 2 13 22.5])
    set(gca, 'YDir','reverse')
    subplot(122)
    plot(VmI,V,'k.','MarkerSize',3)
    grid on
    axis([-0.2 1.3 13 22.5])
    set(gca, 'YDir','reverse')
    
    
    
    PixelData = reshape(sV.Im,[1,numel(sV.Im)]);
    NormPixelData = 1./(1 + exp((-PixelData + (0.7*median(PixelData) + 1*std(PixelData)))./(10*std(PixelData))));
    imNorm = reshape(NormPixelData,size(sV.Im));
    figure(4)
    imshow(imNorm)
    hold on
    plot(sV.Cat(vIdx,sB.Col.XWIN_IMAGE),sV.Cat(vIdx,sB.Col.YWIN_IMAGE),'go')
    
end