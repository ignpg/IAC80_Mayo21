function ImCut(fitspath,x1x2,y1y2)
    s = FITS.read2sim(fitspath);

    for idx = 1:numel(s)
        s(idx).Im = double(s(idx).Im(x1x2(1):x1x2(2),y1y2(1):y1y2(2)));
        s(idx).Header = FITS.cellhead_update(s(idx).Header,'NAXIS1',x1x2(2)-x1x2(1)+1,'length of data axis 1');
        s(idx).Header = FITS.cellhead_update(s(idx).Header,'NAXIS2',y1y2(2)-y1y2(1)+1,'length of data axis 2');

        [filepath,filename,ext] = fileparts(s(idx).ImageFileName);
        sim2fits(s(idx),'OutName',[filepath,'\c',filename,ext]);
    end
end