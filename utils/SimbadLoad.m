function [CatData,idData] = SimbadLoad(s)
    [asRadius,centerRA,centerDec] = ImageRadius(s);
    formatRA = strsplit(StringRightAscension(centerRA*180/pi),':');
    formatDec = strsplit(StringDeclination(centerDec*180/pi),':');
    if centerDec >=0
        hexSign = '%2b';
    else
        hexSign = '%2d';
    end

    url = ['http://simbad.u-strasbg.fr/simbad/sim-coo?Coord=',...
    formatRA{1},'%20',formatRA{2},'%20',formatRA{3},'%20',hexSign,formatDec{1},'%20',formatDec{2},'%20',formatDec{3},...
    '&Radius=',num2str(asRadius),...
    '&Radius.unit=arcsec&submit=submit%20query&output.format=ASCII'];

    datapath = websave('SimbadData.txt',url);
    filestr = fileread(datapath);
    filebyline = regexp(filestr, '\n', 'split');
    lineDelimiter = find(cellfun(@isempty,filebyline));
    [~,gapIndex] = max(diff(lineDelimiter));
    extractedData = filebyline(lineDelimiter(gapIndex)+1:lineDelimiter(gapIndex+1)-2);
    filebyfield = regexp(extractedData, '\|', 'split');
    numfields = cellfun(@length, filebyfield);
    maxfields = max(numfields);
    fieldpattern = repmat({[]}, 1, maxfields);
    firstN = @(S,N) S(1:N);
    filebyfield = cellfun(@(S) firstN([S,fieldpattern], maxfields), filebyfield, 'Uniform', 0);
    fieldarray = vertcat(filebyfield{:});

    idCol = contains(fieldarray(1,:),'identifier');
    cooCol = contains(fieldarray(1,:),'coord');
    UmagCol = contains(fieldarray(1,:),'Mag U');
    BmagCol = contains(fieldarray(1,:),'Mag B');
    VmagCol = contains(fieldarray(1,:),'Mag V');
    RmagCol = contains(fieldarray(1,:),'Mag R');
    ImagCol = contains(fieldarray(1,:),'Mag I');

    fieldarray([1,2],:)=[];
    idData = fieldarray(:,idCol);

    Umag = string(cell2mat(fieldarray(:,UmagCol)));
    Umag(contains(Umag,'~'))="999";
    Umag = erase(Umag," ");

    Bmag = string(cell2mat(fieldarray(:,BmagCol)));
    Bmag(contains(Bmag,'~'))="999";
    Bmag = erase(Bmag," ");

    Vmag = string(cell2mat(fieldarray(:,VmagCol)));
    Vmag(contains(Vmag,'~'))="999";
    Vmag = erase(Vmag," ");

    Rmag = string(cell2mat(fieldarray(:,RmagCol)));
    Rmag(contains(Rmag,'~'))="999";
    Rmag = erase(Rmag," ");

    Imag = string(cell2mat(fieldarray(:,ImagCol)));
    Imag(contains(Imag,'~'))="999";
    Imag = erase(Imag," ");

    Coo = cell2mat(fieldarray(:,cooCol));
    RAref = [];
    Decref = [];
    for nSource = 1:size(Coo,1)
        tempCell = strsplit(Coo(nSource,:),' ');
        RAref = [RAref ; (str2double(tempCell{1}) + str2double(tempCell{2})/60 + str2double(tempCell{3})/3600)*180/12];
        Decref = [Decref ; (abs(str2double(tempCell{4})) + str2double(tempCell{5})/60 + str2double(tempCell{6})/3600)];
        if str2double(tempCell{4}) < 0
            Decref(end) = -Decref(end);
        end
    end


    CatData = [RAref Decref str2double(Umag) str2double(Bmag) str2double(Vmag) str2double(Rmag) str2double(Imag)];
end