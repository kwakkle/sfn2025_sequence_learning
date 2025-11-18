function stdshade_smooth3(amatrix,F,alpha,acolor,smth,shadetype)
%20220913: Zhaoran Zhang: F is the real x value

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end

if exist('smth','var'); if isempty(smth); smth=1; end
else smth=1; %no smoothing by default
end  

if ne(size(F,1),1)
    F=F';
end

% amean = nanmean(amatrix,1); %get man over first dimension
amean = mean(amatrix,'omitnan');
if smth > 1
%     amean = smooth(nanmean(amatrix,1),smth)'; %use boxfilter to smooth data
    amean = smooth(mean(amatrix,1,'omitnan'),smth)';
end

% if nargin < 5
%     shadetype='SEM';
% end
if strcmp(shadetype,'STD')==1
% astd = nanstd(amatrix,[],1); % to get std shading
astd = std(amatrix,[],1,'omitnan'); % to get std shading
elseif strcmp(shadetype,'SEM')==1
% astd = smooth(nanstd(amatrix,[],1)/sqrt(size(amatrix,1)),smth)'; % to get sem shading
astd = smooth(std(amatrix,[],1,'omitnan')/sqrt(size(amatrix,1)),smth)'; % to get sem shading
end

if exist('alpha','var')==0 || isempty(alpha) 
    fill([(F) (fliplr(F))],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else fill([(F) (fliplr(F))],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');    
end

if ishold==0
    check=true; else check=false;
end

hold on;plot((F),amean,'color',acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end


function dataOut = boxFilter(dataIn, fWidth)
% apply 1-D boxcar filter for smoothing

dataStart = cumsum(dataIn(1:fWidth-2),2);
dataStart = dataStart(1:2:end) ./ (1:2:(fWidth-2));
dataEnd = cumsum(dataIn(length(dataIn):-1:length(dataIn)-fWidth+3),2);
dataEnd = dataEnd(end:-2:1) ./ (fWidth-2:-2:1);
dataOut = conv(dataIn,ones(fWidth,1)/fWidth,'full');
dataOut = [dataStart,dataOut(fWidth:end-fWidth+1),dataEnd];

end
