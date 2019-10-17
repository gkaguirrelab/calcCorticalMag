function calcCorticalMag(inferredMapsDirPath, surfPath, outPath, varargin)
%
%
%
% Examples:
%{
    inferredMapsDirPath = '/tmp/flywheel/v0/input/inferredSurfaces/opt/firstOutput';
    surfPath = '/tmp/flywheel/v0/input/structZip/TOME_3045/T1w/TOME_3045/surf/';
    outPath = '/tmp/flywheel/v0/output/';
    calcCorticalMag(inferredMapsDirPath, surfPath, outPath)
%}

%% Parse inputs
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('inferredMapsDirPath',@isstr);
p.addRequired('surfPath',@isstr);
p.addRequired('outPath',@isstr);

% Optional key-value pairs
p.addParameter('hemisphere','rh',@ischar);
p.addParameter('whichSurface','inflated',@ischar); % pial, white, or sphere
p.addParameter('stepSize',1,@isscalar);
p.addParameter('stepTol',4,@isnumeric);
p.addParameter('angleSet',10:40:170,@isnumeric);
p.addParameter('angleTol',10,@isnumeric);
p.addParameter('areaSet',[1],@isnumeric);



% Parse
p.parse(inferredMapsDirPath, surfPath, outPath, varargin{:})

% Pull some variables out of the param parser results
initialStepSize = p.Results.stepSize;
stepTol = p.Results.stepTol;
areaSet = p.Results.areaSet;
angleSet = p.Results.angleSet;
angleTol = p.Results.angleTol;

% The polar angle values will be negative in the right hemisphere
if strcmp(p.Results.hemisphere,'rh')
%    angleSet = -angleSet;
end

%% Load surface files
surfName = fullfile(surfPath,[p.Results.hemisphere '.' p.Results.whichSurface]);
[vert,face] = freesurfer_read_surf(surfName);

% An anonymous function to return the distance map for a given vertex. This
% function uses the Geodesics in Heat algorithm:
%
%   Crane, Keenan, Clarisse Weischedel, and Max Wardetzky. "Geodesics in
%   heat: A new approach to computing distance based on heat flow." ACM
%   Transactions on Graphics (TOG) 32.5 (2013): 152.
%
% as implemented by Alec Jacobson in the gptoolbox:
%   https://github.com/alecjacobson/gptoolbox/
%
getDistanceMap = @(idx) heat_geodesic(vert,face,idx);


%% Load map data file
mapPath = fullfile(inferredMapsDirPath,[p.Results.hemisphere '.inferred_angle.mgz']);
angleMap = squeeze(load_mgh(mapPath));

mapPath = fullfile(inferredMapsDirPath,[p.Results.hemisphere '.inferred_eccen.mgz']);
eccenMap = squeeze(load_mgh(mapPath));

mapPath = fullfile(inferredMapsDirPath,[p.Results.hemisphere '.inferred_varea.mgz']);
[vareaMap, M, mr_parms, volsz] = load_mgh(mapPath);
vareaMap = squeeze(vareaMap);

%% Find all valid vertices
validIdx = find(vareaMap~=0);

cmfMap = nan(size(vareaMap));
for ii=1:length(validIdx)
    neighborIdx = unique(face(any((face==validIdx(ii))'),:));
    thisCortCoord = squeeze(vert(validIdx(ii),:));
    thisVisCoord = [ eccenMap(validIdx(ii)) .* cosd(angleMap(validIdx(ii))), ...
                     eccenMap(validIdx(ii)) .* sind(angleMap(validIdx(ii)))];

    distancesMm = sqrt(sum((vert(neighborIdx,:)-thisCortCoord).^2,2));
    distancesDeg = vecnorm(thisVisCoord - ...
    [ eccenMap(neighborIdx) .* cosd(angleMap(neighborIdx)), ...
                     eccenMap(neighborIdx) .* sind(angleMap(neighborIdx))],2,2);
    cmfMap(validIdx(ii)) = nanmean(nanmean(distancesDeg./distancesMm));
end

mapPathOut = fullfile(outPath,[p.Results.hemisphere '.cmf.mgz']);
save_mgh(reshape(cmfMap,volsz), mapPathOut, M, mr_parms);

end % Main function
