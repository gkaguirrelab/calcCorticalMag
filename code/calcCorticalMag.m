function corMag = calcCorticalMag(inferredMapsDirPath, surfPath, varargin)
%
%
%
% Examples:
%{
    inferredMapsDirPath = '/tmp/flywheel/v0/input/inferredSurfaces/opt/firstOutput';
    surfPath = '/tmp/flywheel/v0/input/structZip/TOME_3045/T1w/TOME_3045/surf/';

    calcCorticalMag(inferredMapsDirPath, surfPath)
%}

%% Parse inputs
p = inputParser; p.KeepUnmatched = false;

% Required
p.addRequired('inferredMapsDirPath',@isstr);
p.addRequired('surfPath',@isstr);

% Optional key-value pairs
p.addParameter('hemisphere','lh',@ischar);
p.addParameter('whichSurface','inflated',@ischar); % pial, white, or sphere
p.addParameter('stepSize',1,@isscalar);
p.addParameter('angleSet',10:40:170,@isnumeric);
p.addParameter('areaSet',[1 2 3],@isnumeric);



% Parse
p.parse(inferredMapsDirPath, surfPath, varargin{:})

% Pull some variables out of the param parser results
stepSize = p.Results.stepSize;
areaSet = p.Results.areaSet;
angleSet = p.Results.angleSet;


%% Load surface files
surfName = fullfile(surfPath,[p.Results.hemisphere '.' p.Results.whichSurface]);
[vert,face] = freesurfer_read_surf(surfName);

% An anonymous function to return the distance map for a given vertex
getDistanceMap = @(idx) heat_geodesic(vert,face,idx);


%% Load map data file
mapPath = fullfile(inferredMapsDirPath,[p.Results.hemisphere '.inferred_angle.mgz']);
angleMap = squeeze(load_mgh(mapPath));

mapPath = fullfile(inferredMapsDirPath,[p.Results.hemisphere '.inferred_eccen.mgz']);
eccenMap = squeeze(load_mgh(mapPath));

mapPath = fullfile(inferredMapsDirPath,[p.Results.hemisphere '.inferred_varea.mgz']);
vareaMap = squeeze(load_mgh(mapPath));


%% Loop over visual areas
for areaIdx = 1:length(areaSet)
    
    % Prepare the angle map for this visual area
    areaAngleMap = angleMap;
    areaAngleMap(vareaMap~=1)=Inf;
    
    % Find the start point for this visual area
    areaVerticies = find(vareaMap == areaIdx);
    [~,minPoint] = min(eccenMap(areaVerticies));
    startPoint = areaVerticies(minPoint);
    
    % Get the distance map from this start point. We save this and use it
    % to force the search to move forward.
    distanceMapStart = getDistanceMap(startPoint);
    distanceMapStart(vareaMap~=areaIdx)=Inf;
    
    
    %% Loop over iso-polar angle paths
    for angleIdx = 1:length(angleSet)
        
        % Copy the start point into the initial position for this path
        positionIdx = 1;
        path(positionIdx,angleIdx,areaIdx) = startPoint;
        pathErrors(positionIdx,angleIdx,areaIdx) = 0;
        pathPositions(positionIdx,angleIdx,areaIdx) = 0;
        
        % Re-initialize the distanceMap with the start map
        distanceMap = distanceMapStart;
        
        %% Perform the path search
        % This is about four minutes on a single laptop core for
        %  stepSize = 1;
        stillWorkingFlag = true;
        while stillWorkingFlag
            
            % Iterate the position counter
            positionIdx = positionIdx+1;
            
            % Anonymous function that provides an error metric
            myObj = @(x) (areaAngleMap-x).^2 + 100.*(distanceMap-stepSize).^2;
            
            % Minimize the error and get the point
            [fVal, idx] = ...
                min(myObj(angleSet(angleIdx)));
            
            % If we have an infinite error value, we have run out of visual
            % area
            if isinf(fVal)
                stillWorkingFlag = false;
            else
                
                % Store the values
                path(positionIdx,angleIdx,areaIdx) = idx;
                pathErrors(positionIdx,angleIdx,areaIdx) = fVal;
                pathPositions(positionIdx,angleIdx,areaIdx) = ...
                    pathPositions(positionIdx-1,angleIdx,areaIdx) + distanceMap(idx);
                
                % Make a new distance map from this point
                distanceMap = getDistanceMap(idx);
                
                % Constrain the distance map to this visual area
                distanceMap(vareaMap~=areaIdx)=Inf;
                
                % Salt the earth that we have already crossed
                distanceMap(distanceMapStart<((positionIdx-1)*stepSize))=Inf;
                
            end
            
        end % Loop over positions
    end % Loop over angle
end % Loop over visual areas

% Fit an exponential function to the eccentricity values
myExp = fittype( @(a,b,c, x) c.* (a.^(x./max(x))-b) );
y = eccenMap(path(:,angleIdx,areaIdx));
x = pathPositions(:,angleIdx,areaIdx);
myFit = fit(x,y,myExp);
yFit = myFit(x);

% This is inverse M
plot(yFit(1:end-1),diff(yFit)./diff(x))

end % Main function
