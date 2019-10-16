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
p.addParameter('hemisphere','rh',@ischar);
p.addParameter('whichSurface','inflated',@ischar); % pial, white, or sphere
p.addParameter('stepSize',1,@isscalar);
p.addParameter('angleSet',10:40:170,@isnumeric);
p.addParameter('areaSet',[1],@isnumeric);



% Parse
p.parse(inferredMapsDirPath, surfPath, varargin{:})

% Pull some variables out of the param parser results
initialStepSize = p.Results.stepSize;
areaSet = p.Results.areaSet;
angleSet = p.Results.angleSet;

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
        pathPositions(positionIdx,angleIdx,areaIdx) = 0;
        pathAngles(positionIdx,angleIdx,areaIdx) = nan;
        
        % Re-initialize the distanceMap with the start map
        distanceMap = distanceMapStart;
        
        % Set the step size back to the initial value
        stepSize = initialStepSize;
        targetAngle = angleSet(angleIdx);
        
        %% Perform the path search
        % This is about four minutes on a single laptop core for
        %  stepSize = 1;
        stillWorkingFlag = true;
        while stillWorkingFlag
            
            % Iterate the position counter
            positionIdx = positionIdx+1;
            
            % Calculate an error metric
            errorMap = (areaAngleMap-targetAngle).^2 + 10.*(distanceMap-stepSize).^2;
            
            % Minimize the error and get the point
            [fVal, idx] =  min(errorMap);

            % Handle the output cases
            if fVal > 1000
                % We have an infinite error (or effectively infinite)
                % which means we have run out of visual area. Time to stop.
                stillWorkingFlag = false;
            elseif fVal > 50
                % We can't find a good next point. Try again with a bigger
                % step size.
                positionIdx = positionIdx-1;
                stepSize = stepSize + initialStepSize;                
            else
                % Store the values
                path(positionIdx,angleIdx,areaIdx) = idx;
                pathPositions(positionIdx,angleIdx,areaIdx) = ...
                    pathPositions(positionIdx-1,angleIdx,areaIdx) + distanceMap(idx);
                pathAngles(positionIdx,angleIdx,areaIdx) = angleMap(idx);
                
                % Make a new distance map from this point
                distanceMap = getDistanceMap(idx);
                
                % Constrain the distance map to this visual area
                distanceMap(vareaMap~=areaIdx)=Inf;
                
                % Salt the earth that we have already crossed
                distanceMap(distanceMapStart<=((positionIdx-1)*stepSize))=Inf;                

                % Re-set the stepSize in case it was changed
                stepSize = initialStepSize;
            end
            
        end % Loop over positions
    end % Loop over angle
end % Loop over visual areas


%% Obtain slope of magnification function per visual area
for areaIdx = 1:length(areaSet)
    p1 = nan(length(areaSet));
    p2 = nan(length(areaSet));
    for angleIdx = 1:length(angleSet)
        y = eccenMap(path(:,angleIdx,areaIdx));
        x = pathPositions(:,angleIdx,areaIdx);
        degPerMm = diff(y)./diff(x);
        lastIdx = find(y>10);
        lastIdx = lastIdx(1);
        myFit = fit(y(1:lastIdx),degPerMm(1:lastIdx),'poly1');
        plot(y(1:lastIdx),degPerMm(1:lastIdx),'x');
        hold on
        plot(y(1:lastIdx),myFit(y(1:lastIdx)),'-');
        p1(angleIdx)=myFit.p1;
        p2(angleIdx)=myFit.p2;
    end    
end


end % Main function

%% LOCAL FUNCTIONS

function plotSurfaceMap(vert, face, srf)

mycolormap = flipud(jet(200));
validIdx = logical(double(~isnan(srf)) .* double(~isinf(srf)));
myvec = linspace(nanmin(srf(validIdx)),nanmax(srf(validIdx)),200);

cmap_vals = repmat(zeros(size(srf))+0.5,1,3);

for ii = 1:length(srf)
    % Find the closest color value to the srf(ii) value
    [~,ind] = min(abs(myvec-srf(ii)));
    if isinf(srf(ii)) || srf(ii)==0
        col4thisvox = [.5 .5 .5]; % set nan and zero to gray
    else
        col4thisvox = mycolormap(ind,:);
    end
    cmap_vals(ii,:) = col4thisvox;
end

brain.vertices = vert;
brain.faces = face;
brain.facevertexcdata = cmap_vals;

% Brain surface
patch(brain,'EdgeColor','none','facecolor','interp','FaceAlpha',1);
daspect([1 1 1]);

% Camera settings
camproj perspective; % orthographic; perspective
lighting phong; % flat; gouraud; phong
material dull; % shiny; metal; dull
camlight('headlight');


end
