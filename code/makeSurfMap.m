function fig = makeSurfMap(dataPath,surfPath, varargin)
% Render and save retinotopic maps on the cortical surface
%
% Syntax:
%  makeSurfMap(dataPath,surfPath)
%
% Description
%
%
% Inputs:
%   dataPath              - Char vector. Full path (including file name) to
%                           a file in MGZ format that contains a
%                           retinotopic mapping result (e.g., eccentricity,
%                           polar angle) in native surface space.
%   surfpath              - Char vector. Path to the directory that
%                           contains various surface maps for the native
%                           space of the subject.
%
% Optional key/value pairs:
%  'mapType'              - Char vector. The type of map to be displayed;
%                           this controls the color gradient. Valid options
%                           include:
%                             {'ecc','pol','rsquared','sigma','areas'}
%  'hemisphere'           - Char vector. The hemisphere to display. Valid
%                           options are: {'lh','rh'}
%  'whichSurface'         - Char vector. The cortical surface on which to
%                           display the map. Valid options are:
%                           {'inflated','pial','sphere','white'}
%  'maxEccentricity'      - Scalar. The maximum eccentricity value
%                           displayed (only used when mapType is ecc)
%  'rsquaredDataPath'     - Char vector. The full path to the map of R^2
%                           values. If set, then the rsquaredThresh value
%                           is used to only display map values that are at
%                           locations with R^2 values above the threshold.
%  'rsquaredThresh'       - Scalar. The R^2 threshold that determines if a
%                           the map value for a vertex is displayed.
%  'alphaVal'             - Scalar. The transparency of the map color.
%  'showCurvature'        - Logical. Whether to show the cortical surface
%                           curvature on the render.
%  'colorRes'             - Scalar. The number of divisions of the color
%                           scale.
%
% Examples:
%{
    dataPath = '/tmp/flywheel/v0/output/maps_nativeMGZ/R_original_eccentricity_map.mgz'
    surfPath = '/tmp/flywheel/v0/input/structZip/TOME_3021/T1w/TOME_3021/surf'
    saveSurfMap(dataPath,surfPath,'hemisphere','rh','maxEccentricity',10)
%}
%{
    dataPath = '/tmp/flywheel/v0/output/maps_nativeMGZ/R_original_angle_map.mgz'
    surfPath = '/tmp/flywheel/v0/input/structZip/TOME_3021/T1w/TOME_3021/surf'
    saveSurfMap(dataPath,surfPath,'hemisphere','rh','mapType','pol')
%}
%{
    dataPath = '/tmp/flywheel/v0/output/maps_nativeMGZ/R_original_R2_map.mgz'
    surfPath = '/tmp/flywheel/v0/input/structZip/TOME_3021/T1w/TOME_3021/surf'
    saveSurfMap(dataPath,surfPath,'hemisphere','rh','mapType','rsquared')
%}
%{
    dataPath = '/tmp/flywheel/v0/output/maps_nativeMGZ/R_original_rfsize_map.mgz'
    surfPath = '/tmp/flywheel/v0/input/structZip/TOME_3021/T1w/TOME_3021/surf'
    saveSurfMap(dataPath,surfPath,'hemisphere','rh','mapType','sigma')
%}



%% Input parser
p = inputParser; p.KeepUnmatched = true;

% Required
p.addRequired('dataPath', @ischar);
p.addRequired('surfPath', @ischar);

% Optional key-value pairs
p.addParameter('mapType', 'ecc', @ischar);
p.addParameter('hemisphere','lh',@ischar);
p.addParameter('whichSurface','inflated',@ischar); % pial, white, or sphere
p.addParameter('maxEccentricity',90,@isscalar);
p.addParameter('maxGain',200,@isscalar);
p.addParameter('rsquaredDataPath','',@ischar);
p.addParameter('rsquaredThresh',0.1,@isscalar);
p.addParameter('alphaVal',0.85,@isscalar);
p.addParameter('showCurvature',true,@islogical);
p.addParameter('colorRes',200,@isscalar);
p.addParameter('visible',true,@islogical);

% parse
p.parse(dataPath, surfPath, varargin{:})


%% Load surface files
surfName = fullfile(surfPath,[p.Results.hemisphere '.' p.Results.whichSurface]);
[vert,face] = freesurfer_read_surf(surfName);


%% Prepare the curvature information
% Load the curvature map
curvName = fullfile(surfPath,[p.Results.hemisphere '.curv']);
[curv,~] = freesurfer_read_curv(curvName);

% Prepare the curvature map and cmap for plotting
if p.Results.showCurvature
    % Reverse the sign for plotting. Not sure why this is needed yet
    curv = -curv;
    
    % Binarize the curvature information
    ind.sulci = curv<0;
    ind.gyri = curv>0;
    ind.medial = curv==0;
    curv(ind.sulci) = .8;
    curv(ind.gyri) = 0.9;
    curv(ind.medial) = 0.7;
    cmap_curv = repmat(curv,1,3);
else
    curv = 0.85 * ones(size(curv));
    cmap_curv = repmat(curv,1,3);
end


%% Place the curvature and map information into a patch structure
brain.vertices = vert;
brain.faces = face;
brain.facevertexcdata = cmap_curv;


%% Load map data file
srf = load_mgh(dataPath);


%% Threshold data
if ~isempty(p.Results.rsquaredDataPath)
    rsquared_srf = load_mgh(p.Results.rsquaredDataPath);
    srf(rsquared_srf<p.Results.rsquaredThresh) = nan;
end


%% Define color map
switch p.Results.mapType
    case 'ecc'
        mapres=[1 p.Results.maxEccentricity p.Results.colorRes];
        mycolormap = make_ecc_colormap(mapres);
        typeLabel = 'Eccentricity [deg]';
        myvec = logspace(log10(mapres(1)),log10(mapres(2)),size(mycolormap,1));
    case 'pol'
        mapres=[-180 180 p.Results.colorRes];
        mycolormap = make_polar_colormap(mapres);
        typeLabel = 'Polar angle [deg]';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    case 'sigma'
        mapres=[0 10 p.Results.colorRes];
        mycolormap = flipud(jet(p.Results.colorRes));
        typeLabel = 'Sigma [deg]';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    case 'rsquared'
        mapres=[0 1 p.Results.colorRes];
        mycolormap = make_grayToRed_colormap(mapres);
        typeLabel = 'R^2';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    case 'hrfshift'
        mapres=[-4 4 p.Results.colorRes];
        mycolormap = make_blueToRed_colormap(mapres);
        typeLabel = 'shift hrf peak time [secs]';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    case 'gain'
        mapres=[0 p.Results.maxGain p.Results.colorRes];
        mycolormap = flipud(jet(p.Results.colorRes));
        typeLabel = 'response gain [T2* units]';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    case 'exponent'
        mapres=[0 2 p.Results.colorRes];
        mycolormap = flipud(jet(p.Results.colorRes));
        typeLabel = 'compressive exponent [au]';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    case 'varea'
        nAreas = length(unique(srf))-1;
        mapres = [1 nAreas nAreas];
        mycolormap = getDistinguishableColors(nAreas,{'w','k'});
        typeLabel = 'Visual area index';
        myvec = linspace(mapres(1),mapres(2),size(mycolormap,1));
    otherwise
        error('Unrecognized mapType');
end


%% Color vertices
cmap_vals = zeros(size(cmap_curv))+0.5;
alpha_vals = zeros(size(cmap_curv,1),1);
for ii = 1:length(srf)
    % Find the closest color value to the srf(ii) value
    [~,ind] = min(abs(myvec-srf(ii)));
    if isnan(srf(ii)) || srf(ii)==0
        col4thisvox = [.8 .8 .8]; % set nan and zero to gray
    else
        col4thisvox = mycolormap(ind,:);
    end
    cmap_vals(ii,:) = col4thisvox;
end

%% Set the alpha transparency
alpha_vals(~isnan(srf)) = p.Results.alphaVal;


%% Make figure
if p.Results.visible
    fig = figure('units','pixels','position',[0 0 1050 469]);
else
    fig = figure('units','pixels','position',[0 0 1050 469],'visible','off');
    
end

%% Plot brain and surface map
smp = brain;
smp.facevertexcdata = cmap_vals;

% Get the plot limits
rl_pa_is_Max = ceil(nanmax(brain.vertices)/10)*10;
rl_pa_is_Min = floor(nanmin(brain.vertices)/10)*10;

% Subplot range
spr = [1 7];

% Set the properties of each of the four view subplots
subSet = {[1 2 3 4 5]};
axisOrder = {[1 2]};

switch p.Results.hemisphere
    case 'lh'
        viewSet = {[45 0]};
    case 'rh'
        viewSet = {[-45 0]};
end

for vv = 1:length(subSet)
    % Set the frame for the plot
    subplot(spr(1),spr(2),subSet{vv});
    
    % Brain surface
    patch(brain,'EdgeColor','none','facecolor','interp','FaceAlpha',1);
    hold on
    
    % Color map
    patch(smp,'EdgeColor','none','facecolor','interp','FaceAlpha','flat',...
        'FaceVertexAlphaData',alpha_vals,'AlphaDataMapping','none');
    
    % Set the view
    view(viewSet{vv});
    
    % Hide the axis and make it tight
    axis tight equal off;
    
    % Set the aspect ratio and limits
    daspect([1 1 1]);
    ao = axisOrder{vv};
    xlim([rl_pa_is_Min(ao(1)) rl_pa_is_Max(ao(1))]);
    ylim([rl_pa_is_Min(ao(2)) rl_pa_is_Max(ao(2))]);
    
    % Camera settings
    camproj perspective; % orthographic; perspective
    lighting phong; % flat; gouraud; phong
    material dull; % shiny; metal; dull
    camlight('headlight');
    
end


%% Add the legend
subplot(spr(1),spr(2),[6 7]);

switch p.Results.mapType
    case 'ecc'
        x = linspace(-1,+1,mapres(3)); [xx,yy] = meshgrid(x);
        [~,r] = cart2pol(xx,yy);
        img = r;
        img(r>1) = nan;
        img = img/nanmax(img(:))*mapres(2);
        img(isnan(img)) = mapres(2);
        imagesc(x,x,img)
        colormap([mycolormap; 1 1 1])
        h = colorbar('southoutside');
        set(gca,'ColorScale','log') 
        h.Ticks = logspace(log10(1), log10(mapres(3)), 8);
        h.TickLabels = num2cell( round(logspace(log10(1), log10(p.Results.maxEccentricity), 8)*2)/2 );
    case 'pol'
        x = linspace(-1,1,mapres(3));
        [xx,yy] = meshgrid(x);
        [img,r] = cart2pol(xx,yy);
        img(r>1) = pi; % Set values outside of circle to black
        % Scale to the -180 to +180 regime
        img = img/nanmax(img(:))*mapres(2);
        % Rotate by 90 deg to get 0 as the upper vertical meridian
        img = rot90(img);
        imagesc(x,x,img)
        caxis([mapres(1) mapres(2)]);
        colormap([mycolormap; 1 1 1])
        h = colorbar('southoutside');
    case 'vareas'
        caxis([mapres(1) mapres(2)]);
        colormap(mycolormap);
        h = colorbar('southoutside');
        h.Ticks = linspace(1, nAreas, nAreas);
        ticLabels = {'V1','V2','V3','hV4','VO1','VO2','LO1','LO2','TO1','TO2','V3b','V3a'};
        h.TickLabels = ticLabels;
    otherwise
        caxis([mapres(1) mapres(2)]);
        colormap(mycolormap);
        h = colorbar('southoutside');
end
xlabel(h, typeLabel)
axis tight equal off


end % Main function


