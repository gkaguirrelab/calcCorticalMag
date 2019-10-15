function calcCorticalMagLocalHook
% calcCorticalMagLocalHook
%
% For use with the ToolboxToolbox.  If you copy this into your
% ToolboxToolbox localToolboxHooks directory (by default,
% ~/localToolboxHooks) and delete "LocalHooksTemplate" from the filename,
% this will get run when you execute tbUse({'calcCorticalMagConfig'}) to set up for
% this project.  You then edit your local copy to match your local machine.
%
% The main thing that this does is define Matlab preferences that specify input and output
% directories.
%
% You will need to edit the project location and i/o directory locations
% to match what is true on your computer.

 
%% Define project
projectName = 'calcCorticalMag';
 

%% Clear out old preferences
if (ispref(projectName))
    rmpref(projectName);
end


%% Specify and save project location
projectBaseDir = tbLocateProject(projectName);
setpref(projectName,'projectBaseDir',projectBaseDir);


%% Set flywheel API key as a preference
flywheelAPIKey='Copy this value from flywheel and paste here';
setpref(projectName,'flywheelAPIKey',flywheelAPIKey);


%% Set the workbench command path
if ismac
setpref(projectName,'wbCommand','/Applications/workbench/bin_macosx64/wb_command');
end


%% Get the userID
[~, userID] = system('whoami');
userID = strtrim(userID);


%% Paths to store flywheel data and scratch space
if ismac
    % Code to run on Mac plaform
    setpref(projectName,'flywheelScratchDir','/tmp/flywheel');
    setpref(projectName,'flywheelRootDir',fullfile('/Users/',userID,'/Documents/flywheel'));
elseif isunix
    % Code to run on Linux plaform
    setpref(projectName,'flywheelScratchDir','/tmp/flywheel');
    setpref(projectName,'flywheelRootDir',fullfile('/home/',userID,'/Documents/flywheel'));
elseif ispc
    % Code to run on Windows platform
    warning('No supported for PC')
else
    disp('What are you using?')
end


%% Check for required Matlab toolboxes
requiredAddOns = {...
    'Optimization Toolbox',...                   % optimization_toolbox
    'Statistics and Machine Learning Toolbox'... % statistics_toolbox
    };
% Given this hard-coded list of add-on toolboxes, we then check for the
% presence of each and issue a warning if absent.
V = ver;
VName = {V.Name};
warnState = warning();
warning off backtrace
for ii=1:length(requiredAddOns)
    if ~any(strcmp(VName, requiredAddOns{ii}))
        warnString = ['The Matlab ' requiredAddOns{ii} ' is missing. ' projectName ' may not function properly.'];
        warning('localHook:requiredMatlabToolboxCheck',warnString);
    end
end
warning(warnState);


%% Ensure that python3 and neuopythy are installed

% Check for some flavor of Python v3
if floor(str2double(pyversion)) ~= 3
	warning('localHook:pythonVersion',['The routines expect Python v3, but v' pyversion ' is installed']);
end

% Make sure that numpy and neuropythy are present
fprintf('Installing neuropythy\n');
command = 'pip install neuropythy > /dev/null';
system(command);

end
