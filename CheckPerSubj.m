%% Check each subject
% 1.Hand trace across all trials (color shade)
% 2.Eye and hand trace of example trials
% 3.Reaction time across all trials
% 4.Hand and eye error across all trials
%%%%%
clear,clc
%%
Load_data = 1;
Hand_trace = 1;
Eye_hand_trace = 1;
RT = 1;
HandEyeErr = 1;

% define working location: either on laptop or in the lab
computer = 'laptop';
switch computer
    case 'laptop' 
        basePath = 'D:\code\lab\Jana eye-hand\';
    case 'lab'
        basePath = 'I:\My Drive\copy MATLAB\';
    case 'opher_lab'
        basePath = 'D:\Opher\OneDrive - Ben Gurion University of the Negev\Qsync\2021 Spring\Eye tracking';
    case 'opher_laptop'
        basePath = 'C:\Users\donch\Qsync\2020 Fall\Grants\Donchin ISF 2019\Matlab\Sivan';
end
functionPath = fullfile(basePath, 'functions');
exportPath = fullfile(basePath,'altmany-export_fig');
handDataPath = fullfile(basePath, 'eyetracking feedback length\handDataNolanmk_test');  % where hand files are 
eyeDataPath = fullfile(basePath, 'eyetracking feedback length\eyeDataNolanmk_test');    % where eye files are
outputFolder = 'D:\code\lab\data\checksubj';  % where figures are saved 
%% Load data
if Load_data
% add relevant paths
addpath(handDataPath);
addpath(eyeDataPath);
dir2Fullpath = @(d) fullfile(d.folder, d.name);

    eyeDirs = dir(eyeDataPath);
    eyeDirs = eyeDirs(~ismember({eyeDirs.name}, {'.', '..'})); % avoid using files . and ..
    eyeDirs = eyeDirs([eyeDirs.isdir]);
    eyeSubjs = string({eyeDirs.name}');
    
    handDir = dir(fullfile(handDataPath, '*.csv'));
    handFileNames = string({handDir.name}');   
    allSubjs = extractBefore(handFileNames, 4);
    disp(strjoin(['Checking data from the following subjects: ' allSubjs']));
    numSubj = length(allSubjs);
    
    % Define scales for data
    tMin = -1.5;
    tMax = 2.5;
    dt = 0.01; 
    tVals = (tMin:dt:tMax)';
    numTVals = length(tVals);    
    % xUNit should be 1 because all units in psychopy are in cm 
    xUnit = 1; 
    yUnit = 1;
    screenWidth = 53;
    screenHeight = 30; 
    % The eye starts from the bottom left of the screen and has normalized
    % coordinates over the screen. 
    x0Eye = 0.5; 
    y0Eye = 0.5;
    xUnitEye = screenWidth;
    yUnitEye = screenHeight;  
    % Determine a confidence threshold for the gaze data
    confidenceThreshold = 0.97;  
    % Determine a dispersion threshold for the fixation data 
    dispersionThresh = 1;  
    % Define indexes as structs
    % Use these structs as a way of using symbolic names to index into arrays
    % (Like a poor man's enums)
    block = struct('baseline', 1, 'rotation', 2, 'aftereffect', 3);
    numBlocks = length(fieldnames(block));
    phase = struct('ring', 1, 'hold', 2, 'move', 3, 'feedback', 4);
    numPhase = length(fieldnames(phase));
    info = struct('startTime', 1, 'mouseStart', 2, 'mouseX', 3, 'mouseY', 4, 'mouseTime', 5); %这几个不太明确理解
    numInfo = length(fieldnames(info));
    coord = struct('x', 1, 'y', 2);
    numCoord = 2;    
    % Define handColumns
    % Define the columns for the hand data .csv file that we will use
    handColumns = strings(numBlocks, numPhase, numInfo);    
    % **Baseline**    
    handColumns(block.baseline, phase.ring, info.startTime) = "Ring_started";
    handColumns(block.baseline, phase.ring, info.mouseStart) = "ring_mouse_started";
    handColumns(block.baseline, phase.ring, info.mouseX) = "ring_mouse_x";
    handColumns(block.baseline, phase.ring, info.mouseY) = "ring_mouse_y";
    handColumns(block.baseline, phase.ring, info.mouseTime) = "ring_mouse_time";
    
    handColumns(block.baseline, phase.hold, info.startTime) = "baseline_fixation_started";
    handColumns(block.baseline, phase.hold, info.mouseStart) = "baseline_mouse_started";
    handColumns(block.baseline, phase.hold, info.mouseX) = "baseline_mouse_x";
    handColumns(block.baseline, phase.hold, info.mouseY) = "baseline_mouse_y";
    handColumns(block.baseline, phase.hold, info.mouseTime) = "baseline_mouse_time";
    
    handColumns(block.baseline, phase.move, info.startTime) = "baseline_fixation2_started";
    handColumns(block.baseline, phase.move, info.mouseStart) = "baseline_move_mouse_started";
    handColumns(block.baseline, phase.move, info.mouseX) = "baseline_move_mouse_x";
    handColumns(block.baseline, phase.move, info.mouseY) = "baseline_move_mouse_y";
    handColumns(block.baseline, phase.move, info.mouseTime) = "baseline_move_mouse_time";
    
    handColumns(block.baseline, phase.feedback, info.startTime) = "baseline_fixation3_started";
    handColumns(block.baseline, phase.feedback, info.mouseStart) = "baseline_feedback_mouse_started";
    handColumns(block.baseline, phase.feedback, info.mouseX) = "baseline_feedback_mouse_x";
    handColumns(block.baseline, phase.feedback, info.mouseY) = "baseline_feedback_mouse_y";
    handColumns(block.baseline, phase.feedback, info.mouseTime) = "baseline_feedback_mouse_time";
    
    % **Rotated**   
    handColumns(block.rotation, phase.ring, info.startTime) = "Ring2_started";
    handColumns(block.rotation, phase.ring, info.mouseStart) = "ring_mouse2_started";
    handColumns(block.rotation, phase.ring, info.mouseX) = "ring_mouse2_x";
    handColumns(block.rotation, phase.ring, info.mouseY) = "ring_mouse2_y";
    handColumns(block.rotation, phase.ring, info.mouseTime) = "ring_mouse2_time";
    
    handColumns(block.rotation, phase.hold, info.startTime) = "rot_hold_fixation_started";
    handColumns(block.rotation, phase.hold, info.mouseStart) = "rot_hold_mouse_started";
    handColumns(block.rotation, phase.hold, info.mouseX) = "rot_hold_mouse_x";
    handColumns(block.rotation, phase.hold, info.mouseY) = "rot_hold_mouse_y";
    handColumns(block.rotation, phase.hold, info.mouseTime) = "rot_hold_mouse_time";
    
    handColumns(block.rotation, phase.move, info.startTime) = "rotated_fixation_2_started";
    handColumns(block.rotation, phase.move, info.mouseStart) = "rotated_move_mouse_started";
    handColumns(block.rotation, phase.move, info.mouseX) = "rotated_move_mouse_x";
    handColumns(block.rotation, phase.move, info.mouseY) = "rotated_move_mouse_y";
    handColumns(block.rotation, phase.move, info.mouseTime) = "rotated_move_mouse_time";
    
    handColumns(block.rotation, phase.feedback, info.startTime) = "rotated_fixation2_started";
    handColumns(block.rotation, phase.feedback, info.mouseStart) = "rotated_feedback_mouse_started";
    handColumns(block.rotation, phase.feedback, info.mouseX) = "rotated_feedback_mouse_x";
    handColumns(block.rotation, phase.feedback, info.mouseY) = "rotated_feedback_mouse_y";
    handColumns(block.rotation, phase.feedback, info.mouseTime) = "rotated_feedback_mouse_time";
    
    % **Aftereffect**
    
    handColumns(block.aftereffect, phase.ring, info.startTime) = "Ring_started";
    handColumns(block.aftereffect, phase.ring, info.mouseStart) = "ring_mouse_started";
    handColumns(block.aftereffect, phase.ring, info.mouseX) = "ring_mouse_x";
    handColumns(block.aftereffect, phase.ring, info.mouseY) = "ring_mouse_y";
    handColumns(block.aftereffect, phase.ring, info.mouseTime) = "ring_mouse_time";
    
    handColumns(block.aftereffect, phase.hold, info.startTime) = "baseline_fixation_started";
    handColumns(block.aftereffect, phase.hold, info.mouseStart) = "baseline_mouse_started";
    handColumns(block.aftereffect, phase.hold, info.mouseX) = "baseline_mouse_x";
    handColumns(block.aftereffect, phase.hold, info.mouseY) = "baseline_mouse_y";
    handColumns(block.aftereffect, phase.hold, info.mouseTime) = "baseline_mouse_time";
    
    handColumns(block.aftereffect, phase.move, info.startTime) = "baseline_fixation2_started";
    handColumns(block.aftereffect, phase.move, info.mouseStart) = "baseline_move_mouse_started";
    handColumns(block.aftereffect, phase.move, info.mouseX) = "baseline_move_mouse_x";
    handColumns(block.aftereffect, phase.move, info.mouseY) = "baseline_move_mouse_y";
    handColumns(block.aftereffect, phase.move, info.mouseTime) = "baseline_move_mouse_time";
    
    handColumns(block.aftereffect, phase.feedback, info.startTime) = "aftereffect_fixation_started";
    handColumns(block.aftereffect, phase.feedback, info.mouseStart) = "aftereffect_mouse_started";
    handColumns(block.aftereffect, phase.feedback, info.mouseX) = "aftereffect_mouse_x";
    handColumns(block.aftereffect, phase.feedback, info.mouseY) = "aftereffect_mouse_y";
    handColumns(block.aftereffect, phase.feedback, info.mouseTime) = "aftereffect_mouse_time";
     
    trialNumberCols(block.baseline) = "trials_thisN";   
    trialNumberCols(block.rotation) = "trials_3_thisN"; 
    trialNumberCols(block.aftereffect) = "trials_after_thisN"; 
    loopNumberCol = "trials_4_thisN";  
    numTrialsPerLoop = 20; 
    numLoops = 4;  
    numBaselineTrials = 60; 
    numRotationTrials = numLoops*numTrialsPerLoop;
     
    % Load hand data 
    emptyVec = zeros(0,1);
    emptyTable = table(emptyVec, emptyVec, emptyVec, 'VariableNames', {'t', 'x', 'y'});
    blockStruct = struct('baseline', {}, 'rotation', {}, 'aftereffect', {});  
    allData = struct();
    for subjNum = 1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;     
        handIndex = find(strncmpi(handFileNames, subjID,3));
        if ~isempty(handIndex)
            handFile = dir2Fullpath(handDir(handIndex));  
            opts = detectImportOptions(handFile);         
            opts.SelectedVariableNames = [handColumns(:); ...
                trialNumberCols(:); loopNumberCol; "target_x"; "target_y"; "group"; ...
                "MoveOnsetTime1"; "MoveEndTime1"; "MoveOnsetTime2"; "MoveEndTime2"];      
            opts = setvartype(opts, trialNumberCols, {'double' 'double' 'double'});  
            handDataFile = readtable(handFile, opts);  
            numRows = size(handDataFile,1);   
            
            % Make list of trials where there is data.
            trialNums = zeros(numRows, 1);
            trialNums(all(isnan(handDataFile{:, trialNumberCols}),2)) = NaN; 
            numTrials = sum(~isnan(trialNums));  
            trialNums(~isnan(trialNums)) = 1:numTrials; 
   
            allData.(sF).group = handDataFile.group(1); 
            allData.(sF).hand = NaN*zeros(numTVals, numTrials, numCoord); 
            allData.(sF).target = NaN*zeros(numTrials, numCoord);
            allData.(sF).block = strings(numTrials, 1);
            allData.(sF).time0 = NaN*zeros(numTrials, 1);  
            allData.(sF).timeFeedback = NaN*zeros(numTrials, 1);
            allData.(sF).timeStart = NaN*zeros(numTrials, 1);
            allData.(sF).timeEnd = NaN*zeros(numTrials, 1);
            allData.(sF).trialNum = NaN*zeros(numTrials, 1);
                
            blockName = "baseline"; % For empty trials, blockName goes to the one before, so we will start with baseline
            dataColumn = "baseline_move";
            for trialNum = 1:numTrials
                thisRow = handDataFile(trialNums == trialNum,:);                 
                if ~isempty(thisRow.aftereffect_mouse_x{1})  
                    blockName = "aftereffect";
                elseif ~isempty(thisRow.rotated_move_mouse_x{1})
                    blockName = "rotation";
                elseif ~isempty(thisRow.baseline_move_mouse_x{1})
                    blockName = "baseline";
                else
                    blockName = "";
                    warning("No data in trial "+string(trialNum)+" for subject "+subjID);
                end                
                if blockName ~= ""
                    thisBlock = block.(blockName);                    
                    
                    thisData = allData.(sF).hand;  
                    
                    trialTableRaw = emptyTable;
                    
                    thisTrialTimeStart = thisRow.(handColumns(thisBlock, phase.ring, info.startTime));  
                    thisTrialTime0 = thisRow.MoveOnsetTime1; 
                    thisTrialEnd =  thisRow.MoveEndTime1;  
                    thisTrialOnset2 = thisRow.MoveOnsetTime2; 
                    thisTrialEnd2 =  thisRow.MoveEndTime2;
                    thisTrialTargetOnset = thisRow.(handColumns(thisBlock, phase.move, info.startTime));  
                    thisTrialHandInOrigin = thisRow.(handColumns(thisBlock, phase.hold, info.startTime)); 
                    for phaseNum = 1:numPhase
                        thisPhaseStart = thisRow.(handColumns(thisBlock, phaseNum, info.startTime)); 
                        if iscell(thisPhaseStart); warning("Converting cell " + thisPhaseStart + " to number"); thisPhaseStart = str2double(thisPhaseStart{1}); end
                        thisPhaseMouseStart = thisRow.(handColumns(thisBlock, phaseNum, info.mouseStart));                         
                        thisTimeString = thisRow.(handColumns(thisBlock, phaseNum, info.mouseTime)){1};                         
                        thisXString = thisRow.(handColumns(thisBlock, phaseNum, info.mouseX)){1};
                        thisYString = thisRow.(handColumns(thisBlock, phaseNum, info.mouseY)){1};
                        thisTable = emptyTable;
                        if ~isempty(thisXString)
                            x = eval(thisXString)'*xUnit;
                            y = eval(thisYString)'*yUnit;
                            % t{phaseNum} = eval(thisTimeString)' + thisPhaseMouseStart + thisPhaseStart - thisTrialTime0;
                            t1{phaseNum} = eval(thisTimeString)'; %#ok<SAGROW>
                            t = t1{phaseNum};
                            
                            posTimes = t > 0;  
                            x = x(posTimes);
                            y = y(posTimes);
                            t = t(posTimes);
                            % t1{phaseNum} = t1{phaseNum} - t1{phaseNum}(1) + thisPhaseStart - thisTrialTime0; %#ok<SAGROW>
                            t = t + (thisPhaseStart - thisTrialTime0);  
                                                                       
                            thisTable = table(t, x, y);  
                        end
                        trialTableRaw = [trialTableRaw; thisTable]; %#ok<AGROW> 
                    end 
                    x = interp1(trialTableRaw.t, trialTableRaw.x, tVals, 'makima', NaN); 
                    y = interp1(trialTableRaw.t, trialTableRaw.y, tVals, 'makima', NaN);
                    
                    thisTrialTimeEnd = thisTrialTimeStart + max(trialTableRaw.t);                
                    if isempty(thisData)
                        thisData = NaN*zeros(numTVals, numBlockTrials, 2);
                    end
                    thisData(:,trialNum,coord.x) = x;
                    thisData(:,trialNum,coord.y) = y;
                    
                    allData.(sF).hand = thisData;
                    
                    thisTargetX = thisRow.target_x*xUnit;  
                    thisTargetY = thisRow.target_y*yUnit;
                    thisTimeFeedback  = thisRow.(handColumns(thisBlock, phase.feedback, info.mouseStart));
                    if ~isnumeric(thisTimeFeedback)
                        thisTimeFeedback = NaN;
                    end
                    
                    thisTrialNumInBlock = thisRow.(trialNumberCols(thisBlock));
                    switch blockName
                        case "baseline" 
                            thisTrialNum = thisTrialNumInBlock + 1;
                            thisTrialStartLabel = string(thisTrialNumInBlock);
                            thisTrialEndLabel = thisTrialStartLabel;
                        case "rotation"
                            thisLoopNumber = thisRow.(loopNumberCol)(1);
                            thisTrialNum = thisTrialNumInBlock + thisLoopNumber*numTrialsPerLoop + numBaselineTrials + 1;
                            thisTrialStartLabel = thisTrialNumInBlock + " loop " + thisLoopNumber;
                            thisTrialEndLabel = thisTrialNumInBlock + "loop " + thisLoopNumber; % note the space
                        case "aftereffect"
                            thisTrialNum = thisTrialNumInBlock + numBaselineTrials + numRotationTrials;
                            thisTrialStartLabel = "60";
                            thisTrialEndLabel = string(thisTrialNumInBlock);
                        otherwise
                            error("Unexpected block name in " + sF + " trial " + trialNum);
                    end
                    
                    allData.(sF).target(trialNum, :) = [thisTargetX thisTargetY]; 
                    allData.(sF).block(trialNum) = blockName;
                    allData.(sF).timeFeedback(trialNum) = thisTimeFeedback - thisTrialTime0;
                    allData.(sF).time0(trialNum) = thisTrialTime0;  
                    allData.(sF).endTime(trialNum) = thisTrialEnd;
                    allData.(sF).timeOnset2(trialNum) = thisTrialOnset2;  
                    allData.(sF).endTime2(trialNum) = thisTrialEnd2;
                    allData.(sF).targetOnset(trialNum) = thisTrialTargetOnset;  
                    allData.(sF).handInOrigin(trialNum) = thisTrialHandInOrigin;  
                    allData.(sF).timeStart(trialNum) = thisTrialTimeStart; 
                    allData.(sF).timeEnd(trialNum) = thisTrialTimeEnd;
                    allData.(sF).trialNum(trialNum) = thisTrialNum; 
                    allData.(sF).trialStartLabel(trialNum) = thisTrialStartLabel; 
                    allData.(sF).trialEndLabel(trialNum) = thisTrialEndLabel;
                    end
                end 
            end
            if size(allData.(sF).target,1) > 160  
                allData.(sF).hand = allData.(sF).hand(:,[1:140 221:240],:);
                allData.(sF).target = allData.(sF).target([1:140 221:240],:);
                allData.(sF).block = allData.(sF).block([1:140 221:240]);
                allData.(sF).timeFeedback = allData.(sF).timeFeedback([1:140 221:240]);
                allData.(sF).time0 = allData.(sF).time0([1:140 221:240]);
                allData.(sF).endTime = allData.(sF).endTime([1:140 221:240]);
                allData.(sF).timeOnset2 = allData.(sF).timeOnset2([1:140 221:240]);
                allData.(sF).endTime2 = allData.(sF).endTime2([1:140 221:240]);
                allData.(sF).targetOnset = allData.(sF).targetOnset([1:140 221:240]);
                allData.(sF).handInOrigin = allData.(sF).handInOrigin([1:140 221:240]);
                allData.(sF).timeStart = allData.(sF).timeStart([1:140 221:240]);
                allData.(sF).timeEnd = allData.(sF).timeEnd([1:140 221:240]);
                allData.(sF).trialNum = allData.(sF).trialNum([1:140 221:240]);
                allData.(sF).trialStartLabel = allData.(sF).trialStartLabel([1:140 221:240]);
                allData.(sF).trialEndLabel = allData.(sF).trialEndLabel([1:140 221:240]);
                end
        end    
    % Load eye data  
    for subjNum = 1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;      
        eyeIndex = find(contains(eyeSubjs, subjID));
        if ~isempty(eyeIndex)
            eyeDir = dir2Fullpath(eyeDirs(eyeIndex));            
            annotationsFile = fullfile(eyeDir, 'annotations.csv');  
            optsAnnotations = detectImportOptions(annotationsFile);
            optsAnnotations.SelectedVariableNames = {'index', 'timestamp','label'};  
            optsAnnotations = setvartype(optsAnnotations,'label',{'string'});
            annotationDataFile = readtable(annotationsFile, optsAnnotations);            
            fixFile = fullfile(eyeDir, 'fixations_on_surface_screen.csv');  
            if exist(fixFile,'file')
                optsFix = detectImportOptions(fixFile);
                if length(optsFix.SelectedVariableNames) < 2
                    fixDataFile = [];
                else
                    optsFix.SelectedVariableNames = {'start_timestamp','duration', ... 
                        'norm_pos_x', 'norm_pos_y', 'dispersion'};
                    fixDataFile = readtable(fixFile, optsFix);
                    allData.(sF).fix = NaN*zeros(numTVals, numTrials, numCoord);
                end
            else
                fixDataFile = [];
            end
            
            gazeFile = fullfile(eyeDir, 'gaze_positions_on_surface_screen.csv');  
            optsGaze = detectImportOptions(gazeFile);
            if length(optsGaze.SelectedVariableNames) < 2
                gazeDataFile = [];
            else
                optsGaze.SelectedVariableNames = {'gaze_timestamp','confidence', ...  
                    'x_norm', 'y_norm'};
                gazeDataFile = readtable(gazeFile);
                allData.(sF).timeGazeStart = NaN*zeros(numTrials, 1);
                allData.(sF).timeGazeEnd = NaN*zeros(numTrials, 1);
                allData.(sF).gaze = NaN*zeros(numTVals, numTrials, numCoord);
            end     
            trialLabel = annotationDataFile.label;  
            timestamp = annotationDataFile.timestamp;  
            trialStartLabels = find(startsWith(trialLabel, "start_trial"));  
            trialEndLabels = find(startsWith(trialLabel, "end_trial"));
            movementStartLabels = find(startsWith(trialLabel, "move_start"));
            movementEndLabels = find(startsWith(trialLabel, "move_end"));          
            numTrialStartLabels = length(trialStartLabels);
            if numTrialStartLabels ~= length(trialEndLabels)
                error("Not the same number of start and end trial label in subject " + sF);
            end            
            if numTrialStartLabels > 160
                trialStartLabels = trialStartLabels([1:140 221:240]);
                trialEndLabels = trialEndLabels([1:140 221:240]);
                movementStartLabels = movementStartLabels([1:140 221:240]);
                movementEndLabels = movementEndLabels([1:140 221:240]);
                numTrialStartLabels = length(movementEndLabels);
            end            
            expStart = timestamp(startsWith(trialLabel, "start_experiment")); 
            trialStarts = zeros(numTrialStartLabels,1);  
            trialEnds = zeros(numTrialStartLabels,1);    
            movementStarts = zeros(numTrialStartLabels,1); 
            movementEnds = zeros(numTrialStartLabels,1);   
            for trialNum = 1:numTrialStartLabels
                thisTrialStartLabel = allData.(sF).trialStartLabel(trialNum);
                thisTrialEndLabel = allData.(sF).trialEndLabel(trialNum);                
                if trialLabel(trialStartLabels(trialNum)) ~= "start_trial " + thisTrialStartLabel ...
                        || trialLabel(trialEndLabels(trialNum)) ~= "end_trial "+thisTrialEndLabel
                    error("The labels in the gaze file of " + trialLabel(trialStartLabels(trialNum)) + ...
                        " and " + trialLabel(trialEndLabels(trialNum)) + " in subject " +sF + ...
                        " are not like they should be");
                end               
                trialStarts(trialNum) = timestamp(trialStartLabels(trialNum)); 
                trialEnds(trialNum) = timestamp(trialEndLabels(trialNum)); 
                
                thisMovementStartIndex = trialStartLabels(trialNum) < movementStartLabels ...  
                    & movementStartLabels < trialEndLabels(trialNum);
                thisMovementEndIndex =  trialStartLabels(trialNum) < movementEndLabels ...
                    & movementEndLabels < trialEndLabels(trialNum);
                if length(thisMovementStartIndex) ~= length(thisMovementEndIndex)
                    error("Number of Movement onset labels is not equal to movement end labels");
                end
                movementStarts(trialNum) = timestamp(movementStartLabels(trialNum));
                movementEnds(trialNum) = timestamp(movementEndLabels(trialNum));
            end            
            allData.(sF).timeGazeExpStart = expStart;  
            allData.(sF).timeGazeStart = trialStarts;  
            allData.(sF).timeGazeEnd = trialEnds;     
            allData.(sF).timeMoveStart = movementStarts;  
            allData.(sF).timeMoveEnd = movementEnds;      
            handTrialTimeStarts = allData.(sF).timeStart;
            handTrialTime0 = allData.(sF).time0;
     
            numTrials = length(trialStarts);
            for trialNum = 1:numTrials
                thisStartTime = trialStarts(trialNum);
                thisEndTime = trialEnds(trialNum);                
                    theseGazeIndexes = find( ...    
                        gazeDataFile.gaze_timestamp > thisStartTime & ...  
                        gazeDataFile.gaze_timestamp < thisEndTime & ...
                        gazeDataFile.confidence > confidenceThreshold);
                    theseGazeT = gazeDataFile.gaze_timestamp(theseGazeIndexes); 
                    theseGazeX = gazeDataFile.x_norm(theseGazeIndexes);
                    theseGazeY = gazeDataFile.y_norm(theseGazeIndexes);
                    badSample = find(diff(theseGazeT) <= 0);
                    while ~isempty(badSample)
                        allSample = 1:length(theseGazeT);
                        keepSample = setdiff(allSample, badSample+1);
                        theseGazeT = theseGazeT(keepSample);
                        theseGazeX = theseGazeX(keepSample);
                        theseGazeY = theseGazeY(keepSample);
                        badSample = find(diff(theseGazeT) <= 0);
                    end                    
                %if ~isempty(theseGazeIndexes)
                if length(theseGazeIndexes) > 1                  
                    theseGazeTZeroed = (theseGazeT - thisStartTime) - (handTrialTime0(trialNum) - handTrialTimeStarts(trialNum));
                    t = tVals;
                    x = interp1(theseGazeTZeroed, theseGazeX, tVals, 'makima', NaN); 
                    y = interp1(theseGazeTZeroed, theseGazeY, tVals, 'makima', NaN);                                
                    allData.(sF).gaze(:,trialNum,coord.x) = (x-x0Eye)*xUnitEye;  
                    allData.(sF).gaze(:,trialNum,coord.y) = (y-y0Eye)*yUnitEye;                    
                    if ~isempty(fixDataFile) % this is because some subjects
                        % have gaze files but no fixation files
                        % So then we fill the data structure with nans
                        %
                        % allData.(sF).fix(:,trialNum,coord.x) = nan(size(t));
                        % allData.(sF).fix(:,trialNum,coord.y) = nan(size(t));
                        %
                        % else                        
                        theseFixIndexes = find( ...
                            fixDataFile.start_timestamp+fixDataFile.duration > thisStartTime & ...
                            fixDataFile.start_timestamp < thisEndTime & ...
                            fixDataFile.dispersion < dispersionThresh );
                        theseFixStart = fixDataFile.start_timestamp(theseFixIndexes);
                        theseFixEnd = theseFixStart + fixDataFile.duration(theseFixIndexes)/1000; % duration data is in ms
                        theseFixX = fixDataFile.norm_pos_x(theseFixIndexes);
                        theseFixY = fixDataFile.norm_pos_y(theseFixIndexes);
                        numFix = length(theseFixIndexes);
                        
                        theseFixStartZeroed = (theseFixStart - thisStartTime) - (handTrialTime0(trialNum) - handTrialTimeStarts(trialNum));
                        theseFixEndZeroed = (theseFixEnd - thisStartTime) - (handTrialTime0(trialNum) - handTrialTimeStarts(trialNum));
                        
                        x = NaN*zeros(size(t));
                        y = NaN*zeros(size(t));
                        for fixNum = 1:numFix
                            x(theseFixStartZeroed(fixNum) <= t & t <= theseFixEndZeroed(fixNum)) = theseFixX(fixNum);
                            y(theseFixStartZeroed(fixNum) <= t & t <= theseFixEndZeroed(fixNum)) = theseFixY(fixNum);
                        end
                        allData.(sF).fix(:,trialNum,coord.x) = (x-x0Eye)*xUnitEye; % from fixations_on_surface_screen file
                        allData.(sF).fix(:,trialNum,coord.y) = (y-y0Eye)*yUnitEye;
                        
                    end
                end
                
                
            end
        end
    end    
    % Zero-correct the eye origin for origin location    
    originMin = -1.5; % -1.5
    originMax = -1.0; % -1.0
    method = 'rloess'; 
    span = 25;  
    tOr = find(originMin < tVals & tVals < originMax);  
    numRows = ceil(sqrt(numSubj));
    numCols = ceil(numSubj / numRows);

    for subjNum = 1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;        
        if isfield(allData.(sF), "gaze")
            if ~isfield(allData.(sF), "steps")
                allData.(sF).steps = struct();
            end
            if ~isfield(allData.(sF).steps, "read") && isfield(allData.(sF),"fix" )
                allData.(sF).steps.read = struct( ...
                    "gaze", allData.(sF).gaze, ...
                    "fix", allData.(sF).fix);
                g = allData.(sF).steps.read.gaze;
                f = allData.(sF).steps.read.fix;
                % Get the median position of the eyes at origin fixation
                gOr = squeeze(median(g(tOr,:,:), 1));  
                % And then smooth it across trials
                sgOr = cat(2, ...
                    smooth(gOr(:,coord.x), span, method), ...
                    smooth(gOr(:,coord.y), span, method));                 
                % Then shift the data so the smoothed 0 is 0
                allData.(sF).gazeOr = sgOr;
                allData.(sF).gaze = g - shiftdim(gOr, -1);   
                allData.(sF).fix = f - shiftdim(sgOr, -1);                
                if length(sgOr) > 160
                    sgOr = sgOr([1:140 221:240],:);
                    gOr = gOr([1:140 221:240],:);
                    allData.(sF).gazeOr = sgOr;  
                    allData.(sF).gaze = allData.(sF).gaze(:,[1:140 221:240],:);
                    allData.(sF).fix = allData.(sF).fix(:,[1:140 221:240],:);
                end
                gShift = squeeze(median(allData.(sF).gaze(tOr,:,:), 1));

            else
                allData.(sF).steps.read = struct( ...
                    "gaze", allData.(sF).gaze);
                g = allData.(sF).steps.read.gaze;
                % Get the median position of the eyes at origin fixation
                gOr = squeeze(median(g(tOr,:,:), 1));
                % And then smooth it across trials
                sgOr = cat(2, ...
                    smooth(gOr(:,coord.x), span, method), ...
                    smooth(gOr(:,coord.y), span, method));
                
                % Then shift the data so the smoothed 0 is 0
                allData.(sF).gazeOr = sgOr;
                allData.(sF).gaze = g - shiftdim(gOr, -1);   
                if size(sgOr,1) > 160
                    sgOr = sgOr([1:140 221:240],:);
                    gOr = gOr([1:140 221:240],:);
                    allData.(sF).gazeOr = sgOr;
                    allData.(sF).gaze = allData.(sF).gaze(:,[1:140 221:240],:);
                end
                gShift = squeeze(median(allData.(sF).gaze(tOr,:,:), 1));
            end
     
        else
            text(0.5, 0.5, "no eye data", 'Units', 'normalized');
        end     
    end    
    % Rotate all movements to vector origin to target 
    
    for subjNum = 1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;
       
        if ~isfield(allData.(sF), "steps")
            if isfield(allData.(sF), "gaze")
                error("Should align hand and eye before rotating");
            end
            allData.(sF).steps = struct();
        end
        if ~isfield(allData.(sF).steps, 'unrotated')
            allData.(sF).steps.unrotated = struct('hand', allData.(sF).hand); 
        end
        if isfield(allData.(sF), 'gaze') % && isfield(allData.(sF), 'fix')
            allData.(sF).steps.unrotated.gaze = allData.(sF).gaze;
            if isfield(allData.(sF), 'fix')
                allData.(sF).steps.unrotated.fix = allData.(sF).fix;
            end
        end
        allData.(sF).steps.unrotated.target = allData.(sF).target;     
        h = allData.(sF).steps.unrotated.hand;      
        or = squeeze(median(h(tOr,:,:), 1));      
        tgtDir = allData.(sF).target - or;
        tgtDir = tgtDir ./ vecnorm(tgtDir,2,2);
        perpDir = [tgtDir(:,2) -tgtDir(:,1)];        
        tar = allData.(sF).steps.unrotated.target;
        tarRotated = [
            dot(tar, perpDir, 2) ...
            dot(tar, tgtDir, 2) ...
            ];
        allData.(sF).target = tarRotated;    
               
        hZeroed = h - shiftdim(or, -1);
        hRotated = cat(3, ...
            dot(hZeroed, repmat(shiftdim(perpDir,-1), [numTVals 1 1]),3), ...
            dot(hZeroed, repmat(shiftdim(tgtDir,-1), [numTVals 1 1]),3) ...
            );
        allData.(sF).hand = hRotated;            
        if isfield(allData.(sF), 'gaze') % && isfield(allData.(sF), 'fix')
            if size(allData.(sF).steps.unrotated.gaze,2) > 160
                g = allData.(sF).steps.unrotated.gaze(:,[1:140 221:240],:);
                if isfield(allData.(sF), 'fix')
                    f = allData.(sF).steps.unrotated.fix(:,[1:140 221:240],:);
                end
            else
                gZeroed = allData.(sF).gaze;
                gRotated = cat(3,...
                    dot(gZeroed, repmat(shiftdim(perpDir,-1), [numTVals 1 1]),3), ...
                    dot(gZeroed, repmat(shiftdim(tgtDir,-1), [numTVals 1 1]),3) ...
                    );
%                 if isfield(allData.(sF), 'fix')
%                     fZeroed =  f - shiftdim(or, -1);
%                     fRotated = cat(3,...
%                         dot(fZeroed, repmat(shiftdim(perpDir,-1), [numTVals 1 1]),3), ...
%                         dot(fZeroed, repmat(shiftdim(tgtDir,-1), [numTVals 1 1]),3) ...
%                         );
%                     allData.(sF).fix = fRotated;
%                 end
                allData.(sF).gaze = gRotated;  
            end
        end
    end    
    % Find beginning and end of movements    
    for subjNum = 1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;        
        if size(allData.(sF).hand,2) > 160
            h = allData.(sF).hand(:,[1:140 221:240],:);
        else 
            h = allData.(sF).hand;
            % Take derivative using the magic of smooth
            smoothH = smoothdata(h, 'loess', 15, 'omitnan');
            dH = diff(smoothH) / dt;
            dH(end+1,:,:) = 0; % Add one more time step to make up for the one that diff took away
            spd = sqrt( sum(dH.^2,3) );
            peakSpd = max(spd);
            thresh = 0.10*peakSpd;
            onset = zeros(numTrials, 1);
            offset = zeros(numTrials, 1);
            for trialNum = 1:numTrials
                thisSpd = spd(:,trialNum);
                thisPeakSpd = peakSpd(trialNum);
                thisThresh = thresh(trialNum);
                halfPeakTime = tVals(find(tVals > 0 & thisSpd > thisPeakSpd/2, 1, 'first'));            
                if isempty(halfPeakTime)
                    lastIndexUnderThresh = nan;
                    onset(trialNum) = nan;
                else
                    lastIndexUnderThresh = find( tVals < halfPeakTime & thisSpd < thisThresh, 1, 'last');
                    onset(trialNum) = tVals(lastIndexUnderThresh+1);
                end
                if isempty(halfPeakTime)
                    firstIndexUnderThresh = nan;
                    offset(trialNum) = nan;
                else
                    firstIndexUnderThresh = find( tVals > halfPeakTime & thisSpd < thisThresh, 1);
                    offset(trialNum) = tVals(firstIndexUnderThresh);
                end
            end
            allData.(sF).onset = onset;  
            allData.(sF).offset = offset;
        end
        
    end


end

%% Reaction time during the experiment
if RT
    rT = nan(numTrials-10);
    for subjNum =  1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;
    rT =  allData.(sF).timeMoveStart(11:160)' - allData.(sF).targetOnset(11:160); 
    rT(rT < -5) = NaN;
    allData.(sF).reactionTime = rT;
%     rT = smoothdata(rT,"movmean");
    figure();
    plot(rT,'color','#ED5736','LineWidth',1);
    ylabel('[s]')
    xlabel('Trial #')
    sgtitle('Reaction time of subject ' + subjID)
    ylim([0 4])
    hline(1,'k--')
    vline([50 130],'k--')
    grid on;

subjFolder = fullfile(outputFolder, subjID);
mkdir(subjFolder);
filename = sprintf('%s_RT.png',subjID);
fullFilename = fullfile(subjFolder, filename);
saveas(gcf, fullFilename);
    end
end
%% Hand and eye error 
if HandEyeErr
errorC = nan(numSubj, numTrials-10);    
for subjNum = 1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;      
        group(subjNum) = allData.(sF).group;        
        if group(subjNum) == 0            
            delay = 0.0;
        elseif group(subjNum) == 1
            delay = 0.1;
        else
            delay = 1.2;
        end        
        offset = zeros(numTrials, 1);
        EndAngle = zeros(numTrials, 1);
        for trialNum = 1:numTrials
            offset = allData.(sF).endTime(trialNum);
            onset = allData.(sF).time0(trialNum);
            offsetTrial = offset - onset;
            endIdx = find(offsetTrial <= tVals,1);
            handEndX = allData.(sF).hand(endIdx,trialNum,1);
            handEndY = allData.(sF).hand(endIdx,trialNum,2);
            tx = allData.(sF).target(trialNum,1);
            ty = allData.(sF).target(trialNum,2);           
            
            [TargetTheta,TargetRadius] = cart2pol(tx,ty); 
            [EndTheta,EndRadius] = cart2pol(handEndX,handEndY);
            if isempty(EndTheta)
                EndAngle(trialNum) = nan;
            else
                EndAngle(trialNum) = wrapTo180(rad2deg(TargetTheta - EndTheta));
            end

        end        
        allData.(sF).moveEndError = EndAngle;

        enderror = allData.(sF).moveEndError(11:end); 
        mBase = mean(enderror(1:49),'omitnan');  
        errorCorrect =  enderror - mBase; 
        % interpolate NaN values
        nanIdx = find(isnan(errorCorrect)); % find NaN values
        tvals = 1:length(errorCorrect); % create time vector
        interpVals = interp1(tvals(~isnan(errorCorrect)), errorCorrect(~isnan(errorCorrect)), tvals(nanIdx), 'spline'); % interpolate NaN values
        errorCorrect(nanIdx) = interpVals; % replace NaN values with interpolated values
       
        errorC(subjNum,:) = errorCorrect;
        allData.(sF).errorCorrect = errorCorrect;
end
 
% Eye error
    numTrials = 160;    
    AngleAtTime = nan(numSubj,numTrials);
    group = nan(numSubj,1);
    
    for subjNum =  1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;
        
        group(subjNum) = allData.(sF).group;     

        if group(subjNum) == 0
            delay = 0.0;
        elseif group(subjNum) == 1
            delay = 0.1;
        else
            delay = 1.2;
        end
        
        if ~isfield(allData.(sF),'gaze')
            continue
        end
        
        for trialNum = 1:numTrials
            
            gazeX = allData.(sF).steps.unrotated.gaze(:,trialNum,1);
            gazeY = allData.(sF).steps.unrotated.gaze(:,trialNum,2);
            
            tx = allData.(sF).steps.unrotated.target(trialNum,1);
            ty = allData.(sF).steps.unrotated.target(trialNum,2);
                       
            [gazeTheta, gazeRadius] = cart2pol(gazeX,gazeY);
            [tht,rhot] = cart2pol(tx,ty);  

            % furthest fixation (in terms of theta) that is more than ~ 50% away from origin
            %         idxAtDist = find(fixRadius > 0.4*rhot);
            idxAtDist = find(gazeRadius > 0.2*rhot & gazeRadius < rhot);      
            thetaAtDist = gazeTheta(idxAtDist);                              
            [maxTheta, idxAtMaxTheta] = max(thetaAtDist);
            [minTheta, idxAtMinTheta] = min(thetaAtDist);
            
            if abs(maxTheta) > abs(minTheta)                                  
                radiusAtMaxTheta = gazeRadius(idxAtDist(idxAtMaxTheta));
            else
                radiusAtMaxTheta = gazeRadius(idxAtDist(idxAtMinTheta));
            end

            % taking the median of the fixations during movement onset 
            thetaTime = median(gazeTheta(245:255));                               
            % discarding eye movements that are further than 45 degrees frm
            % the target area
   
            if rad2deg(tht - thetaTime) >= 35 || rad2deg(tht - thetaTime) <= -35     
                angleTime(trialNum,:) = nan;
            end
            angleTime(trialNum,:) = wrapTo180(rad2deg(tht - thetaTime));  
                                                                                     
        end
        
        % remove outliers (everything above a deviation of 3x the median)
        [~,IDXOut] = rmoutliers(angleTime,'movmean',8);  
        angleTime(IDXOut) = NaN;                        

        
        AngleAtTime(subjNum,:) = smoothdata(angleTime, 'loess', 12,'omitnan');  
%         AngleAtTime(subjNum,:) = smoothdata(angleTime, 'sgolay', 10, 'omitnan');    
        allData.(sF).eyeAngleAtMO = AngleAtTime(subjNum,11:160);  

    figure();
    Herr = plot(errorCorrect,'LineWidth',1);
    hold on
    Eerr = plot(allData.(sF).eyeAngleAtMO,'LineWidth',1);
    ylabel('[°]')
    xlabel('Trial #')
    sgtitle('Hand and eye error across all trials of subject ' + subjID)
    ylim([-60 30])
    xlim([1 150])
    hline(0,'k:')
    hline(-30,'k:')
    vline(50,'k--')
    vline(130,'k--')
    legend([Herr,Eerr],'Hand','Eye')
    grid on;
    set(gca,'GridLineStyle','--')
    
    subjFolder = fullfile(outputFolder, subjID);
    mkdir(subjFolder);
    filename = sprintf('%s_HandEyeError.png',subjID);
    fullFilename = fullfile(subjFolder, filename);
    saveas(gcf, fullFilename);

        end
    end
%% Trace plot
if Hand_trace
    numShowTrialsB = 60;
    numShowTrialsEA = 70;
    numShowTrialsLA = 140;
    numShowTrialsAE = 150;
    numShowTrialsAEL = 160;
    
    selectedTrialsB = 11:5:numShowTrialsB;
    selectedTrialsEA = 61:1:numShowTrialsEA;
    selectedTrialsLA = 71:10:numShowTrialsLA;
    selectedTrialsAE = 141:1:numShowTrialsAE;
    selectedTrialsAEL = 151:1:numShowTrialsAEL;

    for subjNum =  1:numSubj
        subjID = allSubjs(subjNum);
        sF = "s"+subjID;    
    figure();
    set(gcf, 'position', [0 0 1200 800]);

    %% BASELINE
    ax1 = axes;
    set(ax1, 'position', [0.13 0.11 0.775 0.815]);
    for i = 1:length(selectedTrialsB)
        thisTrialNum = selectedTrialsB(i);
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);
        
        B = plot(hx, hy, 'color','#808080', 'LineWidth', 1); %#ffff4d
        hold on;
    end
    axis([-15 5 -2 20]);
    %% Aftereffect  
    ax2 = axes;
    set(ax2, 'position', [0.13 0.11 0.775 0.815]);
    % Initialize an array to store the custom colors for Aftereffect trials
    customColorsAE = zeros(numShowTrialsAE - 141 + 1, 3);

    for i = 1:length(selectedTrialsAE)
        thisTrialNum = selectedTrialsAE(i);
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);
        
        % caculate the value of color of current trial
        % map color value to (0,1)
        value = (thisTrialNum - 141) / (numShowTrialsAE - 141);
        % R、G、B from 0 to 1
        thisColor = [1, value, 0];
        % Store the custom color in the array
        customColorsAE(i, :) = thisColor;

        AE = plot(hx, hy, 'color',thisColor, 'LineWidth', 1); %#77AC30
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;

    %% ADAPTATION (early)
    ax3 = axes;
    set(ax3, 'position', [0.13 0.11 0.775 0.815]);
    % Initialize an array to store the custom colors for Aftereffect trials
    customColorsEA = zeros(numShowTrialsEA - 61 + 1, 3);
    for i = 1:length(selectedTrialsEA)
        thisTrialNum = selectedTrialsEA(i);
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);

        % caculate the value of color of current trial
        % map color value to (0,1)
        value = (thisTrialNum - 61) / (numShowTrialsEA - 61);
        % R、G、B from 0 to 1
        thisColor = [0, value, 1];
        % Store the custom color in the array
        customColorsEA(i, :) = thisColor;

        EA = plot(hx, hy, 'color',thisColor, 'LineWidth', 1); % Use a thicker blue line for late adaptation trials
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;
        %% ADAPTATION (late)
    ax4 = axes;
    set(ax4, 'position', [0.13 0.11 0.773 0.813]);
    for i = 1:length(selectedTrialsLA) % Consider the last 40 trials as early adaptation
        thisTrialNum = selectedTrialsLA(i);
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);
        LA = plot(hx, hy, 'color',"#696969", 'LineWidth', 1,'LineStyle','--'); %#00FFFF
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;
        %% aftereffect (late 10 trials)
    ax5 = axes;
    set(ax5, 'position', [0.1 0.11 0.773 0.813]);
    for i = 1:length(selectedTrialsAEL) % Consider the last 40 trials as early adaptation
        thisTrialNum = selectedTrialsAEL(i);
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);
        AEL = plot(hx, hy, 'color',"#696969", 'LineWidth', 1,'LineStyle',':'); %#00FFFF
        hold on;
    end
    axis([-15 5 -2 20]);
    axis off;
    scatter(0, 12, 100, 'k','filled')
    scatter(0, 0, 100, 'k','filled')
     
    % Set the custom colormap for earlyadapt
    colormap(ax3, customColorsEA);
    % Add colorbar for Aftereffect
    AColorbar = colorbar(ax3, 'Location', 'west');  
    % Modify colorbar ticks and labels
    tickValues = linspace(selectedTrialsEA(1), selectedTrialsEA(end), 10); % You can adjust the number of ticks as needed
    tickLabels = arrayfun(@(x) sprintf('%d', round(x)), tickValues, 'UniformOutput', false);
    AColorbar.Ticks = (tickValues - selectedTrialsEA(1)) / (selectedTrialsEA(end) - selectedTrialsEA(1));
    AColorbar.TickLabels = tickLabels;
    AColorbar.Label.String = 'trialNum-Early adapt';

    % Set the custom colormap for Aftereffect
    colormap(ax2, customColorsAE);
    % Add colorbar for Aftereffect
    hColorbar = colorbar(ax2, 'Location', 'east');  
    % Modify colorbar ticks and labels
    tickValues = linspace(selectedTrialsAE(1), selectedTrialsAE(end), 20); % You can adjust the number of ticks as needed
    tickLabels = arrayfun(@(x) sprintf('%d', round(x)), tickValues, 'UniformOutput', false);
    hColorbar.Ticks = (tickValues - selectedTrialsAE(1)) / (selectedTrialsAE(end) - selectedTrialsAE(1));
    hColorbar.TickLabels = tickLabels;
    hColorbar.Label.String = 'trialNum-washout';

hold off;
% axis equal;
xlim([-15 4]);
ylim([-2 20]);
legend([B, LA, AEL], 'Baseline', 'Late adaptation', 'Late washout',Location='south');
title('Hand trace of subject ' + subjID);
xlabel("x (cm)");
ylabel("y (cm)");

subjFolder = fullfile(outputFolder, subjID);
mkdir(subjFolder);
filename = sprintf('%s_HandTrace.png',subjID);
fullFilename = fullfile(subjFolder, filename);
saveas(gcf, fullFilename);
    end
end
%%
if Eye_hand_trace
numShowTrialsB = 60;
numShowTrialsA = 140;
numShowTrialsAE = 160;
for subjNum =  1:numSubj
    subjID = allSubjs(subjNum);
    sF = "s"+subjID;
    figure('Units', 'pixels', 'Position', [100, 100, 800, 600]);
    i=1;
    % BASELINE
    for thisTrialNum = [20 30 40 50 60] %[20 30 40 50 60]% 11:numShowTrialsEB   
        subplot(2,10,i)
        i=i+1;
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);

        gx = allData.(sF).gaze(:,thisTrialNum,coord.x);
        gy = allData.(sF).gaze(:,thisTrialNum,coord.y);

        nonNan = ~isnan(gx) & ~isnan(gy);
        B = plot(gx, gy, 'color', "#77AC30");
        hold on;
        H = plot(hx, hy, 'color', 'k','LineWidth',1);
        scatter(0, 12, 40, 'k','filled')
        scatter(0, 0, 40, 'k','filled')
        xlim([-20 10]);
        ylim([-2 20]);
        title(sprintf('trial %d',thisTrialNum));
        hold on;
        if i == 2
            ylabel("y (cm)");
        end
    end
    % AE
    for thisTrialNum = [142 145 150 155 160]  % [142 145 150 155 160] % 141:numShowTrialsAE
        subplot(2,10,i)
        i=i+1;
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);
        gx = allData.(sF).gaze(:,thisTrialNum,coord.x);
        gy = allData.(sF).gaze(:,thisTrialNum,coord.y);

        nonNan = ~isnan(gx) & ~isnan(gy);
        AE = plot(gx, gy, 'color',"#A2142F");
        hold on;
        H = plot(hx, hy, 'color', 'k','LineWidth',1);
        scatter(0, 12, 40, 'k','filled')
        scatter(0, 0, 40, 'k','filled')
        xlim([-20 10]);
        ylim([-2 20]);
        title(sprintf('trial %d',thisTrialNum));
        hold on;
    end
    % ADAPTATION
    for thisTrialNum = [61 70 80 90 100 110 120 125 130 140]  % [61 70 80 90 100 110 120 130 140] % 61:numShowTrialsA
        subplot(2,10,i)
        i=i+1;
        hx = allData.(sF).hand(:,thisTrialNum,coord.x);
        hy = allData.(sF).hand(:,thisTrialNum,coord.y);
        gx = allData.(sF).gaze(:,thisTrialNum,coord.x);
        gy = allData.(sF).gaze(:,thisTrialNum,coord.y);
        nonNan = ~isnan(gx) & ~isnan(gy);
        A = plot(gx(nonNan), gy(nonNan), 'color',"#0032BD");
        hold on;
        H = plot(hx, hy, 'color', 'k','LineWidth',1);
        scatter(0, 12, 40, 'k','filled')
        scatter(0, 0, 40, 'k','filled')
        xlim([-20 10]);
        ylim([-2 20]);
        title(sprintf('trial %d',thisTrialNum));
        hold on;
    end

legend([H, B, AE, A], 'Hand', 'Eye Baseline', 'Eye AfterEffect', 'Eye Adaptation', 'Location','northoutside');
sgtitle('Eye and hand trace in sampling trials of subject ' + subjID);
xlabel("x (cm)");
% figure size
figWidth = 1500; % width
figHeight = 800; % height
set(gcf, 'Position', [50, 50, figWidth, figHeight]);

subjFolder = fullfile(outputFolder, subjID);
mkdir(subjFolder);
filename = sprintf('%s_EyeHandTrace.png',subjID);
fullFilename = fullfile(subjFolder, filename);
saveas(gcf, fullFilename);
end
end