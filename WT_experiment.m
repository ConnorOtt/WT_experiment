%--------------------------------------------------------------------------
% This class is built to take in, sort, and average data files from the 
% ITLL wind tunnel at the University of Colorado at Boulder. 
%
%  Properties (immutable) - Can only be modified in constructor
%       
%           atmPressure: Atmospheric Pressure Readings 
%               atmTemp: Atmospheric Temperature
%            atmDensity: Atmospheric Density
%                 V_inf: Free Stream Velocity
%          pitotDynamic: Free Stream(?) Dynamic Pressure 
%            auxDynamic: Probe(?) Dynamic Pressure
%            scanivalve: Static Port Pressures - (Matrix, holds all ports)
%                   AoA: Angle of attack of body on sting
%            stingPitch: (can't remember what this measures)    
%           stingNormal: Sting Balance Normal Force Measurements
%            stingAxial: Sting Balance Axial Force Measurements
%                  ELDx: x-position of ELD Probe
%                  ELDy: y-position of ELD probe
%              fullData: Entire data set before parsing and averaging
%         dataCalibrate: Data gathered from sting calibration tests, will
%                        only be populated if a zero-free stream test is
%                        detected. If populated, calibration data is also
%                        removed from the rest of the data. 
%
%  Propeties (Dependent) - Dependent on other fields in WT_experiment.
%  
%                  lift: Lift values - calibrated w/ zero V trial
%                  drag: Drag values - calibrated
%               V_pitot: Velocity determined by pitot probe(?)
%              C_pPorts: Pressure coefficients of data gathered with 
%                        scanivalve measurements
%
%  Properties - Can be modified anywhere
%
%            sampleSize: Number of samples sent to file (per 'click')
%              testName: Name of data file
%
%  Methods
%           obj = WT_experiment(dataSet, sampleSize) parses, averages, and,
%               if necessary, separates calibration data of matrix 
%               'dataSet' with sample size 'sampleSize'.
%           D = get.drag(obj) Supposedly should calculated drag when the
%               drag property is requested. It's calculating it during
%               construction so maybe I'm misunderstanding that. 
%           L = get.lift(obj) "" - for lift.
%           V = get.V_pitot(obj) "" - for velocity at pitot probe
%               measurement. I actually don't remember if this is the
%               useful velocity or not.
%           C_p = get.C_pPorts(obj) "" - Pressure coefficient from
%               scanivalve port measurements. 
%           ObjMean = mean(tempObj, objCell) calculates the mean of a cell
%               array of WT_experiment objects. 
%               
%  Methods (Static)
%           [] = plot(varargin) is possible a plot wrapper, may trash it,
%               may not.
%
%  Vocab hierarchy that I'm really gonna try to stick to:
%       Experiment > Test > Trial == Sample > Data Point
%   `   - There are sampleSize samples in a trial.
%       - There are length(fullData) / sampleSize trials in a test.
%       - There are probably like 2 or 3 tests in an experiment.
%
% Created: 11/3/17 - Connor Ott
% Last Modified: 11/10/17 - Connor Ott
%--------------------------------------------------------------------------

classdef WT_experiment
    properties (SetAccess = immutable)
        atmPressure
        atmTemp
        atmDensity
        V_inf
        pitotDynamic
        auxDynamic
        scanivalve
        AoA
        stingPitch
        stingNormal
        stingAxial
        ELDx
        ELDy
        fullData
        dataCalibrate
    end
    properties (Dependent)
       lift
       drag
       V_pitot
       C_pPorts
    end
    properties 
        sampleSize
        testName = '';
    end
    methods
        % Constructor
        function obj = WT_experiment(dataSet, sampleSize)
            % Constructor - Goes straight to mean
           
            if nargin > 0
                obj.fullData = dataSet;
                obj.sampleSize = sampleSize;
                [r, ~] = size(dataSet);
                numTrials = r/sampleSize;
               
                % Allocation
                [obj.atmPressure, ...
                 obj.atmTemp, ...
                 obj.atmDensity, ...
                 obj.V_inf, ...
                 obj.pitotDynamic, ...
                 obj.auxDynamic, ...
                 obj.scanivalve, ...
                 obj.AoA, ...
                 obj.stingNormal, ...
                 obj.stingAxial, ...
                 obj.stingPitch, ...
                 obj.ELDx, ...
                 obj.ELDy] = deal(zeros(numTrials, 1));
                 
                 obj.scanivalve = zeros(numTrials, 16);
                
                % Parsing and averaging data between trials
                for i = 1:numTrials
                    tempTrial = dataSet(1 + sampleSize*(i-1):...
                        sampleSize*(i), :);
                    
                    obj.atmPressure(i)    = mean(tempTrial(:, 1));
                    obj.atmTemp(i)        = mean(tempTrial(:, 2));
                    obj.atmDensity(i)     = mean(tempTrial(:, 3));
                    obj.V_inf(i)          = mean(tempTrial(:, 4));
                    obj.pitotDynamic(i)   = mean(tempTrial(:, 5));
                    obj.auxDynamic(i)     = mean(tempTrial(:, 6));
                    obj.AoA(i)            = mean(tempTrial(:, 23));
                    obj.stingNormal(i)    = mean(tempTrial(:, 24));
                    obj.stingAxial(i)     = mean(tempTrial(:, 25));
                    obj.stingPitch(i)     = mean(tempTrial(:, 26));
                    obj.ELDx(i)           = mean(tempTrial(:, 28));
                    obj.ELDy(i)           = mean(tempTrial(:, 27));
                    obj.scanivalve(i, :)  = mean(tempTrial(:, 7:22));
                end
                
                % Attempts to detect a sting calibration trial,
                % characterized by a zero (low) free stream velocity
                
                % Assumes that if a calibration is present, it is the first
                % experiment in the data set. This is typical for ITLL lab
                % procedures. 
                logV = obj.V_inf < 2;
                zeroVidx = find(logV);
                % if detected, separate callibration.
                
                propNames = fieldnames(obj);
                numSep = 13; % Number of things to separate
                if any(logV) 
                    cutOff = zeroVidx(end);
                    for i = 1:numSep
                        name = propNames{i};
                        [~, c] = size(obj.(name));
                        if c == 1 % checks any matrices
                            obj.dataCalibrate.(name) = ...
                                                 obj.(name)(1:cutOff);
                            obj.(name) = obj.(name)(cutOff + 1:end);
                        end
                    end
                    obj.dataCalibrate.scanivalve = ...
                                           obj.scanivalve(1:cutOff, :);
                    obj.scanivalve = obj.scanivalve(cutOff + 1:end, :);
                end
                
            end
        end
        % get drag - calculates drag from input data
        function D = get.drag(obj)
            
            if isempty(obj.dataCalibrate)
                error(['No sting balance calibrations found in' ...
                       ' data set.']);
            else
                % Calibration data
                N_cal = obj.dataCalibrate.stingNormal;
                A_cal = obj.dataCalibrate.stingAxial;
                
                % Length of calibration test should equal length of all
                % tests in dataSet. (If performed correctly).
                cutOff = length(N_cal); 
                numTest = length(obj.stingNormal) / cutOff;
                
                for i = 1:numTest
                    k = cutOff * i;
                    j = 1 + cutOff * (i - 1);
                    D(j:k) = (obj.stingNormal(j:k) - N_cal) .* ...
                                                 sind(obj.AoA(j:k)) + ...
                             (obj.stingAxial(j:k) - A_cal) .* ...
                                                 cosd(obj.AoA(j:k));
                end
                D = D';
            end
        end
        % get lift - calculates lift from input data
        function L = get.lift(obj)
            if isempty(obj.dataCalibrate)
                error(['No sting balance calibrations found in' ...
                       ' data set']);
            else
                % Callibration values from 0 m/s tests.
                N_cal = obj.dataCalibrate.stingNormal;
                A_cal = obj.dataCalibrate.stingAxial;
                
                % Length of calibration test should equal length of all
                % tests in dataSet. (If performed correctly).
                cutOff = length(N_cal); 
                numTests = length(obj.stingNormal) / cutOff;
                
                for i = 1:numTests
                    k = cutOff * i;
                    j = 1 + cutOff * (i - 1);
                    L(j:k) = (obj.stingNormal(j:k) - N_cal) .* ...
                                                 cosd(obj.AoA(j:k)) - ...
                             (obj.stingAxial(j:k) - A_cal) .* ...
                                                 sind(obj.AoA(j:k)); 
                end
                L = L';
            end
        end
        % get Velocity at pitot (good to have, m8)
        function V = get.V_pitot(obj)
            q_inf = obj.pitotDynamic;
            rho_inf = obj.atmDensity;
            V = (2 * q_inf ./ rho_inf).^(1/2);
        end
        % get pressure coefficient at static ports on airfoil
        function C_p = get.C_pPorts(obj)
            q_inf = mean(obj.pitotDynamic);
            P_static = obj.scanivalve; % Pressure difference
            [r, ~] = size(P_static);
            
            C_pRaw = P_static ./ q_inf;
            % Interpolation for the "trailing edge port" (possibly)
            C_p = [C_pRaw(:, 1:8), zeros(r, 1), C_pRaw(:, 9:16)];
        end
        % Take the mean of some WT_experiments. 
        function objMean = mean(tempObj, objCell)
            % To use a method need an input of the type of your
            % class, but I want this to work for cell arrays of my class so
            % 'tempObj' calls my mean function but my mean function doesn't
            % need its stupid shit. 
            
            baseObj = objCell{1}; % Construct mean obj based on this one
            numObj = numel(objCell);
            [rFull, cFull] = size(baseObj.fullData);
            dataSet = zeros(rFull, cFull);
            sSize = baseObj.sampleSize;
            
            for i = 1:numObj
                if ~isempty(objCell{i})
                    dataSet = dataSet + objCell{i}.fullData;
                else
                    numObj = numObj - 1;
                end
            end
            dataSet = dataSet / numObj;
            
            objMean = WT_experiment(dataSet, sSize);
        end
       
    end
    methods (Static)
        function [] = plot(obj, varargin)
            % This is hard, and maybe not as useful as I thought it would
            % be
        end
    end
end

