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
%   `   - There are sampleSize data points in a trial (also called sample).
%       - There are an arbitrary number of trials in a test, probably on
%         order of 10 though (think changing angle of attack).
%       - There are probably like 2 or 3 tests in an experiment (think 
%         changing velocities).
%       The last 2 can change from experiment to experiment, but there will
%       always be sampleSize data points in a trial (or sample).
%
%
% Created: 11/3/17 - Connor Ott
% Last Modified: 11/13/17 - Connor Ott
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
        % Parse or simply separate multiple or just one test from the rest.
        function objCell = testParse(obj, parameterName, pVec)
            %--------------------------------------------------------------
            % Inputs
            %      obj -       object getting parsed
            %      parameter - String of parameter name. 
            %      pVec -      Vector of parameter values by which function
            %                  does the parsing.
            % Outputs
            %      cellArray - cell array of separated WT_experiment
            %                  objects
            %--------------------------------------------------------------
            
            name = parameterName;
            parSepVec = obj.(name); % Parameter used for separation
            
            numTest = length(pVec);
            objCell = cell(numTest, 1);
            L_exp = length(parSepVec); % Length of entire experiment
            L_test = L_exp / numTest; % Lenth of a single test
            sSize = obj.sampleSize; % Sample size
            
            logMat = zeros(L_exp, numTest);
            for i = 1:numTest
                tempPVal = pVec(i);
                logMat(:, i) = parSepVec >= tempPVal*0.9 & ...
                               parSepVec <= tempPVal*1.1;
                
                p_idx = find(logMat(:, i));
                temp1 = p_idx(1);     % first index of given paramter
                tempEnd = p_idx(end); % last index of given parameter
                
                % Pulling out the full data set corresponding to the given
                % parameter. 
                
                % fullData includes calibration data if there is any, don't
                % want that getting mixed up in parsed data, hence the if 
                % statement. 
                if isempty(obj.dataCalibrate) 
                    m = 1 + ((temp1 - 1) * sSize);
                    n = (tempEnd * sSize);
                    parsedDataSet =  obj.fullData(m:n, :);
                else
                    m = 1 + (temp1 + L_test - 1) * sSize;
                    n = (tempEnd + L_test) * sSize;
                    parsedDataSet = obj.fullData(m:n, :);
                end
                
                objCell{i} = WT_experiment(parsedDataSet, sSize);
                objCell{i}.fileName = [obj.fileName, ' - ', name, ...
                                         ' = ', num2str(pVec(i))];
            end
            
        end
        % Combine two objects with odd anmatlad even angles of attack into a
        % single object. 
        function objZip = AoAzip(obj1, obj2)
           %---------------------------------------------------------------
           % To be used only with two objects wherein the only difference 
           % between the two is that one contains odd angles of attack and
           % the other contains the respective even angles of attack that
           % would fit in between the others for the same wind tunnel
           % conditions. 
           %
           % For instance, two groups take data at 25 m/s with the F-16
           % model in its dirty configuration. One group takes data at odd
           % angles of attack, the other at even angles of attack. This
           % function will sort things out and zip it up like a zipper
           % into a single object. 
           %
           % This will not work unless a scenario similar to the above is
           % present (just use your head)
           %---------------------------------------------------------------
           
           isOdd1 = mod(round(obj1.AoA(1)), 2);
           isOdd2 = mod(round(obj2.AoA(1)), 2);
           
           if isOdd1 == isOdd2
              error(['Both even and odd angle of attacked objects must',...
                     ' be present']); 
           else
               if isOdd1
                   oddObj = obj1;
                   evenObj = obj2;
               elseif isOdd2
                   oddObj = obj2;
                   evenObj = obj1;
               end
               
               minOdd = min(oddObj.AoA);
               minEven = min(evenObj.AoA);
               
               
               L_exp = length(oddObj.AoA); % Experiment length
               L_zipped = L_exp * 2; % Both experiments together
               fullDataOdd = oddObj.fullData;
               fullDataEven = evenObj.fullData;
               sSize = oddObj.sampleSize;
               [~, c] = size(fullDataOdd);
               fullDataTotal = zeros(L_zipped, c);
               
               
               for i = 1:L_zipped
                   % For populating fullDataTotal -  i = 1, 2, 3, 4, ...
                   m = 1 + (i-1) * sSize; 
                   n = i * sSize;
                   % For taking from even and odd - j = 1, 1, 2, 2, ...
                   j = ceil(i/2);
                   o = 1 + (j - 1) * sSize;
                   p = j * sSize; 
                   if minOdd < minEven % Start by pulling from odds
                       if mod(i, 2) == 1 % pull from odds.
                           fullDataTotal(m:n, :) = fullDataOdd(o:p, :);
                       else
                           fullDataTotal(m:n, :) = fullDataEven(o:p, :);
                       end
                   else % Start pulling from evens
                       if mod(i, 2) == 1 % pull from evens.
                           fullDataTotal(m:n, :) = fullDataEven(o:p, :);
                       else
                           fullDataTotal(m:n, :) = fullDataOdd(o:p, :);
                       end
                   end
               end
               
               objZip = WT_experiment(fullDataTotal, sSize);
           end
        end
    end
    methods (Static)
        function [] = plot(obj, varargin)
            % This is hard, and maybe not as useful as I thought it would
            % be
        end
    end
end

