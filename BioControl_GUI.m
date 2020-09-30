classdef platform_test_test < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        BioControl                      matlab.ui.Figure
        FileMenu                        matlab.ui.container.Menu
        OpenfolderMenu                  matlab.ui.container.Menu
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        StatusPanel                     matlab.ui.container.Panel
        DataacquisitionLampLabel        matlab.ui.control.Label
        DataacquisitionLamp             matlab.ui.control.Lamp
        ErrorlogTextAreaLabel           matlab.ui.control.Label
        ErrorlogTextArea                matlab.ui.control.TextArea
        ChooseoperationalmodalityButtonGroup  matlab.ui.container.ButtonGroup
        SimulationButton                matlab.ui.control.RadioButton
        SampletestingButton             matlab.ui.control.RadioButton
        ChoosePTRMSdatafilesfolderButton  matlab.ui.control.Button
        FolderSelectedLamp              matlab.ui.control.Lamp
        MonitoringtargetPanel           matlab.ui.container.Panel
        TargetmassEditFieldLabel        matlab.ui.control.Label
        TargetmassEditField             matlab.ui.control.NumericEditField
        mqLabel_4                       matlab.ui.control.Label
        NeighbourhoodEditFieldLabel     matlab.ui.control.Label
        NeighbourhoodEditField          matlab.ui.control.NumericEditField
        mqLabel_5                       matlab.ui.control.Label
        STARTButton                     matlab.ui.control.Button
        STOPButton                      matlab.ui.control.Button
        ElapsedtimeEditFieldLabel       matlab.ui.control.Label
        ElapsedtimeEditField            matlab.ui.control.EditField
        RightPanel                      matlab.ui.container.Panel
        TabGroup                        matlab.ui.container.TabGroup
        CalibrationTab                  matlab.ui.container.Tab
        CalibrationPanel                matlab.ui.container.Panel
        StartnewcalibrationButton       matlab.ui.control.Button
        CalibrationacquiredLampLabel    matlab.ui.control.Label
        CalibrationacquiredLamp         matlab.ui.control.Lamp
        ChoosecalibrationfileButton     matlab.ui.control.Button
        UIAxesCal                       matlab.ui.control.UIAxes
        CalibrationparametersPanel      matlab.ui.container.Panel
        LOG10INTENSITYLabel             matlab.ui.control.Label
        SlopeEditField                  matlab.ui.control.NumericEditField
        xLOG10CONCLabel                 matlab.ui.control.Label
        InterceptEditField              matlab.ui.control.NumericEditField
        SimulationTab                   matlab.ui.container.Tab
        Panel2_3                        matlab.ui.container.Panel
        MockfiletimespanEditFieldLabel  matlab.ui.control.Label
        MockfiletimespanEditField       matlab.ui.control.NumericEditField
        sLabel                          matlab.ui.control.Label
        SimulatedmassEditFieldLabel     matlab.ui.control.Label
        SimulatedmassEditField          matlab.ui.control.NumericEditField
        mqLabel                         matlab.ui.control.Label
        Panel2_5                        matlab.ui.container.Panel
        MockfilegenerationLampLabel     matlab.ui.control.Label
        MockfilegenerationLamp          matlab.ui.control.Lamp
        UIAxesSim                       matlab.ui.control.UIAxes
        MonitoringTab                   matlab.ui.container.Tab
        UIAxesMon1                      matlab.ui.control.UIAxes
        UIAxesMon2                      matlab.ui.control.UIAxes
        HistoryTab                      matlab.ui.container.Tab
        WholeHistory                    matlab.ui.control.UIAxes
        SlidingHistory                  matlab.ui.control.UIAxes
        TimeSlider                      matlab.ui.control.Slider
        ImportoptionsPanel              matlab.ui.container.Panel
        ChoosefolderfileButton          matlab.ui.control.Button
        ImportexperimentButton          matlab.ui.control.Button
        HistoryDataLamp                 matlab.ui.control.Lamp
        ImportGCMScontrolsButton        matlab.ui.control.Button
        ChoosetoimportDropDownLabel     matlab.ui.control.Label
        ChoosetoimportDropDown          matlab.ui.control.DropDown
        GraphoptionsPanel               matlab.ui.container.Panel
        PlottedtimeintervalminSpinnerLabel  matlab.ui.control.Label
        PlottedtimeintervalminSpinner   matlab.ui.control.Spinner
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    
    properties (Access = private)
        TimerMock % Timer object to time the generation of the mockfiles
        TimerAnalysisMock % Timer object to time the acquisition of data from
                          % the mockfiles
        TimerAnalysis % Timer object to time the acquisition of data from
                      % the PTRMS files
        TimerStopWatch % Timer object to account for elapsed time
        TicTimer % Object to indicate the startpoint of the stopwatch
        LastMock % Name of the last mockfile that has been generated
        ErrorNumber % Indicator of the consecutive error to be logged
        CalSlope % Slope of the current calibration
        CalIntercept % Intercept of the current calibration
        DataFolder % Folder containing the data outputted by the PTRMS
        LastAnalysedFile % Name of the last file that has been loaded
        NoAnalysedFiles % Number of files analysed so far
        StartTimerVec % Vector indicating the starting time of the
                     % current acquisition
        HistoryData % Folder containing all of teh files acquired during
                    % an experiment (can be a pointer to both a folder or a
                    % file)
        EndOfSlider % Last time value of the current monitoring history
                    % that has been imported
        TimeToPlotSlide % Array of time values from the imported monitoring
                        % history
        ConcToPlotSlide % Array of conc values from the imported monitoring
                        % history
    end
    
    methods (Access = private)
        
        
%**************************************************************************
%**************************************************************************
% FUNCTION 1                                                              
%**************************************************************************
%**************************************************************************
        function ClearGraphsNoCalFcn(app,~,~)
            % Stores info on current A/R and position
            aspect_ratio1 = app.UIAxesSim.PlotBoxAspectRatio;
            position1 = app.UIAxesSim.Position;
            % Wipes UIAxesSim
            cla(app.UIAxesSim,'reset');
            % Creates new axes
            app.UIAxesSim = uiaxes(app.SimulationTab);
            title(app.UIAxesSim, '')
            xlabel(app.UIAxesSim, 'Time [s]')
            ylabel(app.UIAxesSim, 'Intensity [ions]')
            app.UIAxesSim.PlotBoxAspectRatio = aspect_ratio1;
            app.UIAxesSim.Position = position1;
            
            
            % Stores info on current A/R and position
            aspect_ratio2 = app.UIAxesMon1.PlotBoxAspectRatio;
            position2 = app.UIAxesMon1.Position;
            % Wipes UIAxesMon1
            cla(app.UIAxesMon1,'reset');
            % Creates new axes
            app.UIAxesMon1 = uiaxes(app.MonitoringTab);
            title(app.UIAxesMon1, 'Raw headspace ion count')
            xlabel(app.UIAxesMon1, 'Time [min]')
            ylabel(app.UIAxesMon1, 'Intensity [ions]')
            app.UIAxesMon1.PlotBoxAspectRatio = aspect_ratio2;
            app.UIAxesMon1.Position = position2; 
            
            
            % Stores info on current A/R and position
            aspect_ratio3 = app.UIAxesMon2.PlotBoxAspectRatio;
            position3 = app.UIAxesMon2.Position;
            % Wipes UIAxesMon2
            cla(app.UIAxesMon2,'reset');
            % Creates new axes
            app.UIAxesMon2 = uiaxes(app.MonitoringTab);
            title(app.UIAxesMon2,['2 min rolling average liquid ' ...
                'concentration'])
            xlabel(app.UIAxesMon2, 'Time [min]')
            ylabel(app.UIAxesMon2, 'Concentration [mg/L]')
            app.UIAxesMon2.PlotBoxAspectRatio = aspect_ratio3;
            app.UIAxesMon2.Position = position3;
        end
    


%**************************************************************************
%**************************************************************************
% FUNCTION 2                                                              
%**************************************************************************
%**************************************************************************
        function GenerateMockFcn(app,~,~)
            user_mass= app.SimulatedmassEditField.Value;
            timelength_mock= app.MockfiletimespanEditField.Value;
            app.LastMock.Value = ...
                genaddh5mock(user_mass,timelength_mock);
        end

        
      
%**************************************************************************
%**************************************************************************
% FUNCTION 3                                                              
%**************************************************************************
%**************************************************************************
        function StartFcn(app,~,~)
            ClearGraphsNoCalFcn(app)

            if app.SimulationButton.Value == true
                start(app.TimerAnalysisMock)
            else
                start(app.TimerAnalysis)
                
            end
        end

        
        
%**************************************************************************
%**************************************************************************
% FUNCTION 4
%**************************************************************************
%**************************************************************************
        function ReadMockFcn(app,~,~)
            app.DataacquisitionLamp.Color = [0 1 0];
            
            filename = app.LastMock.Value;
            central_mass = app.TargetmassEditField.Value + 1;
            neighbourhood = app.NeighbourhoodEditField.Value;
            [cumpeakprof,~,times] = geth5mrcumpeaks( ...
                filename, ...
                central_mass, neighbourhood);
            
            % =============================================================
            % Saves time and central ion profile to plot over time
            % =============================================================
            cumpeakprof = cumpeakprof';
            times = times';
            
            mocks = dir('PTRMSmocksequence*.h5');
            n_mocks = length(mocks);
            
            if n_mocks <= 1
                x = times;
                y = cumpeakprof;
                RawInt = y;
                Time = x;
            else
                load('whole_simulation.mat')
                x = [x; times];
                y = [y; cumpeakprof];
                RawInt = y;
                Time = x;
            end
            
            save('whole_simulation.mat','x','y','Time','RawInt')
            
            
            % =============================================================
            % Plots ions peaks distribution across elapsed time since the
            % first mockfile that has been generated
            % =============================================================
            plot(app.UIAxesSim,times,cumpeakprof)
            plot(app.UIAxesMon1,x,y)
            
            rolling_rng = 360; % Timespan for rolling range graph (~=2mins)
            y = (y - app.CalIntercept.Value)./app.CalSlope.Value;
            
            len = length(x);
            if len > rolling_rng
                X = x(rolling_rng:end);
                Y = zeros((len - (rolling_rng - 1)),1);
                for i = rolling_rng:len
                    Y(i - (rolling_rng - 1)) = ...
                        sum(y((i - (rolling_rng - 1)):i))/rolling_rng;
                end
                
                plot(app.UIAxesMon2,X,Y)
            end
        
            % =============================================================
            % Gets rid of past, useless mockfiles (leaves 5 buffer files)
            % =============================================================
            if n_mocks > 5
                delete(mocks(1).name)
            end
    
        end
    
    
    
%**************************************************************************
%**************************************************************************
% FUNCTION 5
%**************************************************************************
%**************************************************************************
        function ReadFcn(app,~,~)
            app.DataacquisitionLamp.Color = [0 1 0];
            
            cd(app.DataFolder.Value)
            filelist = dir('*.h5');

            numeric_dates = zeros(length(filelist),1);
            for j = 1:length(filelist)
                numeric_dates(j) = datenum(filelist(j).date);
            end
            [~, index]= max(numeric_dates);
            filename = filelist(index).name;
            
            
            check = strcmp(filename(1), '2');
            if ~check
                return
            end
            
            
            if strcmp(filename, app.LastAnalysedFile.Value)
                return
            end

            
            app.NoAnalysedFiles.Value = app.NoAnalysedFiles.Value + 1;
                
            central_mass = app.TargetmassEditField.Value + 1;
            neighbourhood = app.NeighbourhoodEditField.Value;
            [cumpeakprof,~,times] = geth5mrcumpeaks( ...
                filename, ...
                central_mass, neighbourhood);

            % =============================================================
            % Saves time and central ion profile to plot over time
            % =============================================================
            cumpeakprof = cumpeakprof';
            times = times';
            
            FormatInput = 'yyyymmdd_HHMMSS';
            [~,DateStringCurr,~] = fileparts(filename);

            if isempty(app.LastAnalysedFile.Value)
                app.StartTimerVec.Value = ...
                    datevec(DateStringCurr,FormatInput);
                x = times;
                y = cumpeakprof;
                RawInt = y;
                Time = x;
            else
                load('whole_simulation.mat')
                endtime_current_file = datevec(DateStringCurr,FormatInput);

                [~,DateStringPast,~] = ...
                    fileparts(app.LastAnalysedFile.Value);
                endtime_past_file = datevec(DateStringPast,FormatInput);
                

                t_elapsed_from_beg = ...
                    etime(endtime_past_file, app.StartTimerVec.Value) ...
                    + times(end);
                t_between_files = ...
                    etime(endtime_current_file, endtime_past_file);
                addition = ...
                    t_elapsed_from_beg + t_between_files - times(end);
            
                times = times + addition;
                
                x = [x; times];
                y = [y; cumpeakprof];
                RawInt = y;
                Time = x;
                
            end
            
            save('whole_simulation.mat','x','y','Time','RawInt')
            
            
            % =============================================================
            % Plots ions peaks distribution across elapsed time since the
            % first mockfile that has been generated
            % =============================================================
            if length(x) < 1000
                plot(app.UIAxesMon1, x./60, y)
                xlim(app.UIAxesMon1, [0 5])
                xlim(app.UIAxesMon2, [0 5])
            else
                plot(app.UIAxesMon1, x((end-999):end)./60,y((end-999):end))
                xlim(app.UIAxesMon1, [(x(end)./60 - 5) x(end)./60])
                xlim(app.UIAxesMon2, [(x(end)./60 - 5) x(end)./60])
            end
            
            rolling_rng = 400; % Timespan for rolling range graph (~=2mins)
            y = 0.858*0.001*(10.^...
                ((log10(y) - app.CalIntercept.Value)./app.CalSlope.Value));
            
            len = length(x);
            if len > rolling_rng && len <= 1500
                X = x(rolling_rng:end)./60;
                Y = [];
                for i = rolling_rng:len
                    Y = [Y; ...
                        sum(y((i - (rolling_rng-1)):i))/rolling_rng];
                end
                
                plot(app.UIAxesMon2,X,Y)
            elseif len > 1500
                X = x((end-999):end)./60;
                Y = [];
                for i = (len-999):len
                    Y = [Y ...
                        sum(y((i - (rolling_rng-1)):i))/rolling_rng];
                end
                
                plot(app.UIAxesMon2,X,Y)
            end
    
            app.LastAnalysedFile.Value = filename;
        end

        
%**************************************************************************
%**************************************************************************
% FUNCTION 6
%**************************************************************************
%**************************************************************************
        function StopWatchStartFcn(app,~,~)
            TocTimer = 0;
            TocTimer = seconds(TocTimer);
            TocTimer.Format = 'hh:mm:ss';
            text = char(TocTimer);
            app.ElapsedtimeEditField.Value = text;
        end

        
        
%**************************************************************************
%**************************************************************************
% FUNCTION 7
%**************************************************************************
%**************************************************************************        
        function StopWatchFcn(app,~,~)
            TocTimer = toc;
            TocTimer = seconds(TocTimer);
            TocTimer.Format = 'hh:mm:ss';
            text = char(TocTimer);
            app.ElapsedtimeEditField.Value = text;
        end        

        
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            format long
            
            app.ErrorNumber.Value = 0;
            app.LastMock.Value = '';
            app.DataFolder.Value = '';
            app.LastAnalysedFile.Value = '';
            app.NoAnalysedFiles.Value = 0;
            app.StartTimerVec.Value = [];
            app.HistoryData.Value = '';
            app.EndOfSlider.Value = 0;
            app.TimeToPlotSlide.Value = [];
            app.ConcToPlotSlide.Value = [];
            
            app.MockfilegenerationLamp.Color = [.8 .8 .8];
            app.DataacquisitionLamp.Color = [.8 .8 .8];
            app.CalibrationacquiredLamp.Color = [.8 .8 .8];
            app.FolderSelectedLamp.Color = [.8 .8 .8];
            
            app.MockfiletimespanEditField.Value = 30;
            app.SimulatedmassEditField.Value = 154.25;
            app.TargetmassEditField.Value = 154.25;
            app.NeighbourhoodEditField.Value = 0.2;
            
            app.ElapsedtimeEditField.Value = '00:00:00';
                
            
            
            delete 'PTRMSmocksequence*.h5' 'whole_simulation*'
            
            period = 30;    %Period for timer (in seconds)
            
            app.TimerMock = timer(...
                'ExecutionMode', 'fixedRate', ...   
                'Period', period, ...                
                'BusyMode', 'queue');
            app.TimerMock.TimerFcn = @(~,~) app.GenerateMockFcn;
            app.TimerMock.StartFcn = @(~,~) app.StartFcn;
            
            app.TimerAnalysisMock = timer(...
                'StartDelay', 2, ...
                'ExecutionMode', 'fixedRate', ...    
                'Period', period, ...           
                'BusyMode', 'queue');
            app.TimerAnalysisMock.TimerFcn = @(~,~) app.ReadMockFcn;  
            
            app.TimerAnalysis = timer(...
                'ExecutionMode', 'fixedRate', ... 
                'Period', 5, ...        
                'BusyMode', 'queue');
            app.TimerAnalysis.TimerFcn = @(~,~) app.ReadFcn;  
            
            app.TimerStopWatch = timer(...
                'ExecutionMode', 'fixedRate', ...   
                'Period', 0.25, ...              
                'BusyMode', 'queue');
            app.TimerStopWatch.TimerFcn = @(~,~) app.StopWatchFcn;
            app.TimerStopWatch.StartFcn = @(~,~) app.StopWatchStartFcn;
        end

        % Button pushed function: STARTButton
        function STARTButtonPushed(app, event)
            % Checks that all the necessary inputs are non-null
            % Checks on calibration to be acquired
            check = app.CalibrationacquiredLamp.Color == [0 1 0];
            if ~check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Need to import calibration first.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            
            % Checks on target mass
            check = app.TargetmassEditField.Value <= 1;
            if check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Target mass needs to be higher than 1.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            %Checks on neighbourhood
            check = app.NeighbourhoodEditField.Value <= 0;
            if check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Neighbourhood needs to be a positive number.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            %Checks on simulation timespan
            simulation = app.SimulationButton.Value == true;
            check = app.MockfiletimespanEditField.Value <= 5;
            if simulation && check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['The timespan covered by each mockfile ' ...
                    'can''t be lower than 5s.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            %Checks on simulated mass
            simulation = app.SimulationButton.Value == true;
            check = app.SimulatedmassEditField.Value <= 1;
            if simulation && check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Simulated mass needs to be higher than 1.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            %Checks on selected folder for PTRMS datafiles
            if app.SimulationButton.Value == false
                check = isempty(app.DataFolder.Value);
                if check
                    app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                    heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                    text = ['Need to select the folder containing ' ...
                        'PTRMS datafiles when operating in ''Sample ' ...
                        'Testing'' mode.'];
                    app.ErrorlogTextArea.Value = ...
                        [heading ...
                        newline text ...
                        newline ' ' ...
                        newline ' '];
                    return
                end
            end
            
            
            close all
            delete 'PTRMSmocksequence*.h5' 'whole_simulation*'
            
            app.LastAnalysedFile.Value = '';
            app.StartTimerVec.Value = [];
            app.NoAnalysedFiles.Value = 0;
            
            if app.SimulationButton.Value == true
                start(app.TimerMock)
            else
                StartFcn(app)
            end
            
            StopwatchStatus = get(app.TimerStopWatch,'Running');
            if strcmp(StopwatchStatus, 'on')
                stop(app.TimerStopWatch)
            end
            start(app.TimerStopWatch)
            
            app.MockfilegenerationLamp.Color = [0 1 0];
            
            tic
        end

        % Button pushed function: STOPButton
        function STOPButtonPushed(app, event)
            load('whole_simulation.mat', 'x', 'y')
            
            check = isempty(x) | isempty(y);
            if ~check
                time = datetime('now','format','yyyy.MM.dd-HH.mm.ss');
                temp = datestr(time,'yyyy.mm.dd-HH.MM.ss');
                name = ['monitoring-' temp '.mat'];
                filter = {'*.mat'};
                [file,path] = uiputfile(filter,'BioControl 1.0',name);
                
                check = file == 0 || path == 0;
                if check
                    msg = 'Sure about not saving your acquisition?';
                    title = 'BioControl 1.0';
                    selection = uiconfirm(app.BioControl,msg,title,...
                       'Options',{'Don''t save','Cancel'},...
                       'DefaultOption',2,'CancelOption',2, ...
                       'Icon','warning');
                    switch selection
                        case 'Don''t save'
                            stop(app.TimerMock)
                            stop(app.TimerAnalysis)
                            stop(app.TimerStopWatch)
                
                            app.MockfilegenerationLamp.Color = [.8 .8 .8];
                            app.DataacquisitionLamp.Color = [.8 .8 .8];
                        case 'Cancel'
                            return
                    end
                else
                    Time = x;
                    RawInt = y;
                    directory = fullfile(path, file);
                    save(directory,"Time","RawInt")
                    
                    stop(app.TimerMock)
                    stop(app.TimerAnalysis)
                    stop(app.TimerStopWatch)
        
                    app.MockfilegenerationLamp.Color = [.8 .8 .8];
                    app.DataacquisitionLamp.Color = [.8 .8 .8];
                end
            else
                stop(app.TimerMock)
                stop(app.TimerAnalysis)
                stop(app.TimerStopWatch)
    
                app.MockfilegenerationLamp.Color = [.8 .8 .8];
                app.DataacquisitionLamp.Color = [.8 .8 .8];
            end
                    
            
        end

        % Close request function: BioControl
        function BioControlCloseRequest(app, event)
            stop(app.TimerMock)
            stop(app.TimerAnalysis)
            stop(app.TimerStopWatch)
            
            close all
            delete 'PTRMSmocksequence*.h5' 'whole_simulation*'
            
            delete(app)
        end

        % Button pushed function: StartnewcalibrationButton
        function StartnewcalibrationButtonPushed(app, event)
            cal
        end

        % Button pushed function: ChoosecalibrationfileButton
        function ChoosecalibrationfileButtonPushed(app, event)
            filter = '.xlsx';
            title = 'BioControl 1.0';
            [file, path] = uigetfile(filter, title);
            
            check = ischar(file) & ischar(path);
            if ~check
                return
            end
            
            filename = fullfile(path, file);
            
            fit_parameters = readcell(filename,'Range','A6:B6');
            fit_parameters = cell2mat(fit_parameters);
            exp_points = readcell(filename,'Range','A8:B20');
            for i = 1:length(exp_points)
                check = ismissing(exp_points{i});
                if check
                    exp_points = exp_points(1:(i-1),:);
                end
            end
            exp_points = cell2mat(exp_points);
            
            X = [exp_points(1,1) exp_points(end,1)];
            X = log10(X);
            Y = X.*fit_parameters(1) + fit_parameters(2);
            plot(app.UIAxesCal,10.^X,10.^Y,'LineWidth', 2,'Color','b');
            
            hold(app.UIAxesCal, 'on')
            
            scatter(app.UIAxesCal,exp_points(:,1),exp_points(:,2), ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'r', ...
                'SizeData', 50);
            
            set(app.UIAxesCal,'YScale','log')
            set(app.UIAxesCal,'XScale','log')
            grid(app.UIAxesCal, 'on')
            

            app.CalSlope.Value = fit_parameters(1);
            app.CalIntercept.Value = fit_parameters(2);
            
            app.InterceptEditField.Value = app.CalIntercept.Value;
            app.SlopeEditField.Value = app.CalSlope.Value;
            
            app.CalibrationacquiredLamp.Color = [0 1 0];
        end

        % Button pushed function: ChoosePTRMSdatafilesfolderButton
        function ChoosePTRMSdatafilesfolderButtonPushed(app, event)
            title = 'BioControl 1.0';
            selpath = uigetdir(title);
            
            if selpath ~= 0
                app.FolderSelectedLamp.Color = [0 1 0];
                app.DataFolder.Value = selpath;
            end
        end

        % Button pushed function: ChoosefolderfileButton
        function ChoosefolderfileButtonPushed(app, event)
            switch app.ChoosetoimportDropDown.Value
                case 'Raw data'
                    title = 'BioControl 1.0';
                    selpath = uigetdir(title);
                    
                    if selpath ~= 0
                        app.HistoryDataLamp.Color = [0 1 0];
                        app.HistoryData.Value = selpath;
                    end
                case 'Pre-processed data'
                    filter = '.mat';
                    title = 'BioControl 1.0';
                    [file, path] = uigetfile(filter, title);
                    
                    check = ischar(file) & ischar(path);
                    if ~check
                        return
                    end
                    
                    app.HistoryData.Value = fullfile(path, file);
                    app.HistoryDataLamp.Color = [0 1 0];
            end
        end

        % Value changed function: ChoosetoimportDropDown
        function ChoosetoimportDropDownValueChanged(app, event)
            value = app.ChoosetoimportDropDown.Value;
            app.HistoryDataLamp.Color = [.8 .8 .8];
        end

        % Button pushed function: ImportexperimentButton
        function ImportexperimentButtonPushed(app, event)
            % Checks that all the necessary inputs are non-null
            % Checks on folder/file pointer to be present
            check = app.HistoryDataLamp.Color == [0 1 0];
            if ~check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Need to import folder/file first.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            % Checks on calibration to be acquired
            check = app.CalibrationacquiredLamp.Color == [0 1 0];
            if ~check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Need to import calibration first.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            
            % Checks on target mass
            check = app.TargetmassEditField.Value <= 1;
            if check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Target mass needs to be higher than 1.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            %Checks on neighbourhood
            check = app.NeighbourhoodEditField.Value <= 0;
            if check
                app.ErrorNumber.Value = app.ErrorNumber.Value + 1;
                heading = ['ERROR NO.' num2str(app.ErrorNumber.Value)];
                text = ['Neighbourhood needs to be a positive number.'];
                app.ErrorlogTextArea.Value = ...
                    [heading ...
                    newline text ...
                    newline ' ' ...
                    newline ' '];
                return
            end
            
            switch app.ChoosetoimportDropDown.Value
                case 'Raw data'
                    cd(app.HistoryData.Value)
                    filelist = dir('*.h5');
                    
                    Dialog = uiprogressdlg(app.BioControl, ...
                        'Title', 'BioControl 1.0', ...
                        'Message', 'Loading files...');
                    
                    files_count = length(filelist);
                    
                    for i = 1:files_count
                        Dialog.Value = i/files_count;
                        curr_filename = filelist(i).name;
                        
                        central_mass = app.TargetmassEditField.Value + 1;
                        neighbourhood = app.NeighbourhoodEditField.Value;
                        [cumpeakprof,~,times] = geth5mrcumpeaks( ...
                            curr_filename, ...
                            central_mass, neighbourhood);


                        cumpeakprof = cumpeakprof';
                        times = times';
            
                        FormatInput = 'yyyymmdd_HHMMSS';
                        [~,DateStringCurr,~] = fileparts(curr_filename);
            
                        if i == 1
                            app.StartTimerVec.Value = ...
                                datevec(DateStringCurr,FormatInput);
                            x = times;
                            y = cumpeakprof;
                        else
                            endtime_current_file = ...
                                datevec(DateStringCurr,FormatInput);

                            [~,DateStringPast,~] = ...
                                fileparts(past_filename);
                            endtime_past_file = ...
                                datevec(DateStringPast,FormatInput);
                

                            t_elapsed_from_beg = ...
                                etime(...
                                endtime_past_file, ...
                                app.StartTimerVec.Value) ...
                                + times(end);
                            t_between_files = ...
                                etime(...
                                endtime_current_file, endtime_past_file);
                            addition = ...
                                t_elapsed_from_beg + t_between_files ...
                                - times(end);
                        
                            times = times + addition;
                            
                            x = [x; times];
                            y = [y; cumpeakprof];
                        end
            
                        past_filename = curr_filename;
                    end
        
                    Time = x;
                    RawInt = y;
                    save('monitoring_history.mat',"Time","RawInt")
        
                    app.TimeSlider.Enable = 'on';
                    
                    app.EndOfSlider.Value = x(end)/60; % in minutes
                    
                    X = x(400:end);
                    app.TimeToPlotSlide.Value = X;
                
                    y = 0.858*0.001*...
                        (10.^...
                        ((log10(y)-app.CalIntercept.Value)./...
                        app.CalSlope.Value));
                    
                    Y = [];
                    for i = 400:length(x)
                    Y = [Y; ...
                        sum(y((i - (399)):i))/400];
                    end
                    app.ConcToPlotSlide.Value = Y;
                    
                    plot(app.WholeHistory, ...
                        X(1:300:end)./3600, Y(1:300:end))
                    xlim(app.WholeHistory, [x(1)/3600 x(end)/3600])
                
                    plot(app.SlidingHistory, ...
                        X(1:(app.PlottedtimeintervalminSpinner.Value*200))./60,...
                        Y(1:(app.PlottedtimeintervalminSpinner.Value*200)))
                    xlim(app.SlidingHistory, ...
                        [0 app.PlottedtimeintervalminSpinner.Value])
                    
                    close(Dialog)
                    
                case 'Pre-processed data'
                    load(app.HistoryData.Value)
                    x = Time;
                    y = RawInt;
                    
                    app.TimeSlider.Enable = 'on';
                    
                    app.EndOfSlider.Value = x(end)/60; % in minutes

                    X = x(400:end);
                    app.TimeToPlotSlide.Value = X;
                
                    y = 0.858*0.001*...
                        (10.^...
                        ((log10(y)-app.CalIntercept.Value)./...
                        app.CalSlope.Value));
                    
                    Y = [];
                    for i = 400:length(x)
                    Y = [Y; ...
                        sum(y((i - (399)):i))/400];
                    end
                    app.ConcToPlotSlide.Value = Y;
                    
                    plot(app.WholeHistory, ...
                        X(1:300:end)./3600, Y(1:300:end))
                    xlim(app.WholeHistory, [x(1)/3600 x(end)/3600])
                
                    plot(app.SlidingHistory, ...
                        X(1:(app.PlottedtimeintervalminSpinner.Value*200))./60,...
                        Y(1:(app.PlottedtimeintervalminSpinner.Value*200)))
                    xlim(app.SlidingHistory, ...
                        [0 app.PlottedtimeintervalminSpinner.Value])
            end
        end

        % Value changing function: TimeSlider
        function TimeSliderValueChanging(app, event)
            changingValue = event.Value;
            
            tot_range = ...
                app.EndOfSlider.Value ...
                - app.PlottedtimeintervalminSpinner.Value;
            plot_range = [(changingValue*tot_range), ...
                (changingValue*tot_range ...
                + app.PlottedtimeintervalminSpinner.Value)];

            indexes = ...
                find((app.TimeToPlotSlide.Value./60) > plot_range(1) & ...
                (app.TimeToPlotSlide.Value./60) < plot_range(2));
            
            plot(app.SlidingHistory, ...
                app.TimeToPlotSlide.Value(indexes)./60, ...
                app.ConcToPlotSlide.Value(indexes))
            xlim(app.SlidingHistory, plot_range)
            
            ylim(app.SlidingHistory, ...
                [min(app.ConcToPlotSlide.Value(indexes)), ...
                max(app.ConcToPlotSlide.Value(indexes))])
            
        end

        % Button pushed function: ImportGCMScontrolsButton
        function ImportGCMScontrolsButtonPushed(app, event)
            condition = true;
            while condition
                prompt = {'Enter time value [min]:', ...
                    'Enter concentration value [mg/L];'};
                dlgtitle = 'BioControl 1.0';
                dims = [1 50];
                answer = inputdlg(prompt, dlgtitle, dims);
                
                check = isempty(answer);
                if check
                    return
                end
                
                x = str2double(answer{1});
                y = str2double(answer{2});
                
                check = isempty(answer{1}) | isempty(answer{2}) | ...
                    isnan(x) | isnan(y);
                if check
                    uialert(app.BioControl, ...
                        'Entered inputs are not valid!', ...
                        'BioControl 1.0', 'Icon', 'error');
                    condition = false;
                else
                    hold(app.WholeHistory, 'on')
                    scatter(app.WholeHistory,x./60,y, ...
                        'Marker', 'o', ...
                        'MarkerFaceColor', 'r', ...
                        'SizeData', 50);
                    hold(app.WholeHistory, 'off')
                end
            end
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.BioControl.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {655, 655};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {236, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create BioControl and hide until all components are created
            app.BioControl = uifigure('Visible', 'off');
            app.BioControl.AutoResizeChildren = 'off';
            app.BioControl.Position = [100 100 1076 655];
            app.BioControl.Name = 'BioControl 1.0';
            app.BioControl.CloseRequestFcn = createCallbackFcn(app, @BioControlCloseRequest, true);
            app.BioControl.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create FileMenu
            app.FileMenu = uimenu(app.BioControl);
            app.FileMenu.Text = 'File';

            % Create OpenfolderMenu
            app.OpenfolderMenu = uimenu(app.FileMenu);
            app.OpenfolderMenu.Text = 'Open folder';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.BioControl);
            app.GridLayout.ColumnWidth = {236, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create StatusPanel
            app.StatusPanel = uipanel(app.LeftPanel);
            app.StatusPanel.Title = 'Status';
            app.StatusPanel.Position = [25 15 189 294];

            % Create DataacquisitionLampLabel
            app.DataacquisitionLampLabel = uilabel(app.StatusPanel);
            app.DataacquisitionLampLabel.HorizontalAlignment = 'right';
            app.DataacquisitionLampLabel.Position = [12 237 92 22];
            app.DataacquisitionLampLabel.Text = 'Data acquisition';

            % Create DataacquisitionLamp
            app.DataacquisitionLamp = uilamp(app.StatusPanel);
            app.DataacquisitionLamp.Position = [144 238 20 20];
            app.DataacquisitionLamp.Color = [0.8 0.8 0.8];

            % Create ErrorlogTextAreaLabel
            app.ErrorlogTextAreaLabel = uilabel(app.StatusPanel);
            app.ErrorlogTextAreaLabel.HorizontalAlignment = 'right';
            app.ErrorlogTextAreaLabel.Position = [14 168 55 22];
            app.ErrorlogTextAreaLabel.Text = 'Error log:';

            % Create ErrorlogTextArea
            app.ErrorlogTextArea = uitextarea(app.StatusPanel);
            app.ErrorlogTextArea.Position = [14 10 162 156];

            % Create ChooseoperationalmodalityButtonGroup
            app.ChooseoperationalmodalityButtonGroup = uibuttongroup(app.LeftPanel);
            app.ChooseoperationalmodalityButtonGroup.Title = 'Choose operational modality';
            app.ChooseoperationalmodalityButtonGroup.Position = [24 496 189 143];

            % Create SimulationButton
            app.SimulationButton = uiradiobutton(app.ChooseoperationalmodalityButtonGroup);
            app.SimulationButton.Text = 'Simulation';
            app.SimulationButton.Position = [12 90 78 22];
            app.SimulationButton.Value = true;

            % Create SampletestingButton
            app.SampletestingButton = uiradiobutton(app.ChooseoperationalmodalityButtonGroup);
            app.SampletestingButton.Text = 'Sample testing';
            app.SampletestingButton.Position = [12 60 102 22];

            % Create ChoosePTRMSdatafilesfolderButton
            app.ChoosePTRMSdatafilesfolderButton = uibutton(app.ChooseoperationalmodalityButtonGroup, 'push');
            app.ChoosePTRMSdatafilesfolderButton.ButtonPushedFcn = createCallbackFcn(app, @ChoosePTRMSdatafilesfolderButtonPushed, true);
            app.ChoosePTRMSdatafilesfolderButton.Position = [12 15 102 36];
            app.ChoosePTRMSdatafilesfolderButton.Text = {'Choose PTRMS'; 'datafiles folder'};

            % Create FolderSelectedLamp
            app.FolderSelectedLamp = uilamp(app.ChooseoperationalmodalityButtonGroup);
            app.FolderSelectedLamp.Position = [145 23 20 20];
            app.FolderSelectedLamp.Color = [0.8 0.8 0.8];

            % Create MonitoringtargetPanel
            app.MonitoringtargetPanel = uipanel(app.LeftPanel);
            app.MonitoringtargetPanel.Title = 'Monitoring target';
            app.MonitoringtargetPanel.Position = [25 387 189 100];

            % Create TargetmassEditFieldLabel
            app.TargetmassEditFieldLabel = uilabel(app.MonitoringtargetPanel);
            app.TargetmassEditFieldLabel.HorizontalAlignment = 'right';
            app.TargetmassEditFieldLabel.Position = [6 49 70 22];
            app.TargetmassEditFieldLabel.Text = 'Target mass';

            % Create TargetmassEditField
            app.TargetmassEditField = uieditfield(app.MonitoringtargetPanel, 'numeric');
            app.TargetmassEditField.Position = [102 49 46 22];

            % Create mqLabel_4
            app.mqLabel_4 = uilabel(app.MonitoringtargetPanel);
            app.mqLabel_4.Position = [152 49 27 22];
            app.mqLabel_4.Text = 'm/q';

            % Create NeighbourhoodEditFieldLabel
            app.NeighbourhoodEditFieldLabel = uilabel(app.MonitoringtargetPanel);
            app.NeighbourhoodEditFieldLabel.HorizontalAlignment = 'right';
            app.NeighbourhoodEditFieldLabel.Position = [6 16 89 22];
            app.NeighbourhoodEditFieldLabel.Text = 'Neighbourhood';

            % Create NeighbourhoodEditField
            app.NeighbourhoodEditField = uieditfield(app.MonitoringtargetPanel, 'numeric');
            app.NeighbourhoodEditField.Position = [102 16 46 22];

            % Create mqLabel_5
            app.mqLabel_5 = uilabel(app.MonitoringtargetPanel);
            app.mqLabel_5.Position = [152 16 27 22];
            app.mqLabel_5.Text = 'm/q';

            % Create STARTButton
            app.STARTButton = uibutton(app.LeftPanel, 'push');
            app.STARTButton.ButtonPushedFcn = createCallbackFcn(app, @STARTButtonPushed, true);
            app.STARTButton.FontWeight = 'bold';
            app.STARTButton.FontColor = [0.4667 0.6745 0.1882];
            app.STARTButton.Position = [25 336 88 28];
            app.STARTButton.Text = 'START';

            % Create STOPButton
            app.STOPButton = uibutton(app.LeftPanel, 'push');
            app.STOPButton.ButtonPushedFcn = createCallbackFcn(app, @STOPButtonPushed, true);
            app.STOPButton.FontWeight = 'bold';
            app.STOPButton.FontColor = [1 0 0];
            app.STOPButton.Position = [128 336 86 28];
            app.STOPButton.Text = 'STOP';

            % Create ElapsedtimeEditFieldLabel
            app.ElapsedtimeEditFieldLabel = uilabel(app.LeftPanel);
            app.ElapsedtimeEditFieldLabel.HorizontalAlignment = 'right';
            app.ElapsedtimeEditFieldLabel.Position = [39 219 75 22];
            app.ElapsedtimeEditFieldLabel.Text = 'Elapsed time';

            % Create ElapsedtimeEditField
            app.ElapsedtimeEditField = uieditfield(app.LeftPanel, 'text');
            app.ElapsedtimeEditField.Position = [141 219 60 22];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.RightPanel);
            app.TabGroup.Position = [7 6 828 643];

            % Create CalibrationTab
            app.CalibrationTab = uitab(app.TabGroup);
            app.CalibrationTab.Title = 'Calibration';

            % Create CalibrationPanel
            app.CalibrationPanel = uipanel(app.CalibrationTab);
            app.CalibrationPanel.Tooltip = {'Provide a preliminary calibration of the target compouns(s). Either call the calibration tool or select the appropriate calibration file'};
            app.CalibrationPanel.Title = 'Calibration';
            app.CalibrationPanel.Position = [576 295 222 196];

            % Create StartnewcalibrationButton
            app.StartnewcalibrationButton = uibutton(app.CalibrationPanel, 'push');
            app.StartnewcalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @StartnewcalibrationButtonPushed, true);
            app.StartnewcalibrationButton.Position = [52 118 118 42];
            app.StartnewcalibrationButton.Text = {'Start'; 'new calibration'};

            % Create CalibrationacquiredLampLabel
            app.CalibrationacquiredLampLabel = uilabel(app.CalibrationPanel);
            app.CalibrationacquiredLampLabel.Position = [21 13 64 27];
            app.CalibrationacquiredLampLabel.Text = {'Calibration'; 'acquired'};

            % Create CalibrationacquiredLamp
            app.CalibrationacquiredLamp = uilamp(app.CalibrationPanel);
            app.CalibrationacquiredLamp.Position = [179 15 23 23];
            app.CalibrationacquiredLamp.Color = [0.8 0.8 0.8];

            % Create ChoosecalibrationfileButton
            app.ChoosecalibrationfileButton = uibutton(app.CalibrationPanel, 'push');
            app.ChoosecalibrationfileButton.ButtonPushedFcn = createCallbackFcn(app, @ChoosecalibrationfileButtonPushed, true);
            app.ChoosecalibrationfileButton.Position = [52 62 118 39];
            app.ChoosecalibrationfileButton.Text = {'Choose'; 'calibration file'};

            % Create UIAxesCal
            app.UIAxesCal = uiaxes(app.CalibrationTab);
            title(app.UIAxesCal, 'Calibration graph')
            xlabel(app.UIAxesCal, 'Target compound concentration [ppb vol]')
            ylabel(app.UIAxesCal, 'Intensity [ions/s]')
            app.UIAxesCal.PlotBoxAspectRatio = [1 1.04312114989733 1];
            app.UIAxesCal.Position = [15 16 536 564];

            % Create CalibrationparametersPanel
            app.CalibrationparametersPanel = uipanel(app.CalibrationTab);
            app.CalibrationparametersPanel.Title = 'Calibration parameters';
            app.CalibrationparametersPanel.Position = [567 177 240 100];

            % Create LOG10INTENSITYLabel
            app.LOG10INTENSITYLabel = uilabel(app.CalibrationparametersPanel);
            app.LOG10INTENSITYLabel.HorizontalAlignment = 'right';
            app.LOG10INTENSITYLabel.Position = [2 45 124 22];
            app.LOG10INTENSITYLabel.Text = 'LOG10 (INTENSITY) =';

            % Create SlopeEditField
            app.SlopeEditField = uieditfield(app.CalibrationparametersPanel, 'numeric');
            app.SlopeEditField.HorizontalAlignment = 'center';
            app.SlopeEditField.Position = [6 21 59 22];

            % Create xLOG10CONCLabel
            app.xLOG10CONCLabel = uilabel(app.CalibrationparametersPanel);
            app.xLOG10CONCLabel.HorizontalAlignment = 'right';
            app.xLOG10CONCLabel.Position = [63 22 112 22];
            app.xLOG10CONCLabel.Text = 'x LOG10 (CONC) + ';

            % Create InterceptEditField
            app.InterceptEditField = uieditfield(app.CalibrationparametersPanel, 'numeric');
            app.InterceptEditField.HorizontalAlignment = 'center';
            app.InterceptEditField.Position = [175 21 60 22];

            % Create SimulationTab
            app.SimulationTab = uitab(app.TabGroup);
            app.SimulationTab.Title = 'Simulation';

            % Create Panel2_3
            app.Panel2_3 = uipanel(app.SimulationTab);
            app.Panel2_3.AutoResizeChildren = 'off';
            app.Panel2_3.Title = 'Simulation parameters';
            app.Panel2_3.Position = [593 341 200 94];

            % Create MockfiletimespanEditFieldLabel
            app.MockfiletimespanEditFieldLabel = uilabel(app.Panel2_3);
            app.MockfiletimespanEditFieldLabel.HorizontalAlignment = 'right';
            app.MockfiletimespanEditFieldLabel.Position = [2 41 103 22];
            app.MockfiletimespanEditFieldLabel.Text = 'Mockfile timespan';

            % Create MockfiletimespanEditField
            app.MockfiletimespanEditField = uieditfield(app.Panel2_3, 'numeric');
            app.MockfiletimespanEditField.Position = [117 41 46 22];

            % Create sLabel
            app.sLabel = uilabel(app.Panel2_3);
            app.sLabel.Position = [167 41 25 22];
            app.sLabel.Text = 's';

            % Create SimulatedmassEditFieldLabel
            app.SimulatedmassEditFieldLabel = uilabel(app.Panel2_3);
            app.SimulatedmassEditFieldLabel.HorizontalAlignment = 'right';
            app.SimulatedmassEditFieldLabel.Position = [1 12 91 22];
            app.SimulatedmassEditFieldLabel.Text = 'Simulated mass';

            % Create SimulatedmassEditField
            app.SimulatedmassEditField = uieditfield(app.Panel2_3, 'numeric');
            app.SimulatedmassEditField.Position = [118 12 46 22];

            % Create mqLabel
            app.mqLabel = uilabel(app.Panel2_3);
            app.mqLabel.Position = [168 12 27 22];
            app.mqLabel.Text = 'm/q';

            % Create Panel2_5
            app.Panel2_5 = uipanel(app.SimulationTab);
            app.Panel2_5.AutoResizeChildren = 'off';
            app.Panel2_5.Title = 'Simulation state';
            app.Panel2_5.Position = [592 235 200 71];

            % Create MockfilegenerationLampLabel
            app.MockfilegenerationLampLabel = uilabel(app.Panel2_5);
            app.MockfilegenerationLampLabel.HorizontalAlignment = 'right';
            app.MockfilegenerationLampLabel.Position = [2 17 111 22];
            app.MockfilegenerationLampLabel.Text = 'Mockfile generation';

            % Create MockfilegenerationLamp
            app.MockfilegenerationLamp = uilamp(app.Panel2_5);
            app.MockfilegenerationLamp.Position = [169 18 20 20];
            app.MockfilegenerationLamp.Color = [0.8 0.8 0.8];

            % Create UIAxesSim
            app.UIAxesSim = uiaxes(app.SimulationTab);
            title(app.UIAxesSim, 'Last randomly-generated datachunk')
            xlabel(app.UIAxesSim, 'Time [s]')
            ylabel(app.UIAxesSim, 'Intensity [ions]')
            app.UIAxesSim.PlotBoxAspectRatio = [1 1.04857444561774 1];
            app.UIAxesSim.XTick = [0 0.2 0.4 0.6 0.8 1];
            app.UIAxesSim.Position = [41 28 524 552];

            % Create MonitoringTab
            app.MonitoringTab = uitab(app.TabGroup);
            app.MonitoringTab.Title = 'Monitoring';

            % Create UIAxesMon1
            app.UIAxesMon1 = uiaxes(app.MonitoringTab);
            title(app.UIAxesMon1, 'Raw headspace ion count')
            xlabel(app.UIAxesMon1, 'Time [min]')
            ylabel(app.UIAxesMon1, 'Intensity [ions]')
            app.UIAxesMon1.PlotBoxAspectRatio = [3.03846153846154 1 1];
            app.UIAxesMon1.Position = [29 313 760 290];

            % Create UIAxesMon2
            app.UIAxesMon2 = uiaxes(app.MonitoringTab);
            title(app.UIAxesMon2, '2 min rolling average liquid concentration')
            xlabel(app.UIAxesMon2, 'Time [min]')
            ylabel(app.UIAxesMon2, 'Concentration [mg/L]')
            app.UIAxesMon2.PlotBoxAspectRatio = [3.03846153846154 1 1];
            app.UIAxesMon2.Position = [29 19 760 290];

            % Create HistoryTab
            app.HistoryTab = uitab(app.TabGroup);
            app.HistoryTab.Title = 'History';

            % Create WholeHistory
            app.WholeHistory = uiaxes(app.HistoryTab);
            title(app.WholeHistory, 'Entire acquisition')
            xlabel(app.WholeHistory, 'Time [h]')
            ylabel(app.WholeHistory, 'Concentration [mg/L]')
            app.WholeHistory.PlotBoxAspectRatio = [2.69444444444444 1 1];
            app.WholeHistory.Position = [39 330 590 272];

            % Create SlidingHistory
            app.SlidingHistory = uiaxes(app.HistoryTab);
            title(app.SlidingHistory, 'Zoomed graph')
            xlabel(app.SlidingHistory, 'Time [min]')
            ylabel(app.SlidingHistory, 'Concentration [mg/L]')
            app.SlidingHistory.PlotBoxAspectRatio = [2.59821428571429 1 1];
            app.SlidingHistory.Position = [39 42 590 280];

            % Create TimeSlider
            app.TimeSlider = uislider(app.HistoryTab);
            app.TimeSlider.Limits = [0 1];
            app.TimeSlider.MajorTicks = [];
            app.TimeSlider.ValueChangingFcn = createCallbackFcn(app, @TimeSliderValueChanging, true);
            app.TimeSlider.MinorTicks = [];
            app.TimeSlider.Enable = 'off';
            app.TimeSlider.Position = [126 19 457 3];

            % Create ImportoptionsPanel
            app.ImportoptionsPanel = uipanel(app.HistoryTab);
            app.ImportoptionsPanel.Title = 'Import options';
            app.ImportoptionsPanel.Position = [633 374 183 198];

            % Create ChoosefolderfileButton
            app.ChoosefolderfileButton = uibutton(app.ImportoptionsPanel, 'push');
            app.ChoosefolderfileButton.ButtonPushedFcn = createCallbackFcn(app, @ChoosefolderfileButtonPushed, true);
            app.ChoosefolderfileButton.Position = [14 81 109 36];
            app.ChoosefolderfileButton.Text = {'Choose'; 'folder/file'};

            % Create ImportexperimentButton
            app.ImportexperimentButton = uibutton(app.ImportoptionsPanel, 'push');
            app.ImportexperimentButton.ButtonPushedFcn = createCallbackFcn(app, @ImportexperimentButtonPushed, true);
            app.ImportexperimentButton.Position = [14 46 152 23];
            app.ImportexperimentButton.Text = 'Import experiment';

            % Create HistoryDataLamp
            app.HistoryDataLamp = uilamp(app.ImportoptionsPanel);
            app.HistoryDataLamp.Position = [140 89 20 20];
            app.HistoryDataLamp.Color = [0.8 0.8 0.8];

            % Create ImportGCMScontrolsButton
            app.ImportGCMScontrolsButton = uibutton(app.ImportoptionsPanel, 'push');
            app.ImportGCMScontrolsButton.ButtonPushedFcn = createCallbackFcn(app, @ImportGCMScontrolsButtonPushed, true);
            app.ImportGCMScontrolsButton.Position = [14 13 152 22];
            app.ImportGCMScontrolsButton.Text = 'Import GCMS controls';

            % Create ChoosetoimportDropDownLabel
            app.ChoosetoimportDropDownLabel = uilabel(app.ImportoptionsPanel);
            app.ChoosetoimportDropDownLabel.Position = [14 149 102 23];
            app.ChoosetoimportDropDownLabel.Text = 'Choose to import:';

            % Create ChoosetoimportDropDown
            app.ChoosetoimportDropDown = uidropdown(app.ImportoptionsPanel);
            app.ChoosetoimportDropDown.Items = {'Raw data', 'Pre-processed data'};
            app.ChoosetoimportDropDown.ValueChangedFcn = createCallbackFcn(app, @ChoosetoimportDropDownValueChanged, true);
            app.ChoosetoimportDropDown.Position = [13 129 153 22];
            app.ChoosetoimportDropDown.Value = 'Raw data';

            % Create GraphoptionsPanel
            app.GraphoptionsPanel = uipanel(app.HistoryTab);
            app.GraphoptionsPanel.Title = 'Graph options';
            app.GraphoptionsPanel.Position = [633 144 183 85];

            % Create PlottedtimeintervalminSpinnerLabel
            app.PlottedtimeintervalminSpinnerLabel = uilabel(app.GraphoptionsPanel);
            app.PlottedtimeintervalminSpinnerLabel.Position = [10 36 142 23];
            app.PlottedtimeintervalminSpinnerLabel.Text = 'Plotted time interval [min]';

            % Create PlottedtimeintervalminSpinner
            app.PlottedtimeintervalminSpinner = uispinner(app.GraphoptionsPanel);
            app.PlottedtimeintervalminSpinner.Position = [10 12 140 22];
            app.PlottedtimeintervalminSpinner.Value = 100;

            % Show the figure after all components are created
            app.BioControl.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = platform_test_test

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.BioControl)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.BioControl)
        end
    end
end