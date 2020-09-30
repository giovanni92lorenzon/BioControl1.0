classdef cal < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        BioControl10CalibrationtabUIFigure  matlab.ui.Figure
        UIAxes                         matlab.ui.control.UIAxes
        PaircalibrationfileButton      matlab.ui.control.Button
        CompoundconcentrationppbLabel  matlab.ui.control.Label
        ComponentConcentrationEditField  matlab.ui.control.NumericEditField
        CleargraphButton               matlab.ui.control.Button
        SlopeEditFieldLabel            matlab.ui.control.Label
        SlopeEditField                 matlab.ui.control.NumericEditField
        LinearfitButton                matlab.ui.control.Button
        StorecalibrationButton         matlab.ui.control.Button
        InterceptEditFieldLabel        matlab.ui.control.Label
        InterceptEditField             matlab.ui.control.NumericEditField
        R2EditFieldLabel               matlab.ui.control.Label
        R2EditField                    matlab.ui.control.NumericEditField
        GeneratemockcalibrationButton  matlab.ui.control.Button
        CompoundcentralmassDaLabel     matlab.ui.control.Label
        ComponentCentralMassEditField  matlab.ui.control.NumericEditField
        CompoundmassneighbourhoodDaLabel  matlab.ui.control.Label
        ComponentMassNeighbourhoodEditField  matlab.ui.control.NumericEditField
    end


    properties (Access = private)
        
        TempFolderPath % Path of the temporary folder where mock
                       % calibrations points are stored
        NoPointsOnGraph % Number of points already plotted on the graph
        DataToFit % Array containing concentration&intensity values for
                  % each point added to the calibration
        FitParameters % Structure storing the last fitted parameters
    end
    
    methods (Access = private)
        
        function DeleteH5Content(app,~,~)
            filePattern = fullfile( ...
                app.TempFolderPath.Value, 'PTRMSmockcalsequence*.h5');
            files = dir(filePattern);
            for k = 1:length(files)
              baseFileName = files(k).name;
              fullFileName = fullfile( ...
                  app.TempFolderPath.Value, baseFileName);
              delete(fullFileName);
            end
        end
        
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            format long
            hold(app.UIAxes, 'on')
            
            app.TempFolderPath.Value = '';
            app.NoPointsOnGraph.Value = 0;
            app.DataToFit.Value = [];
            app.FitParameters.Value = [];
        end

        % Button pushed function: PaircalibrationfileButton
        function PaircalibrationfileButtonPushed(app, event)
            check = isempty(app.ComponentConcentrationEditField.Value);
            if check
                uialert(app.BioControl10CalibrationtabUIFigure, ...
                ['Insert a concentration value before proceding ' ...
                'to pair up file and concentration.'], ...
                'BioControl 1.0', 'Icon', 'error');
            
                return
            end
            
            check = app.ComponentCentralMassEditField.Value == 0;
            if check
                uialert(app.BioControl10CalibrationtabUIFigure, ...
                ['Insert target mass value before proceding ' ...
                'to pair up file and concentration.'], ...
                'BioControl 1.0', 'Icon', 'error');
            
                return
            end
            
            check = app.ComponentMassNeighbourhoodEditField.Value == 0;
            if check
                uialert(app.BioControl10CalibrationtabUIFigure, ...
                ['Insert a mass neighbourhood value before proceding ' ...
                'to pair up file and concentration.'], ...
                'BioControl 1.0', 'Icon', 'error');
            
                return
            end
            
            filter = '.h5';
            title = 'BioControl 1.0 - Select one or more calibration files';
            [filename,path] = uigetfile(filter, title,'MultiSelect','on');
            
            multiple_entrance_check = iscell(filename);
            if ~multiple_entrance_check
                single_entrance_check = ischar(filename);
                if ~single_entrance_check
                    return
                end
            end
            
            if multiple_entrance_check
                ext_check = false;
                for i = 1:length(filename)
                    [~,~,extension] = fileparts(filename{i});
                    ext_check = ext_check | ~strcmpi(extension, filter);
                    if ext_check
                        uialert(app.BioControl10CalibrationtabUIFigure, ...
                        ['The selected file doesn''t have the ' ...
                        'expected extension (<.h5>).'], ...
                        'BioControl 1.0', 'Icon', 'error');
                    
                        return
                    end
                end
            else
                [~,~,extension] = fileparts(filename);
                ext_check = ~strcmpi(extension, filter);
                if ext_check
                    uialert(app.BioControl10CalibrationtabUIFigure, ...
                    ['The selected file doesn''t have the ' ...
                    'expected extension (<.h5>).'], ...
                    'BioControl 1.0', 'Icon', 'error');
                
                    return
                end
            end
            
            
            cd(path);
            if multiple_entrance_check
                pruned_data = [];
                for i = 1:length(filename)
                    [cumpeakprof,~,~] = geth5mrcumpeaks( ...
                        filename{i}, ...
                        (app.ComponentCentralMassEditField.Value + 1), ...
                        app.ComponentMassNeighbourhoodEditField.Value);
                    
                    pruned_data = [pruned_data rmoutliers(cumpeakprof)];
                end
            else
                [cumpeakprof,~,~] = geth5mrcumpeaks( ...
                    filename, ...
                    (app.ComponentCentralMassEditField.Value + 1), ...
                    app.ComponentMassNeighbourhoodEditField.Value);
                    
                pruned_data = rmoutliers(cumpeakprof);
            end

            avg_output = mean(pruned_data);
            error = std(pruned_data);
            
            
            errorbar(app.UIAxes, ...
                app.ComponentConcentrationEditField.Value, ...
                avg_output, ...
                error, ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'b');
            
            set(app.UIAxes,'YScale','log')
            set(app.UIAxes,'XScale','log')
            grid(app.UIAxes, 'on')
            
            app.NoPointsOnGraph.Value = app.NoPointsOnGraph.Value + 1;
            
            new_point = ...
                [app.ComponentConcentrationEditField.Value; ...
                avg_output; ...
                error];
            app.DataToFit.Value = [app.DataToFit.Value new_point];
        end

        % Button pushed function: CleargraphButton
        function CleargraphButtonPushed(app, event)
            aspect_ratio = app.UIAxes.PlotBoxAspectRatio;
            position = app.UIAxes.Position;
            
            cla(app.UIAxes,'reset');
            
            % Create UIAxes
            app.UIAxes = uiaxes(app.BioControl10CalibrationtabUIFigure);
            title(app.UIAxes, 'Calibration data')
            xlabel(app.UIAxes, 'Component concentration [ppb vol]')
            ylabel(app.UIAxes, 'Average intensity [ions/s] ')
            app.UIAxes.PlotBoxAspectRatio = aspect_ratio;
            app.UIAxes.Position = position;
            
            hold(app.UIAxes, 'on')
            
            app.NoPointsOnGraph.Value = 0;
            app.DataToFit.Value = [];
            app.FitParameters.Value = [];
        end

        % Button pushed function: LinearfitButton
        function LinearfitButtonPushed(app, event)
            check = app.NoPointsOnGraph.Value >= 2;
            
            if ~check
                uialert(app.BioControl10CalibrationtabUIFigure, ...
                ['Can''t fit less than 2 experimental points. ' ...
                'Add more points to the graph.'], ...
                'BioControl 1.0', 'Icon', 'error');
            
                return
            end
            
                
            X = app.DataToFit.Value(1,:);
            Y = app.DataToFit.Value(2,:);
            logx = log10(X);
            logy = log10(Y);
            
            
            [app.FitParameters.Value,S] = polyfit(logx,logy,1);
            %%%
            R2 = 1 - (S.normr^2)/(norm(logy-mean(logy))^2);
            app.R2EditField.Value = R2;
            %%%
            app.SlopeEditField.Value = app.FitParameters.Value(1);
            app.InterceptEditField.Value = app.FitParameters.Value(2);
            
            A = [min(logx) max(logx)];
            B = polyval(app.FitParameters.Value,A);
            loglog(app.UIAxes,10.^A,10.^B,'Color','r','LineWidth',1)
        end

        % Button pushed function: StorecalibrationButton
        function StorecalibrationButtonPushed(app, event)
            check = isempty(app.DataToFit.Value);
            
            if check
                uialert(app.BioControl10CalibrationtabUIFigure, ...
                ['No data points available for storing. ' ...
                'Plot data points before proceeding.'], ...
                'BioControl 1.0', 'Icon', 'error');
            
                return
            end

            time = datetime('now','format','yyyy.MM.dd-HH.mm.ss');
            temp = datestr(time,'yyyy.mm.dd-HH.MM.ss');
            name = ['calibration-' temp '.xlsx'];
            
            
            filter = {'*.xlsx';'*.csv'};
            [file,path] = uiputfile(filter,'BioControl 1.0',name);
            
            check = file == 0;
            if check
                return
            end
            
            if strcmpi(app.TempFolderPath.Value, path(1:(end-1)))
                path = path(1:(end-5));
            end
            
            name = fullfile(path, file);

            prompt = 'Comments:';
            dlgtitle = 'BioControl 1.0';
            dims = [7 70];
            comments = inputdlg(prompt,dlgtitle,dims);
            if ~isempty(char(comments))
                textlog{1} = 'Comments:';
                textlog{end + 1} = char(comments);
                textlog{end + 1} = '';
            else
                textlog{1} = 'Comments:';
                textlog{end + 1} = 'none';
                textlog{end + 1} = '';   
            end
            writecell(textlog, name)
            
            if isempty(app.FitParameters.Value)
                % Stores just data points
                datalog = {'Date','Time','';
                datestr(datetime('now','format','dd/MM/yyyy'),'dd/mm/yyyy'),...
                datestr(datetime('now','format','HH:mm:ss'),'HH:MM:ss'),'';
                    'LINEAR FIT ON LOG10 VALUES','','';
                    'Slope','Intercept','R2';
                    '','','';
                    'Concentration [ppm vol]','Intensity [ions/s]','STD Error []'};
            else
                % Stores data points and calibration data
                datalog = {'Date','Time','';
                datestr(datetime('now','format','dd/MM/yyyy'),'dd/mm/yyyy'),...
                datestr(datetime('now','format','HH:mm:ss'),'HH:MM:ss'),'';
                    'LINEAR FIT ON LOG10 VALUES','','';
                    'Slope','Intercept','R2';
                app.FitParameters.Value(1),app.FitParameters.Value(2),app.R2EditField.Value;
                    'Concentration [ppb vol]','Intensity [ions/s]','STD Error []'};
            end

            writecell(datalog, name, 'Range', 'A2:C7')

            co = (app.DataToFit.Value(1,:));
            av = (app.DataToFit.Value(2,:));
            er = (app.DataToFit.Value(3,:));
            data = [co' av' er'];
            
            textpoints = 'datalog_points = {';
            for i = 1:length(co)
                textpoints = [textpoints 'data(' num2str(i) ',1),data(' num2str(i) ',2),data(' num2str(i) ',3)'];
                if i ~= length(co)
                    textpoints = [textpoints ';'];
                end
            end
            textpoints = [textpoints '};'];
            
            eval(textpoints);

            writecell(datalog_points, name, 'WriteMode', 'append')
        end

        % Close request function: BioControl10CalibrationtabUIFigure
        function BioControl10CalibrationtabUIFigureCloseRequest(app, event)
            if ~isempty(app.TempFolderPath.Value)
                if 7 == exist(app.TempFolderPath.Value, 'dir')
                    cd(app.TempFolderPath.Value)
                    DeleteH5Content(app)
                    if isempty(ls)
                        cd ..
                        rmdir(app.TempFolderPath.Value)
                    end
                end
            end
            
            clear parameters
            global parameters
            parameters = app.FitParameters.Value;
            
            delete(app);
            fclose('all');
        end

        % Button pushed function: GeneratemockcalibrationButton
        function GeneratemockcalibrationButtonPushed(app, event)
            prompt = 'Mass [Da]';
            dlgtitle = 'Enter numerical value of mocked mass';
            dims = [1 100];
            output = inputdlg(prompt,dlgtitle,dims);
            
            check = isempty(output);
            if check
                return
            end
            %---------------------------
            CleargraphButtonPushed(app, event)
            %---------------------------
            if isempty(app.TempFolderPath.Value)
                mkdir('Temp')
                cd('Temp')
                app.TempFolderPath.Value = cd;
            else
                cd(app.TempFolderPath.Value)
                DeleteFolderContent(app)
            end
            %----------------------------
            user_mass = str2double(output{1});
            timelength_mock = 60;
            signal_intensity = 25;
            conc = 0;
            
            for i = 1:5
                filename = genaddh5calmock( ...
                    user_mass,timelength_mock,signal_intensity);
                
                [cumpeakprof,~,~] = geth5mrcumpeaks( ...
                    filename, ...
                    (user_mass + 1), ...
                    0.3);
                
                avg_output = mean(cumpeakprof);
                error = std(cumpeakprof);

                errorbar(app.UIAxes, ...
                conc, ...
                avg_output, ...
                error, ...
                'Marker', 'o', ...
                'MarkerFaceColor', 'b');
                
                new_point = [conc; avg_output; error];
                app.DataToFit.Value = [app.DataToFit.Value new_point];
                
                signal_intensity = signal_intensity + 500;
                conc = conc + 10;
            end
            %----------------------------
            app.NoPointsOnGraph.Value = app.NoPointsOnGraph.Value + 5;
            
            LinearfitButtonPushed(app, event)

        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create BioControl10CalibrationtabUIFigure and hide until all components are created
            app.BioControl10CalibrationtabUIFigure = uifigure('Visible', 'off');
            app.BioControl10CalibrationtabUIFigure.AutoResizeChildren = 'off';
            app.BioControl10CalibrationtabUIFigure.Position = [330 130 700 480];
            app.BioControl10CalibrationtabUIFigure.Name = 'BioControl 1.0 - Calibration tab';
            app.BioControl10CalibrationtabUIFigure.Resize = 'off';
            app.BioControl10CalibrationtabUIFigure.CloseRequestFcn = createCallbackFcn(app, @BioControl10CalibrationtabUIFigureCloseRequest, true);

            % Create UIAxes
            app.UIAxes = uiaxes(app.BioControl10CalibrationtabUIFigure);
            title(app.UIAxes, 'Calibration data')
            xlabel(app.UIAxes, 'Target compound concentration [ppb vol]')
            ylabel(app.UIAxes, 'Average intensity [ions/s] ')
            app.UIAxes.PlotBoxAspectRatio = [1.24148606811146 1 1];
            app.UIAxes.Position = [13 79 449 378];

            % Create PaircalibrationfileButton
            app.PaircalibrationfileButton = uibutton(app.BioControl10CalibrationtabUIFigure, 'push');
            app.PaircalibrationfileButton.ButtonPushedFcn = createCallbackFcn(app, @PaircalibrationfileButtonPushed, true);
            app.PaircalibrationfileButton.Position = [531 244 100 36];
            app.PaircalibrationfileButton.Text = {'Pair'; 'calibration file'};

            % Create CompoundconcentrationppbLabel
            app.CompoundconcentrationppbLabel = uilabel(app.BioControl10CalibrationtabUIFigure);
            app.CompoundconcentrationppbLabel.Position = [486 402 110 28];
            app.CompoundconcentrationppbLabel.Text = {'Compound'; 'concentration [ppb]'};

            % Create ComponentConcentrationEditField
            app.ComponentConcentrationEditField = uieditfield(app.BioControl10CalibrationtabUIFigure, 'numeric');
            app.ComponentConcentrationEditField.Position = [613 402 61 28];

            % Create CleargraphButton
            app.CleargraphButton = uibutton(app.BioControl10CalibrationtabUIFigure, 'push');
            app.CleargraphButton.ButtonPushedFcn = createCallbackFcn(app, @CleargraphButtonPushed, true);
            app.CleargraphButton.Position = [531 131 100 22];
            app.CleargraphButton.Text = 'Clear graph';

            % Create SlopeEditFieldLabel
            app.SlopeEditFieldLabel = uilabel(app.BioControl10CalibrationtabUIFigure);
            app.SlopeEditFieldLabel.HorizontalAlignment = 'right';
            app.SlopeEditFieldLabel.Position = [171 29 36 22];
            app.SlopeEditFieldLabel.Text = 'Slope';

            % Create SlopeEditField
            app.SlopeEditField = uieditfield(app.BioControl10CalibrationtabUIFigure, 'numeric');
            app.SlopeEditField.Editable = 'off';
            app.SlopeEditField.Position = [214 29 74 22];

            % Create LinearfitButton
            app.LinearfitButton = uibutton(app.BioControl10CalibrationtabUIFigure, 'push');
            app.LinearfitButton.ButtonPushedFcn = createCallbackFcn(app, @LinearfitButtonPushed, true);
            app.LinearfitButton.Position = [44 29 100 22];
            app.LinearfitButton.Text = 'Linear fit';

            % Create StorecalibrationButton
            app.StorecalibrationButton = uibutton(app.BioControl10CalibrationtabUIFigure, 'push');
            app.StorecalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @StorecalibrationButtonPushed, true);
            app.StorecalibrationButton.Position = [530 104 102 22];
            app.StorecalibrationButton.Text = 'Store calibration';

            % Create InterceptEditFieldLabel
            app.InterceptEditFieldLabel = uilabel(app.BioControl10CalibrationtabUIFigure);
            app.InterceptEditFieldLabel.HorizontalAlignment = 'right';
            app.InterceptEditFieldLabel.Position = [321 29 52 22];
            app.InterceptEditFieldLabel.Text = 'Intercept';

            % Create InterceptEditField
            app.InterceptEditField = uieditfield(app.BioControl10CalibrationtabUIFigure, 'numeric');
            app.InterceptEditField.Editable = 'off';
            app.InterceptEditField.Position = [380 29 74 22];

            % Create R2EditFieldLabel
            app.R2EditFieldLabel = uilabel(app.BioControl10CalibrationtabUIFigure);
            app.R2EditFieldLabel.HorizontalAlignment = 'right';
            app.R2EditFieldLabel.Position = [499 29 25 22];
            app.R2EditFieldLabel.Text = 'R2';

            % Create R2EditField
            app.R2EditField = uieditfield(app.BioControl10CalibrationtabUIFigure, 'numeric');
            app.R2EditField.Editable = 'off';
            app.R2EditField.Position = [536 29 74 22];

            % Create GeneratemockcalibrationButton
            app.GeneratemockcalibrationButton = uibutton(app.BioControl10CalibrationtabUIFigure, 'push');
            app.GeneratemockcalibrationButton.ButtonPushedFcn = createCallbackFcn(app, @GeneratemockcalibrationButtonPushed, true);
            app.GeneratemockcalibrationButton.Position = [531 203 100 36];
            app.GeneratemockcalibrationButton.Text = {'Generate mock'; 'calibration'};

            % Create CompoundcentralmassDaLabel
            app.CompoundcentralmassDaLabel = uilabel(app.BioControl10CalibrationtabUIFigure);
            app.CompoundcentralmassDaLabel.Position = [488 363 108 28];
            app.CompoundcentralmassDaLabel.Text = {'Compound central'; 'mass (Da)'};

            % Create ComponentCentralMassEditField
            app.ComponentCentralMassEditField = uieditfield(app.BioControl10CalibrationtabUIFigure, 'numeric');
            app.ComponentCentralMassEditField.Position = [613 363 61 28];

            % Create CompoundmassneighbourhoodDaLabel
            app.CompoundmassneighbourhoodDaLabel = uilabel(app.BioControl10CalibrationtabUIFigure);
            app.CompoundmassneighbourhoodDaLabel.Position = [488 323 112 28];
            app.CompoundmassneighbourhoodDaLabel.Text = {'Compound mass'; 'neighbourhood (Da)'};

            % Create ComponentMassNeighbourhoodEditField
            app.ComponentMassNeighbourhoodEditField = uieditfield(app.BioControl10CalibrationtabUIFigure, 'numeric');
            app.ComponentMassNeighbourhoodEditField.Position = [613 323 61 28];

            % Show the figure after all components are created
            app.BioControl10CalibrationtabUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = cal

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.BioControl10CalibrationtabUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.BioControl10CalibrationtabUIFigure)
        end
    end
end