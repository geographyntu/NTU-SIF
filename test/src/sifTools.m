classdef sifTools
    %SIFTOOLS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Data
        CalibrationFile
        LMPCalibrationFile
        DarkSpectrum
        CalibratedSpectrum
    end
    
    methods
        function obj = sifTools(currentpath)
            %SIFTOOLS Construct an instance of this class
            %   Detailed explanation goes here
            [D,dir]=uigetfile('*.txt','Choose Data','MultiSelect','on');
            selectfiles = size(D',1); Data = [];
            for i = 1:selectfiles
            if selectfiles == 1, filename{i,1} = [dir,D];else filename{i,1} = [dir,D{i}];end
            fid = fopen(filename{i}, 'r');
            % Initialize variables to store information above ">>>>>Begin Spectral Data<<<<<"
            theData.filename = D{i};
            theData.dateInfo = '';
            theData.userInfo = '';
            theData.spectrometerInfo = '';
            theData.triggerMode = '';
            theData.integrationTime = [];
            theData.scansToAverage = '';
            theData.darkCorrectionEnabled = '';
            theData.nonlinearityCorrectionEnabled = '';
            theData.boxcarWidth = '';
            theData.xAxisMode = '';
            theData.numPixels = '';
            % Read and store information above the spectral data
            while true
                line = fgetl(fid);
                if ischar(line)
                    if contains(line, '>>>>>Begin Spectral Data<<<<<')
                        break;  % Exit the loop when you reach the data section
                    end
                    if contains(line, 'Date:')
                        line=split(line,'Date: ');
                        theData.dateInfo = line{2};
                    elseif contains(line, 'User:')
                        theData.userInfo = line;
                    elseif contains(line, 'Spectrometer:')
                        theData.spectrometerInfo = line;
                    elseif contains(line, 'Trigger mode:')
                        theData.triggerMode = line;
                    elseif contains(line, 'Integration Time (sec):')
                        line=split(line,':');
                        theData.integrationTime = str2num(line{2});
                    elseif contains(line, 'Scans to average:')
                        theData.scansToAverage = line;
                    elseif contains(line, 'Electric dark correction enabled:')
                        theData.darkCorrectionEnabled = line;
                    elseif contains(line, 'Nonlinearity correction enabled:')
                        theData.nonlinearityCorrectionEnabled = line;
                    elseif contains(line, 'Boxcar width:')
                        theData.boxcarWidth = line;
                    elseif contains(line, 'XAxis mode:')
                        theData.xAxisMode = line;
                    elseif contains(line, 'Number of Pixels in Spectrum:')
                        theData.numPixels = line;
                    end
                else
                    break;  % Exit the loop when end of file is reached
                end
            end
            % Now, you can read and process the spectral data
            theData.Raw = fscanf(fid, '%f %f', [2, inf])';
            Data=[Data;{theData}];
            end
            obj.Data = Data;
            fclose(fid);
            cd(currentpath);
        end
        function outputArg = callCalibrationFile(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [D,dir]=uigetfile('*.cal','Choose Data','MultiSelect','on');
            selectfiles = size(D',1); CalibrationFile = [];
            for i = 1:selectfiles
            if selectfiles == 1, filename{i,1} = [dir,D];else filename{i,1} = [dir,D{i}];end
            fid = fopen(filename{i}, 'r');
            % Initialize variables to store information above ">>>>>Begin Spectral Data<<<<<"
            theData.filename = D{i};
            theData.dateInfo = '';
            theData.userInfo = '';
            theData.spectrometerInfo = '';
            theData.triggerMode = '';
            theData.integrationTime = '';
            theData.scansToAverage = '';
            theData.darkCorrectionEnabled = '';
            theData.nonlinearityCorrectionEnabled = '';
            theData.boxcarWidth = '';
            theData.xAxisMode = '';
            theData.numPixels = '';
            % Read and store information above the spectral data
            while true
                line = fgetl(fid);
                if ischar(line)
                    if contains(line, '[uJoule/count]')
                        break;  % Exit the loop when you reach the data section
                    end
                    if contains(line, 'Date:')
                        theData.dateInfo = line;
                    elseif contains(line, 'User:')
                        theData.userInfo = line;
                    elseif contains(line, 'Spectrometer:')
                        theData.spectrometerInfo = line;
                    elseif contains(line, 'Trigger mode:')
                        theData.triggerMode = line;
                    elseif contains(line, 'Integration Time (sec):')
                        theData.integrationTime = line;
                    elseif contains(line, 'Scans to average:')
                        theData.scansToAverage = line;
                    elseif contains(line, 'Electric dark correction enabled:')
                        theData.darkCorrectionEnabled = line;
                    elseif contains(line, 'Nonlinearity correction enabled:')
                        theData.nonlinearityCorrectionEnabled = line;
                    elseif contains(line, 'Boxcar width:')
                        theData.boxcarWidth = line;
                    elseif contains(line, 'XAxis mode:')
                        theData.xAxisMode = line;
                    elseif contains(line, 'Number of Pixels in Spectrum:')
                        theData.numPixels = line;
                    end
                else
                    break;  % Exit the loop when end of file is reached
                end
            end
            % Now, you can read and process the spectral data
            theData.Raw = fscanf(fid, '%f %f', [2, inf])';
            CalibrationFile=[CalibrationFile;{theData}];
            end
            obj.CalibrationFile = CalibrationFile;
            fclose(fid);
            outputArg=obj;
        end
        function outputArg = callLMPCalibrationFile(obj)
            [D,dir]=uigetfile('*.lmp','Choose Data','MultiSelect','off');
            filename = [dir,D];
            fid = fopen(filename, 'r');
            theData.Raw = fscanf(fid, '%f %f', [2, inf])';
            obj.LMPCalibrationFile = theData;
            fclose(fid);
            outputArg=obj;
        end
        function outputArg = callDarkSpectrum(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            [D,dir]=uigetfile('*.txt','Choose Data','MultiSelect','off');
            filename = [dir,D];
            fid = fopen(filename, 'r');
            % Initialize variables to store information above ">>>>>Begin Spectral Data<<<<<"
            theData.filename = D;
            theData.dateInfo = '';
            theData.userInfo = '';
            theData.spectrometerInfo = '';
            theData.triggerMode = '';
            theData.integrationTime = '';
            theData.scansToAverage = '';
            theData.darkCorrectionEnabled = '';
            theData.nonlinearityCorrectionEnabled = '';
            theData.boxcarWidth = '';
            theData.xAxisMode = '';
            theData.numPixels = '';
            % Read and store information above the spectral data
            while true
                line = fgetl(fid);
                if ischar(line)
                    if contains(line, '>>>>>Begin Spectral Data<<<<<')
                        break;  % Exit the loop when you reach the data section
                    end
                    if contains(line, 'Date:')
                        theData.dateInfo = line;
                    elseif contains(line, 'User:')
                        theData.userInfo = line;
                    elseif contains(line, 'Spectrometer:')
                        theData.spectrometerInfo = line;
                    elseif contains(line, 'Trigger mode:')
                        theData.triggerMode = line;
                    elseif contains(line, 'Integration Time (sec):')
                        theData.integrationTime = line;
                    elseif contains(line, 'Scans to average:')
                        theData.scansToAverage = line;
                    elseif contains(line, 'Electric dark correction enabled:')
                        theData.darkCorrectionEnabled = line;
                    elseif contains(line, 'Nonlinearity correction enabled:')
                        theData.nonlinearityCorrectionEnabled = line;
                    elseif contains(line, 'Boxcar width:')
                        theData.boxcarWidth = line;
                    elseif contains(line, 'XAxis mode:')
                        theData.xAxisMode = line;
                    elseif contains(line, 'Number of Pixels in Spectrum:')
                        theData.numPixels = line;
                    end
                else
                    break;  % Exit the loop when end of file is reached
                end
            end
            % Now, you can read and process the spectral data
            theData.Raw = fscanf(fid, '%f %f', [2, inf])';
            figure, plot(theData.Raw(:,1),theData.Raw(:,2));title('DarkSpec');
            obj.DarkSpectrum = theData;
            fclose(fid);
            outputArg=obj;
        end
        function outputArg = callCalibrateSpectrum(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            dat1=obj.Data;
            selectedfiles = length(obj.Data);
            for i = 1:selectedfiles
            Filename=split(obj.Data{i, 1}.filename,'_');
            Filename = strtrim(Filename);
            Filename=Filename{2,1};
            % Sample
            targetdat=i; %whichO2=mod(targetdat,2);
            Sample = obj.Data{targetdat, 1}.Raw; % leaf1
            mx = max(Sample(:,1));
            mn = min(Sample(:,1));
            
            % Dark Spectrum
            if (mn<760 && mx>760)
                % Resample of calibration
                int_calibrate = interp1(obj.CalibrationFile{1, 1}.Raw(:,1), obj.CalibrationFile{1, 1}.Raw(:,2), Sample(:,1), 'spline');
                %T = obj.Data{targetdat, 1}.integrationTime;
            else
                % Resample of calibration
                int_calibrate = interp1(obj.CalibrationFile{2, 1}.Raw(:,1), obj.CalibrationFile{2, 1}.Raw(:,2), Sample(:,1), 'spline');
                %T = obj.Data{targetdat-1, 1}.integrationTime;
                %T = obj.Data{targetdat, 1}.integrationTime;
            end
            
            % Dark spec
            DarkSpec = obj.DarkSpectrum.Raw;
            DarkSpec = interp1(DarkSpec(:,1),DarkSpec(:,2),Sample(:,1),'spline');% Resample
            DarkSpec=[Sample(:,1),DarkSpec];
            
            % Lamp spec
            LampSpec = obj.LMPCalibrationFile.Raw;
            LampSpec = interp1(LampSpec(:,1),LampSpec(:,2),Sample(:,1),'spline');% Resample
            LampSpec=[Sample(:,1),LampSpec];
            
            % Subtraction with dark
            int_raw = Sample(:,2) - DarkSpec(:,2);
            %int_raw = Sample(:,2) - LampSpec(:,2);
            
            % Wavelength
            wl_sample=Sample(:,1);
            
            % Calibrate
            T = obj.Data{targetdat, 1}.integrationTime; % s
            d = 355;                                    % um
            A = pi*(d*0.0001/2)^2;                      % cm^2
            dLp = (wl_sample(end) - wl_sample(1))/length(wl_sample);
            fit_AbsIrr = (int_raw.*int_calibrate)/(T*A*dLp);
            
            Calibrated = [wl_sample,fit_AbsIrr];
            if (mn<760 && mx>760)
                figure, plot(Calibrated(:,1),Calibrated(:,2)); ylim([0,1.2]); xlim([730,790]);
                title(obj.Data{i, 1}.filename);
            else
                figure, plot(Calibrated(:,1),Calibrated(:,2)); ylim([0,1.2]); xlim([660,720]);
                title(obj.Data{i, 1}.filename);
            end
            dat1{i,1}.Raw=Calibrated;
            end
            
            obj.CalibratedSpectrum = dat1;
            outputArg = obj;
        end
        function outputArg = callSIF(obj)
            
            outputArg=obj;
        end
    end
    methods(Static)
        function result = AddCalibrationFile(Raw)
            [file, path] = uigetfile('*.mat', 'Select calibration data file',"MultiSelect","on");
            if isequal(file, 0) || isequal(path, 0)
                result='Nan';
                return;
            end
            selectedFile = length(file);
            for i = 1:selectedFile
                lbl = split(file{i},["_","."]); lbl = lbl{2};
                loadedData = load(fullfile(path, file{i}));
                calFile(:,1) = loadedData.wavelengthData;
                calFile(:,2) = loadedData.IpCalData;
                Data.(lbl) = calFile;
            end
            Data.infoData = loadedData.infoData;
            DataMentah = Raw;
            Devices = fieldnames(DataMentah);

            wavSamplePlot=[];intSamplePlot=[];
            for i = 1:length(Devices)
                reSample(:,1) = DataMentah.(Devices{i})(:,1);
                reSample(:,2) = interp1(Data.(Devices{i})(:,1),Data.(Devices{i})(:,2),reSample(:,1),'spline');
                Data.(Devices{i}) = reSample;
                wavSamplePlot = [wavSamplePlot,reSample(:,1)];
                intSamplePlot = [intSamplePlot,reSample(:,2)];
                infoSamplePlot{i}=['Calibration file ',Devices{i}];
            end
            result = Data;
        end
        function result = AddDiffuserFile(Raw)
            [file, path] = uigetfile('*.*', 'Select diffuser data file'); 
            if isequal(file, 0) || isequal(path, 0)
                result='Nan';
                return;
            end
            loadedData = load(fullfile(path, file));
            DataMentah = Raw;
            Devices = fieldnames(DataMentah);
            wavSamplePlot=[];intSamplePlot=[];
            for i = 1:length(Devices)
                reSample(:,1) = DataMentah.(Devices{i})(:,1);
                reSample(:,2) = (interp1(loadedData(:,1),loadedData(:,2),reSample(:,1),'spline'))*0.01;
                Data.(Devices{i}) = reSample;
                wavSamplePlot = [wavSamplePlot,reSample(:,1)];
                intSamplePlot = [intSamplePlot,reSample(:,2)];
                infoSamplePlot{i}=['Diffuser ',Devices{i}];
            end
            result = Data;
        end
        function result = AddGratingFile(Raw)
            [file, path] = uigetfile('*.*', 'Select grating data file',"MultiSelect","on"); 
            if isequal(file, 0) || isequal(path, 0)
                result='Nan';
                return;
            end
            selectedFile = length(file);
            for i = 1:selectedFile
                lbl = split(file{i},'.'); lbl = lbl{1, 1}  ;
                loadedData = load(fullfile(path, file{i}));
                if contains(file{i}, 'H6')
                    dvcname = ['QEP03935','QEP03916'];
                    Data.QEP03935 = loadedData;
                    Data.QEP03916 = loadedData;
                elseif contains(file{i}, 'H11') % Perlu edit lg kalo ada 4 device
                    dvcname = ['QEP03410','QEP03375'];
                    Data.QEP03410 = loadedData;
                    Data.QEP03375 = loadedData;
                end
            end
            
            DataMentah = Raw;
            Devices = fieldnames(DataMentah);
            wavSamplePlot=[];intSamplePlot=[];
            for i = 1:length(Devices)
                reSample(:,1) = DataMentah.(Devices{i})(:,1);
                reSample(:,2) = interp1(Data.(Devices{i})(:,1),Data.(Devices{i})(:,2),reSample(:,1),'spline');
                Data.(Devices{i}) = reSample;
                wavSamplePlot = [wavSamplePlot,reSample(:,1)];
                intSamplePlot = [intSamplePlot,reSample(:,2)];
                infoSamplePlot{i}=['Grating ',Devices{i}];
            end
            result = Data;
        end
        function result = AddQuantumEffFile(Raw)
            [file, path] = uigetfile('*.*', 'Select quantum eff. data file'); 
            if isequal(file, 0) || isequal(path, 0)
                result='Nan';
                return;
            end
            loadedData = load(fullfile(path, file));
            DataMentah = Raw;
            Devices = fieldnames(DataMentah);
            wavSamplePlot=[];intSamplePlot=[];
            for i = 1:length(Devices)
                reSample(:,1) = DataMentah.(Devices{i})(:,1);
                reSample(:,2) = interp1(loadedData(:,1),loadedData(:,2),reSample(:,1),'spline');
                Data.(Devices{i}) = reSample;
                wavSamplePlot = [wavSamplePlot,reSample(:,1)];
                intSamplePlot = [intSamplePlot,reSample(:,2)];
                infoSamplePlot{i}=['QE ',Devices{i}];
            end
            result = Data;
        end

        function result = dataResample(loadedData)
            % Resample
            wl_ref  = loadedData.wavelengthData;
            int_ref = loadedData.spectralData;
            if size(loadedData.infoData,2) == 3
                tools = loadedData.infoData(2,:);
                for i = 1:size(loadedData.infoData,2) % O2A
                    if strcmpi(tools{i},'QEP03935')
                        index = find(cellfun(@(x) isequal(x, 'QEP03375'), tools));
                        int_ref(:,i-1) = interp1(wl_ref(:,i-1),int_ref(:,i-1),wl_ref(:,index-1),'spline');% Resample
                        wl_ref(:,i-1) = wl_ref(:,index-1);
                    elseif strcmpi(tools{i},'QEP03916') % O2B
                        index = find(cellfun(@(x) isequal(x, 'QEP03410'), tools));
                        int_ref(:,i-1) = interp1(wl_ref(:,i-1),int_ref(:,i-1),wl_ref(:,index-1),'spline');% Resample
                        wl_ref(:,i-1) = wl_ref(:,index-1);
                    end
                end
                loadedData.wavelengthData = wl_ref;
                loadedData.spectralData = int_ref;
            elseif size(loadedData.infoData,2) == 5
                tools = loadedData.infoData(2,:);
                indexQEP03935 = find(cellfun(@(x) isequal(x, 'QEP03935'), tools));
                indexQEP03375 = find(cellfun(@(x) isequal(x, 'QEP03375'), tools));
                int_ref(:,indexQEP03935-1) = interp1(wl_ref(:,indexQEP03935-1),int_ref(:,indexQEP03935-1),wl_ref(:,indexQEP03375-1),'spline');% Resample
                wl_ref(:,indexQEP03935-1) = wl_ref(:,indexQEP03375-1);
                indexQEP03916 = find(cellfun(@(x) isequal(x, 'QEP03916'), tools));
                indexQEP03410 = find(cellfun(@(x) isequal(x, 'QEP03410'), tools));
                int_ref(:,indexQEP03916-1) = interp1(wl_ref(:,indexQEP03916-1),int_ref(:,indexQEP03916-1),wl_ref(:,indexQEP03410-1),'spline');% Resample
                wl_ref(:,indexQEP03916-1) = wl_ref(:,indexQEP03410-1);
                % save
                loadedData.wavelengthData = wl_ref;
                loadedData.spectralData = int_ref;
            end
            result=loadedData;
        end

        % METHODS-function level 1
        function result = methodsfld_O2A(L,E,varargin)
            sleft = varargin{1}; sright = varargin{2};
            pE = [E(sleft,:);E(sright,:)];
            pL = [L(sleft,:);L(sright,:)];
            PIL = sifTools.determinePointsFLD_O2A(L,sleft,sright);
            PIE = sifTools.determinePointsFLD_O2A(E,sleft,sright);
            PointsFLD = sifTools.sFLD(PIE,PIL);
            %figure, plot(L(:,1),L(:,2),'-b'); hold on; title('sFLD');
            %plot(E(:,1),E(:,2),'-r');
            %plot(pE(:,1),pE(:,2),'*g'); plot(pL(:,1),pL(:,2),'*g');
            %plot(PIE(:,1),PIE(:,2),'sk','MarkerFaceColor','k');
            %plot(PIL(:,1),PIL(:,2),'sk','MarkerFaceColor','k');
            %legend('L','E','pointE','pointL','Point');
            result.L = L;
            result.E = E;
            result.pL = pL;
            result.pE = pE;
            result.PIL = PIL;
            result.PIE = PIE;
            result.PointsFLD = PointsFLD;
        end
        function result = methods3fld_O2A(L,E,varargin)
            sleft = varargin{1}; sright = varargin{2};
            pE = [E(sleft,:);E(sright,:)];
            pL = [L(sleft,:);L(sright,:)];
            [PIE,LE] = sifTools.determinePoint3FLD_O2A(E,sleft,sright);
            [PIL,LL] = sifTools.determinePoint3FLD_O2A(L,sleft,sright);
            PointsFLD = sifTools.sFLD(PIE,PIL);
            %figure, plot(L(:,1),L(:,2),'-b'); hold on; title('3FLD');
            %plot(E(:,1),E(:,2),'-r');
            %plot(pE(:,1),pE(:,2),'*g'); plot(pL(:,1),pL(:,2),'*g');
            %plot(PIE(:,1),PIE(:,2),'sk','MarkerFaceColor','k');
            %plot(PIL(:,1),PIL(:,2),'sk','MarkerFaceColor','k');
            %plot(LE(:,1),LE(:,2),'m'); plot(LL(:,1),LL(:,2),'m');
            %legend('L','E','pointE','pointL','Point');
            %result = PointsFLD;
            result.L = L;
            result.E = E;
            result.pL = pL;
            result.pE = pE;
            result.PIL = PIL;
            result.PIE = PIE;
            result.LE = LE;
            result.LL = LL;
            result.Point3FLD = PointsFLD;
        end
        function result = methodsifld_O2A(L,E,varargin)
            sleft = varargin{1}; sright = varargin{2};
            pE = [E(sleft,:);E(sright,:)];
            pL = [L(sleft,:);L(sright,:)];
            [PIE,LE] = sifTools.determinePointiFLD_O2A(E,sleft,sright);
            [PIL] = sifTools.determinePointsFLD_O2A(L,sleft,sright);
            [ppe] = sifTools.determinePointsFLD_O2A(E,sleft,sright);
            [oRap, sRap] = sifTools.Rap(L,E);
            PointiFLD = sifTools.iFLD(L,E,PIE,LE,PIL,oRap,sleft,sright);
            %figure, plot(L(:,1),L(:,2),'-b'); hold on; title('iFLD');
            %plot(E(:,1),E(:,2),'-r');
            %plot(pE(:,1),pE(:,2),'*g'); plot(pL(:,1),pL(:,2),'*g');
            %plot(PIE(:,1),PIE(:,2),'sk','MarkerFaceColor','k');
            %plot(PIL(:,1),PIL(:,2),'sk','MarkerFaceColor','k');
            %plot(LE(:,1),LE(:,2),'m'); plot(ppe(1,1),ppe(1,2),'sk','MarkerFaceColor','k');
            %legend('L','E','pointE','pointL','Point');
            %result = {PointiFLD,oRap,sRap};
            result.L = L;
            result.E = E;
            result.pL = pL;
            result.pE = pE;
            result.PIL = PIL;
            result.PIE = PIE;
            result.LE = LE;
            result.ppe = ppe;
            result.oRap = oRap;
            result.sRap = sRap;
            result.PointiFLD = PointiFLD;
        end
        %__________________________________________________________________
        function result = methodsfld_O2B(L,E,varargin)
            sleft = varargin{1}; sright = varargin{2};
            pE = [E(sleft,:);E(sright,:)];
            pL = [L(sleft,:);L(sright,:)];
            PIE = sifTools.determinePointsFLD_O2B(E,sleft,sright);
            PIL = sifTools.determinePointsFLD_O2B(L,sleft,sright);
            PointsFLD = sifTools.sFLD(PIE,PIL);
            %figure, plot(L(:,1),L(:,2),'-b'); hold on; title('sFLD');
            %plot(E(:,1),E(:,2),'-r');
            %plot(pE(:,1),pE(:,2),'*g'); plot(pL(:,1),pL(:,2),'*g');
            %plot(PIE(:,1),PIE(:,2),'sk','MarkerFaceColor','k');
            %plot(PIL(:,1),PIL(:,2),'sk','MarkerFaceColor','k');
            %legend('L','E','pointE','pointL','Point');
            result.L = L;
            result.E = E;
            result.pL = pL;
            result.pE = pE;
            result.PIL = PIL;
            result.PIE = PIE;
            result.PointsFLD = PointsFLD;
        end
        function result = methods3fld_O2B(L,E,varargin)
            sleft = varargin{1}; sright = varargin{2};
            pE = [E(sleft,:);E(sright,:)];
            pL = [L(sleft,:);L(sright,:)];
            [PIE,LE] = sifTools.determinePoint3FLD_O2B(E,sleft,sright);
            [PIL,LL] = sifTools.determinePoint3FLD_O2B(L,sleft,sright);
            PointsFLD = sifTools.sFLD(PIE,PIL);
            %figure, plot(L(:,1),L(:,2),'-b'); hold on; title('3FLD');
            %plot(E(:,1),E(:,2),'-r');
            %plot(pE(:,1),pE(:,2),'*g'); plot(pL(:,1),pL(:,2),'*g');
            %plot(PIE(:,1),PIE(:,2),'sk','MarkerFaceColor','k');
            %plot(PIL(:,1),PIL(:,2),'sk','MarkerFaceColor','k');
            %plot(LE(:,1),LE(:,2),'m'); plot(LL(:,1),LL(:,2),'m');
            %legend('L','E','pointE','pointL','Point');
            result.L = L;
            result.E = E;
            result.pL = pL;
            result.pE = pE;
            result.PIL = PIL;
            result.PIE = PIE;
            result.LE = LE;
            result.LL = LL;
            result.Point3FLD = PointsFLD;
        end
        function result = methodsifld_O2B(L,E,varargin)
            sleft = varargin{1}; sright = varargin{2};
            pE = [E(sleft,:);E(sright,:)];
            pL = [L(sleft,:);L(sright,:)];
            [PIE,LE] = sifTools.determinePointiFLD_O2B(E,sleft,sright);
            [PIL] = sifTools.determinePointsFLD_O2B(L,sleft,sright);
            [ppe] = sifTools.determinePointsFLD_O2B(E,sleft,sright);
            [oRap, sRap] = sifTools.Rap(L,E);
            PointiFLD = sifTools.iFLD(L,E,PIE,LE,PIL,sRap,sleft,sright);
            %figure, plot(L(:,1),L(:,2),'-b'); hold on; title('iFLD'); 
            %plot(E(:,1),E(:,2),'-r'); xlim([670,720]);ylim([0,1]);
            %plot(pE(:,1),pE(:,2),'*g'); plot(pL(:,1),pL(:,2),'*g');
            %plot(PIE(:,1),PIE(:,2),'sk','MarkerFaceColor','k');
            %plot(PIL(:,1),PIL(:,2),'sk','MarkerFaceColor','k');
            %plot(LE(:,1),LE(:,2),'m'); plot(ppe(1,1),ppe(1,2),'sk','MarkerFaceColor','k');
            %legend('L','E','pointE','pointL','Point');
            result.L = L;
            result.E = E;
            result.pL = pL;
            result.pE = pE;
            result.PIL = PIL;
            result.PIE = PIE;
            result.LE = LE;
            result.ppe = ppe;
            result.oRap = oRap;
            result.sRap = sRap;
            result.PointiFLD = PointiFLD;
        end
        % Sub-Function level 2
        function resl = determinePointsFLD_O2A(E,sleft,sright)
            pup = E(sleft,:); idx = find(pup(:,2) == max(pup(:,2)));
            pup = E(idx+(sleft(1)-1),:);
            idx = find(E(:,2) == min(E(sleft(end):sright(1),2)));
            pbtm = E(idx,:);
            resl = [pup;pbtm];
        end
        function [resl,linear] = determinePoint3FLD_O2A(E,sleft,sright)
            pupl = E(sleft,:); idx = find(pupl(:,2) == max(pupl(:,2)));
            pupl = E(idx+(sleft(1)-1),:);
            pupr = E(sright,:); idx = find(pupr(:,2) == max(pupr(:,2)));
            pupr = E(idx+(sright(1)-1),:);
            dat  = [pupl;pupr];
            %p = polyfit(dat(:,1),dat(:,2),1); h = polyval(p,E(:,1)); h = [E(:,1),h];
            %---------------------------------------------------------------------------------GOOD------------------------------------------------------------------------
            h = E(:,1); h(:,2) = interp1(dat(:,1),dat(:,2),E(:,1),'linear'); % 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', or 'spline'
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            idx = find(E(:,2) == min(E(sleft(end):sright(1),2)));
            pbtm = E(idx,:);
            pup = h(idx,:);
            resl = [pup;pbtm];
            linear = h;
        end
        function [resl,linear] = determinePointiFLD_O2A(E,sleft,sright)
            pupl = E(sleft,:); idx = find(pupl(:,2) == max(pupl(:,2)));
            pupl = E(idx+(sleft(1)-1),:);
            pupr = E(sright,:); idx = find(pupr(:,2) == max(pupr(:,2)));
            pupr = E(idx+(sright(1)-1),:);
            dat  = [pupl;pupr];
            %p = polyfit(dat(:,1),dat(:,2),1); h = polyval(p,E(:,1)); h = [E(:,1),h];
            %---------------------------------------------------------------------------------GOOD------------------------------------------------------------------------
            h = E(:,1); h(:,2) = interp1(dat(:,1),dat(:,2),E(:,1),'linear'); % 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', or 'spline'
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------
            
            idx = find(E(:,2) == min(E(sleft(end):sright(1),2)));
            pbtm = E(idx,:);
            pup = h(idx,:);
            resl = [pup;pbtm];
            linear = h;
        end
        %__________________________________________________________________
        function resl = determinePointsFLD_O2B(E,sleft,sright)
            pup = E(sleft,:); idx = find(pup(:,2) == max(pup(:,2)));
            pup = E(idx+(sleft(1)-1),:);
            idx = find(E(:,2) == min(E(sleft(end):sright(1),2)));
            pbtm = E(idx,:);
            resl = [pup;pbtm];
        end
        function [resl,linear] = determinePoint3FLD_O2B(E,sleft,sright)
            pupl = E(sleft,:); idx = find(pupl(:,2) == max(pupl(:,2)));
            pupl = E(idx+(sleft(1)-1),:);
            pupr = E(sright,:); idx = find(pupr(:,2) == max(pupr(:,2)));
            pupr = E(idx+(sright(1)-1),:);
            %dat  = [pupl;pupr];
            dat  = [E(sleft,:);E(sright,:)];
            %h = E(:,1); h(:,2) = spline(dat(:,1),dat(:,2),E(:,1)); % pchip, spline, makima %
            %h = E(:,1); h(:,2) = pchip(dat(:,1),dat(:,2),E(:,1)); % pchip, spline, makima %
            %---------------------------------------------------------------------------------GOOD------------------------------------------------------------------------
            h = E(:,1); h(:,2) = interp1(dat(:,1),dat(:,2),E(:,1),'linear'); % 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', or 'spline'
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------
            %p = polyfit(dat(:,1),dat(:,2),4); h = polyval(p,E(:,1)); h = [E(:,1),h];
            idx = find(E(:,2) == min(E(sleft(end):sright(1),2)));
            pbtm = E(idx,:);
            pup = h(idx,:);
            resl = [pup;pbtm];
            linear = h;
        end
        function [resl,linear] = determinePointiFLD_O2B(E,sleft,sright)
            pupl = E(sleft,:); idx = find(pupl(:,2) == max(pupl(:,2)));
            pupl = E(idx+(sleft(1)-1),:);
            pupr = E(sright,:); idx = find(pupr(:,2) == max(pupr(:,2)));
            pupr = E(idx+(sright(1)-1),:);
            %dat  = [pupl;pupr];
            dat  = [E(sleft,:);E(sright,:)];
            %h = E(:,1); h(:,2) = spline(dat(:,1),dat(:,2),E(:,1)); % pchip, spline, makima %
            %h = E(:,1); h(:,2) = pchip(dat(:,1),dat(:,2),E(:,1)); % pchip, spline, makima %
            %---------------------------------------------------------------------------------GOOD------------------------------------------------------------------------
            h = E(:,1); h(:,2) = interp1(dat(:,1),dat(:,2),E(:,1),'linear'); % 'linear', 'nearest', 'next', 'previous', 'pchip', 'cubic', 'v5cubic', 'makima', or 'spline'
            %-------------------------------------------------------------------------------------------------------------------------------------------------------------
            % p = polyfit(dat(:,1),dat(:,2),1); h = polyval(p,E(:,1)); h = [E(:,1),h];
            idx = find(E(:,2) == min(E(sleft(end):sright(1),2)));
            pbtm = E(idx,:);
            pup = h(idx,:);
            resl = [pup;pbtm];
            linear = h;
        end
        % Sub-Function level 3
        function PointiFLD = iFLD(L,E,PIE,LE,PIL,sRap,varargin)
            % L     = Radiance Up
            % E     = Irradiance Down
            % PIE   = E - Shoulder interpolated and Bottom
            % LE    = interpolated-array shoulder
            % PIL   = L - Shoulder and Bottom
            % sRap  = scalled Apparent Reflectance
            rap = [sRap(varargin{1},:);sRap(varargin{2},:)];
            h = E(:,1); h(:,2) = pchip(rap(:,1),rap(:,2),E(:,1)); % pchip, spline, makima %
            
            %p = polyfit(rap(:,1),rap(:,2),2); h = polyval(p,sRap(:,1)); h = [sRap(:,1),h];
            
            %figure, plot(sRap(:,1),sRap(:,2)); hold on; title('Apparent Reflectance'); 
            %idx = find(L(:,2) == PIL(1,2));
            %plot(rap(:,1),rap(:,2),'og');
            %plot(sRap(idx,1),sRap(idx,2),'*r');
            %idx = find(L(:,2) == PIL(2,2));
            %plot(h(idx,1),h(idx,2),'sr');
            %plot(h(:,1),h(:,2),'-r');
            %legend('Rap','Point','Shoulder','Bottom','Regression');
            rappout = find(L(:,2) == PIL(1,2));
            rappin = find(L(:,2) == PIL(2,2));
            Sh = sifTools.determinePointsFLD_O2A(E,varargin{1},varargin{2});
            ein = find(E(:,2) == PIE(2,2));
            
            AlfaR = sRap(rappout,2)/h(rappin,2);
            AlfaF = AlfaR*(Sh(1,2)/LE(ein,2));
            
            subtractedRap = sRap(:,1);
            subtractedRap(:,2) = sRap(:,2) - h(:,2);
            pointE = sifTools.determinePointsFLD_O2A(E,varargin{1},varargin{2});
            pointL = sifTools.determinePointsFLD_O2A(L,varargin{1},varargin{2});
            PointiFLD = sifTools.lagiiFLD(pointE,pointL,AlfaR,AlfaF);
            PointiFLD = {PointiFLD,subtractedRap};
        end
        function resl = lagiiFLD(E,L,AlfaR,AlfaF)
            F = ((AlfaR*E(1,2)*L(2,2))-(L(1,2)*E(2,2)))/((AlfaR*E(1,2))-(E(2,2)*AlfaF));
            R = pi*((L(2,2)-F)/E(2,2));
            resl=[F;R];
        end
        function resl = sFLD(E,L)
            F = ((E(1,2)*L(2,2))-(L(1,2)*E(2,2)))/(E(1,2)-E(2,2));
            R = pi*((L(1,2)-L(2,2))/(E(1,2)-E(2,2)));
            resl=[F;R];
        end
        function [oRap, sRap] = Rap(L,E)
            % Apparent Reflectance
            RapO2A(:,1) = L(:,1);
            RapO2A(:,2) = L(:,2)./E(:,2);
            oRap = RapO2A;
            %figure, plot(RapO2A(:,1),RapO2A(:,2),'-b'); title('Apparent Reflectance'); legend('Rap');
            % Rescaled
            %RapO2A(:,2) = rescale(RapO2A(:,2), 0, 0.02);
            %figure, plot(RapO2A(:,1),RapO2A(:,2),'-b'); title('Apparent Reflectance (scaled)'); legend('Rap scaled');
            sRap = RapO2A;
        end
        
    end
end


