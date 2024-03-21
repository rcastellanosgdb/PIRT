classdef TestPIRT < matlab.unittest.TestCase
    properties
        Thot
        Tcold
        HFS
        Conditions
        extradata
    end
    methods (TestClassSetup)
        % Shared setup for the entire test class
    end

    methods (TestMethodSetup)
        % Setup for each test
        function loaddata(testCase)
            addpath("src/")
            addpath("utils/")
            addpath("tests/")
            % Load the test testCase.Conditions
            load('TestConditions.mat');
            load('Resolution.mat');
            % Load the temperature maps
            load('Tcold.mat');
            load('Thot.mat');

            % Crop the images: Could be done in PIRT by adding PIRT('Crop',Crop)
            Crop = [86 242;12 186]; %[x1 x2, y1, y2]
            testCase.Thot = Timage_hot(Crop(2,1):Crop(2,2),Crop(1,1):Crop(1,2),:);
            testCase.Tcold = Timage_cold(Crop(2,1):Crop(2,2),Crop(1,1):Crop(1,2),:);

            % Aquisition frequency
            testCase.extradata.dt = dt;
            testCase.extradata.dx = 0.001/dx;
            testCase.extradata.dy = 0.001/dy;

            % Some data is not included in .mat and was introduced manually
            testCase.Conditions.V = V;
            testCase.Conditions.I = I;
            testCase.Conditions.Uinf = Uinf;
            testCase.Conditions.L = 17e-3;
            testCase.Conditions.Tamb = Tamb([1,3]);

            % testCase.HFS.Type = 'PCB';
            testCase.HFS.s = 0.2e-3;
            testCase.HFS.rho = 8000;
            testCase.HFS.cp = 500;
            testCase.HFS.W = 0.15;
            testCase.HFS.H = 0.15;
            testCase.HFS.epsilon = 0.95;
            testCase.HFS.lambdax = 1.6;
            testCase.HFS.lambday = 3.1; %[W/mk]
        end
    end

    methods (Test)
        % Test methods

        function test_Conditions_data_save(testCase)
            object = PIRT('Thot',testCase.Thot,'CalculateHeatTransfer','h',...
                'Conditions',testCase.Conditions);
            Conditions = testCase.Conditions;
            Conditions.L_char = Conditions.L;
            Conditions = rmfield(Conditions,'L');

            testCase.verifyEqual(object.HeatTransfer_params.conditions,...
                Conditions,'Something ocurred when saving the data')
        end

        function test_HFS_data_save(testCase)
            object = PIRT('Thot',testCase.Thot,'CalculateHeatTransfer','h',...
                'HFS',testCase.HFS);
            HFS = testCase.HFS;
            HFS.A = HFS.W*HFS.H;
            HFS.m_plate = HFS.rho*HFS.s*HFS.A;
            HFS.Type = 'Foil';

            testCase.verifyEqual(object.HeatTransfer_params.HFS,...
                HFS,'Something ocurred when saving the data')
        end
    end

end