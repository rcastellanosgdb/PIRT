classdef TestPIRTSynthetic < matlab.unittest.TestCase
    % TESTPIRTSYNTHETIC unit tests of the PIRT MATLAB toolbox on synthetic data.
    %
    %   Run from the "matlab" folder with:
    %       results = runtests('tests/TestPIRTSynthetic.m')
    %
    %   The tests mirror the pytest suite of the Python version
    %   (python/tests) so that both implementations are checked against the
    %   same analytical references.
    %
    %   NOTE: this file was written during the reorganisation of the
    %   repository (v1.1) on a machine without a valid MATLAB licence and has
    %   therefore NOT been executed. Please report any failure.

    methods (TestClassSetup)
        function addPaths(testCase)
            here = fileparts(mfilename('fullpath'));
            addpath(fullfile(here,'..','src'));
            addpath(fullfile(here,'..','utils'));
        end
    end

    methods (Test)

        function test_sgolay32_quadratic_field(testCase)
            dx = 0.2; dy = 0.3; dt = 0.05;
            [X,Y,T] = meshgrid((0:29)*dx,(0:19)*dy,(0:9)*dt); % rows=y, cols=x
            F = 3 + 2*X - Y + 0.5*T + X.*Y + 4*X.^2 - 2*Y.^2 + 3*T.^2 + X.*T;
            [Ff,dFdt,d2Fdx2,d2Fdy2] = sgolay32_filter(F,[5,5,3],[dx,dy,dt]);
            testCase.verifySize(Ff,[20,30,8]);
            inner = {3:18, 3:28, 1:8};
            testCase.verifyEqual(Ff(inner{:}), F(3:18,3:28,2:9), 'AbsTol', 1e-9);
            testCase.verifyEqual(d2Fdx2(inner{:}), 8*ones(16,26,8), 'AbsTol', 1e-8);
            testCase.verifyEqual(d2Fdy2(inner{:}), -4*ones(16,26,8), 'AbsTol', 1e-8);
            expected = 0.5 + 6*T(:,:,2:9) + X(:,:,2:9);
            testCase.verifyEqual(dFdt(inner{:}), expected(3:18,3:28,:), 'AbsTol', 1e-8);
        end

        function test_derivative_fd_convention(testCase)
            % x along the columns (dx), y along the rows (dy)
            dx = 0.2; dy = 0.3; dt = 0.05;
            [X,Y,T] = meshgrid((0:14)*dx,(0:11)*dy,(0:7)*dt);
            F = 1 + 4*X.^2 - 2*Y.^2 + 0.5*T + 3*T.^2 + X.*T;
            [d2x,d2y,dFdt] = Derivative_FD(F,1,1,dx,dy,dt);
            testCase.verifyEqual(d2x, 8*ones(size(F)), 'AbsTol', 1e-8);
            testCase.verifyEqual(d2y, -4*ones(size(F)), 'AbsTol', 1e-8);
            testCase.verifyEqual(dFdt, 0.5 + 6*T + X, 'AbsTol', 1e-8);
        end

        function test_pod_filter_low_rank(testCase)
            rng(0);
            A = randn(30,40,3); B = randn(3,50);
            X = reshape(reshape(A,[],3)*B,[30,40,50]) + 5;
            Xf = POD_filter(X,'Nmod',3);
            testCase.verifyEqual(Xf, X, 'AbsTol', 1e-9);
        end

        function test_svht_coefficient(testCase)
            % Gavish & Donoho (2014): lambda*(1) = 4/sqrt(3), omega(1) ~ 2.858
            testCase.verifyEqual(optimal_SVHT_coef(1,1), 4/sqrt(3), 'RelTol', 1e-10);
            testCase.verifyEqual(optimal_SVHT_coef(1,0), 2.858, 'AbsTol', 2e-3);
        end

        function test_wiener3_constant_field(testCase)
            X = 300*ones(12,12,4);
            testCase.verifyEqual(wiener3(X,[3 3 1],0.1), X, 'AbsTol', 1e-9);
        end

        function test_cutoff_kernels_symmetric_and_complementary(testCase)
            for sz = {[20,30],[21,31],[20,31]}
                n = sz{1}(1); m = sz{1}(2);
                Hl = lowpass_kernel(n,m,0.4,0.6); Hh = highpass_kernel(n,m,0.4,0.6);
                testCase.verifyTrue(all(Hl(:)+Hh(:) >= 1));
                % DC bin of fftshift kept by the low-pass, removed by the high-pass
                testCase.verifyEqual(Hl(floor(n/2)+1,floor(m/2)+1), single(1));
                testCase.verifyEqual(Hh(floor(n/2)+1,floor(m/2)+1), single(0));
                % Hermitian symmetry H(k) = H(-k)
                Hs = circshift(flip(flip(Hl,1),2),[mod(n+1,2), mod(m+1,2)]);
                testCase.verifyEqual(Hs, Hl);
            end
            % a low-pass filtered real field is real (not folded by abs)
            [xx,~] = meshgrid(0:95,0:63);
            F = 10 + cos(2*pi*2*xx/96) + cos(2*pi*40*xx/96);
            Ff = Spatial_Cutoff_Filter(F,0.5);
            testCase.verifyEqual(Ff, 10 + cos(2*pi*2*xx/96), 'AbsTol', 1e-9);
            Fh = Spatial_Cutoff_Filter(F,0.5,'high');
            testCase.verifyEqual(Fh, cos(2*pi*40*xx/96), 'AbsTol', 1e-9);
        end

        function test_heat_transfer_uniform_case(testCase)
            HFS.s = 1e-5; HFS.rho = 7900; HFS.cp = 460; HFS.k = 14.7e-5;
            HFS.H = 0.276; HFS.W = 0.1; HFS.epsilon = 0.95; HFS.sides = 2;
            HFS.s_paint = 42e-6; HFS.cp_paint = 3061.5; HFS.rho_paint = 1261.175; HFS.lambda_paint = 1.38;
            Conditions.V = 2.227; Conditions.I = 7.5; Conditions.Uinf = 1.233;
            Conditions.Tamb = [294.75 294.75]; Conditions.L = 0.01; Conditions.dt = 1/253;
            Conditions.dx = 2.27e-4; Conditions.dy = 2.27e-4; Conditions.rhoinf = 1.2; Conditions.cpinf = 1005;
            Thot = 300*ones(5,6,4); Tcold = 294.35*ones(5,6);
            obj = PIRT('Thot',Thot,'Tcold',Tcold,'CalculateHeatTransfer','h','Nu','St', ...
                'TimeDer','SpatialDer','HFS',HFS,'Conditions',Conditions);
            obj = obj.go();
            qj = 2.227*7.5/(0.276*0.1);
            qr = 2*5.67e-8*0.95*(300^4-294.75^4);
            h = (qj-qr)/(300-294.35);
            testCase.verifyEqual(obj.result.h, h*ones(5,6,4), 'RelTol', 1e-10);
            Tfilm = (300+294.75)/2;
            kair = 1.5207E-11*Tfilm^3-4.8574E-08*Tfilm^2+1.0184E-04*Tfilm-3.9333E-04;
            testCase.verifyEqual(obj.result.Nu, h*0.01/kair*ones(5,6,4), 'RelTol', 1e-10);
            testCase.verifyEqual(obj.result.St, h/(1.2*1005*1.233)*ones(5,6,4), 'RelTol', 1e-10);
        end
    end
end
