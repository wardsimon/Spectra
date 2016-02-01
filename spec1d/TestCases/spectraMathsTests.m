classdef spectraMathsTests < matlab.unittest.TestCase
    %SPECTRAMATHSTESTS Summary of this class goes here
    %   Detailed explanation goes here
    properties (Constant)
        N = 1E2;
    end
    properties
        s
        ref
        p
    end
    
    methods (Test)
        function plusSolution(testCase)
            s_test = struct(plus(testCase.s,testCase.p));
            testCase.verifyEqual(s_test.y,testCase.ref.plusS(1).y);
            testCase.verifyEqual(s_test.e,testCase.ref.plusS(1).e);
            s_test = struct(plus(testCase.p,testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.plusS(2).y);
            testCase.verifyEqual(s_test.e,testCase.ref.plusS(2).e);
            s_test = struct(plus(testCase.s,testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.plusS(3).y);
            testCase.verifyEqual(s_test.e,testCase.ref.plusS(3).e);
        end
        function timesSolution(testCase)
            s_test = struct(times(testCase.s,testCase.p));
            testCase.verifyEqual(s_test.y,testCase.ref.timesS(1).y);
            testCase.verifyEqual(s_test.e,testCase.ref.timesS(1).e);
            s_test = struct(times(testCase.p,testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.timesS(2).y);
            testCase.verifyEqual(s_test.e,testCase.ref.timesS(2).e);
            s_test = struct(times(testCase.s,testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.timesS(3).y);
            testCase.verifyEqual(s_test.e,testCase.ref.timesS(3).e);
        end
        function divideSolution(testCase)
            s_test = struct(rdivide(testCase.s,testCase.p));
            testCase.verifyEqual(s_test.y,testCase.ref.divideS(1).y);
            testCase.verifyEqual(s_test.e,testCase.ref.divideS(1).e);
            s_test = struct(rdivide(testCase.p,testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.divideS(2).y);
            testCase.verifyEqual(s_test.e,testCase.ref.divideS(2).e);
            s_test = struct(rdivide(testCase.s,testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.divideS(3).y);
            testCase.verifyEqual(s_test.e,testCase.ref.divideS(3).e);
        end
        function expSolution(testCase)
            s_test = struct(exp(testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.expS(1).y);
            testCase.verifyEqual(s_test.e,testCase.ref.expS(1).e);
        end
        function logSolution(testCase)
            s_test = struct(log(testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.logS(1).y);
            testCase.verifyEqual(s_test.e,testCase.ref.logS(1).e);
        end
        function log10Solution(testCase)
            s_test = struct(log10(testCase.s));
            testCase.verifyEqual(s_test.y,testCase.ref.log10S(1).y);
            testCase.verifyEqual(s_test.e,testCase.ref.log10S(1).e);
        end
        function minSolution(testCase)
            testCase.verifyEqual(min(testCase.s),testCase.ref.minS);
        end
        function maxSolution(testCase)
            testCase.verifyEqual(max(testCase.s),testCase.ref.maxS);
        end
        function absSolution(testCase)
            s_test = struct(abs(testCase.s - testCase.p));
            testCase.verifyEqual(s_test.y,testCase.ref.absS.y);
            testCase.verifyEqual(s_test.e,testCase.ref.absS.e);
        end
        function sumSolution(testCase)
            [su, eu] = sum(testCase.s);
            testCase.verifyEqual(su,testCase.ref.sumS.y);
            testCase.verifyEqual(eu,testCase.ref.sumS.e);
        end
        function meanSolution(testCase)
            [su, eu] = mean(testCase.s,'method','mean');
            testCase.verifyEqual(su,testCase.ref.meanS.y(1));
            testCase.verifyEqual(eu,testCase.ref.meanS.e(1));
            [su, eu] = mean(testCase.s,'method','counts');
            testCase.verifyEqual(su,testCase.ref.meanS.y(2));
            testCase.verifyEqual(eu,testCase.ref.meanS.e(2));
            [su, eu] = mean(testCase.s,'method','weight');
            testCase.verifyEqual(su,testCase.ref.meanS.y(3));
            testCase.verifyEqual(eu,testCase.ref.meanS.e(3));
        end
    end
    
    methods (TestMethodSetup)
        function makeIN(testCase)
            rng(314152, 'twister');
            x = linspace(0,10,testCase.N);
            testCase.s = spec1d(x,rand(testCase.N,1),0.1*rand(testCase.N,1));
            testCase.p = rand(1);
        end
        
        function loadReference(testCase)
            testCase.ref  = load('references_Maths.mat');
        end
    end
end

