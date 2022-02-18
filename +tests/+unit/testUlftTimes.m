%% Requirements:
%  1. Ulft.times shall perform the .* multiplication operation 
%     between two LFTs.
%     1.1. If both LFTs have conformable dimensions, the .* operation shall 
%          replicate the mtimes (*) operation.
%     1.2. If only one of the LFTs has a size of 1 x 1, then the .* operation
%          shall multiply each element of the non-scalar LFT by the scalar LFT.
%
%  2. If the user provides LFTs whose dimensions do not satisfy the requirements
%     to 1.1 or 1.2, then Ulft.times shall throw an error.
%
%  3. If the user provides two LFTs that have conformable dimensions, Ulft.times
%     shall issue a warning that it is replicating the behavior of mtimes.
%
%  4. Ulft.times shall be capable of taking one input argument that is not 
%     a Ulft object, as long as the input argument can be converted to a 
%     Ulft. If the non-Ulft input is not convertible, an error shall be 
%     thrown. Objects convertible to Ulfts are:
%       - doubles
%       - Delta objects
%       - ss objects
%
%  5. If the Ulft.times operands have different horizon_periods, 
%     Ulft.times shall ensure that the output Ulft has a resulting 
%     horizon_period that is consistent with the horizon_periods of both 
%     operands.
%

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for Ulft.
classdef testUlftTimes < matlab.unittest.TestCase
    
methods (TestMethodSetup)
    function seedAndReportRng(testCase)
        seed = floor(posixtime(datetime('now')));
        rng('default');
        rng(seed, 'twister');
        diagnose_str = ...
            sprintf(['Random inputs may be regenerated by calling: \n',...
                     '>> rng(%10d) \n',...
                     'before running the remainder of the test''s body'],...
                    seed);
        testCase.onFailure(@() fprintf(diagnose_str));
    end    
end
    
methods (Test)

function testTimesCorrectness(testCase)
    % Conformable dimensions
    ss_left = drss;
    ss_right = drss(randi([1, 10]), size(ss_left, 2), randi([1, 10]));
    ss_out = ss_left * ss_right;
    lft_left = toLft(ss_left);
    lft_right = toLft(ss_right);
    % Make sure warning is thrown when doing typical mtimes
    warning_state = warning;
    warning('on', 'Ulft:times')
    testCase.verifyWarning(@()  lft_left .* lft_right, 'Ulft:times');
    warning(warning_state)
    lft_out = lft_left .* lft_right;
    % Check correctness of .*
    testCase.verifyEqual(ss_out.a, lft_out.a{1});
    testCase.verifyEqual(ss_out.b, lft_out.b{1});
    testCase.verifyEqual(ss_out.c, lft_out.c{1});
    testCase.verifyEqual(ss_out.d, lft_out.d{1});
    
    % When lft_left is scalar (need to compare via Linf norm)
    scalar_left = drss(1, 1, 1);
    scalar_left.a = scalar_left.a * 0.9;
    ss_right = drss(randi([1, 10]), randi([2, 10]), randi([2, 10]));
    ss_right.a = ss_right.a * 0.9;
    ss_out = scalar_left * ss_right;
    lft_left = toLft(scalar_left);
    lft_right = toLft(ss_right);
    lft_out = lft_left .* lft_right;
    sys_diff = norm(ss_out - lftToSs(lft_out), 'inf');
    testCase.verifyLessThan(sys_diff, 1e-8);
    
    % When lft_right is scalar (need to compare via Linf norm)
    ss_left = drss(randi([1, 10]), randi([2, 10]), randi([2, 10]));
    ss_left.a = ss_left.a * 0.9;
    scalar_right = drss(1, 1, 1);
    scalar_right.a = scalar_right.a * 0.9;
    ss_out = ss_left * scalar_right;
    lft_left = toLft(ss_left);
    lft_right = toLft(scalar_right);
    lft_out = lft_left .* lft_right;
    sys_diff = norm(ss_out - lftToSs(lft_out), 'inf');
    testCase.verifyLessThan(sys_diff, 1e-8);
end

function testTimesNonLftRandom(testCase)
    % REQUIREMENT 4
    % Ensure failure of multiplication of non-lft ojects
        
    % Define lfts
    % Discrete time lft
    z = DeltaDelayZ();
    a = {0};
    b = {1};
    c = {1};
    d = {0};
    lft1 = Ulft(a, b, c, d, z);
    
    % - LTI
    lfti_l = Ulft.random('req_deltas', {DeltaIntegrator()});
    lfti_r = Ulft.random('req_deltas', {DeltaIntegrator()}, ...
                         'dim_out', 1);
    
    lfti_l_disc = Ulft.random('req_deltas', {DeltaDelayZ()});
    lfti_r_disc = Ulft.random('req_deltas', {DeltaDelayZ()}, ...
                              'dim_out', 1);
    
    lfti_l_size = size(lfti_l, 2);
    lfti_l_disc_size = size(lfti_l_disc, 2);
    % Convertible non-lft objects
    lft_double = ones(lfti_l_size, 1);
    lft_delta = DeltaSlti('a');
    lft_ss_cont = ss(1, 1, ones(lfti_l_size, 1), ...
        ones(lfti_l_size, 1));
    lft_ss_disc = ss(1, 1, ones(lfti_l_disc_size, 1), ...
        ones(lfti_l_disc_size, 1), -1);
    
    obj = 'non-convertible';
    
    % non-lft object is right
    lft_a1 = times(lfti_l, lft_double);
    lft_a3 = times(lft1, lft_delta);
    lft_a4 = times(lfti_l, lft_ss_cont);
    lft_a5 = times(lfti_l_disc, lft_ss_disc);
    
    verifyClass(testCase, lft_a1, 'Ulft') 
    verifyClass(testCase, lft_a3, 'Ulft') 
    verifyClass(testCase, lft_a4, 'Ulft') 
    verifyClass(testCase, lft_a5, 'Ulft') 
    
    % non-lft object is left
    lft_b1 = times(lft_double, lfti_r);
    lft_b3 = times(lft_delta, lft1);
    lft_b4 = times(lft_ss_disc, lfti_r_disc);
    lft_b5 = times(lft_ss_cont, lfti_r);
    
    verifyClass(testCase, lft_b1, 'Ulft') 
    verifyClass(testCase, lft_b3, 'Ulft') 
    verifyClass(testCase, lft_b4, 'Ulft') 
    verifyClass(testCase, lft_b5, 'Ulft') 
    
    % non-convertible non-lft object - right
    verifyError(testCase, @() times(lfti_l, obj), ?MException)
    % non-convertible non-lft object - left
    verifyError(testCase, @() times(obj, lfti_r), ?MException)    
end

function testTimesCommonHp(testCase)
    % REQUIREMENT 5
    z = DeltaDelayZ();
    a = {0};
    b = {1};
    c = {1};
    d = {0};
    a1 = {1, 1, 1, 1};
    b1 = {1, 1, 1, 1};
    c1 = {1, 1, 1, 1};
    d1 = {0, 0, 0, 0};
    lft1 = Ulft(a, b, c, d, z, 'horizon_period', [0 1]);
    lft2 = Ulft(a1, b1, c1, d1, z, 'horizon_period', [2 2]);
    correct_hp = [2 2];
    lftm = times(lft1, lft2);
    verifyEqual(testCase, lftm.horizon_period, correct_hp)
end
end
end

%%  CHANGELOG
% Nov. 29, 2021: Added after v0.6.0 - Micah Fry (micah.fry@ll.mit.edu)