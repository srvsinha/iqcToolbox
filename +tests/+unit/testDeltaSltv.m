%% Requirements:
%  1. DeltaSltv shall be defined by it's name, in/out dimensions, 
%      lower/upper bounds, and horizon_period.
%  2. Upon construction, and when queried by user, it shall display the
%      information described in (1).
%
%  3. If dimension and/or bound information is not provided by the user, by
%      default the object shall be 1 x 1, with bounds [-1, 1], and 
%      horizon_period [0, 1].
%
%  4. If the user provides no name, or the name is not a string, DeltaSltv 
%      shall throw an exception
%  5. If the user provides an in/out dimension that is not a natural number
%      DeltaSltv shall throw an exception
%
%  6. The in/out dimensions of DeltaSltv shall be equal.
%  7. If the user provides a lower bound which is greater than the upper
%      bound, DeltaSltv shall throw an exception
%  8. If the user provides +/- inf or Nan for either of the bounds, 
%      DeltaSltv shall throw an exception
%
%  9. DeltaSltv shall ensure that it's properties are consistent with its
%      current horizon_period property
%  10.DeltaSltv shall be capable of changing it's properties to match a
%      newly input horizon_period, as long as the new horizon_period is
%      consistent with the prior horizon_period
%
%  11.DeltaSltv shall return the necessary mappings and DeltaSltv
%      object to be used for normalizing an LFT with a DeltaSltv
%      uncertainty.
%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for DeltaSltv.
classdef testDeltaSltv < matlab.unittest.TestCase

methods (TestMethodSetup)
function seedAndReportRng(testCase)
    seed = floor(posixtime(datetime('now')));
    rng('default');
    rng(seed);
    diagnose_str = ...
        sprintf(['Random inputs may be regenerated by calling: \n',...
                 '>> rng(%10d) \n',...
                 'before running the remainder of the test''s body'],...
                seed);
    testCase.onFailure(@() fprintf(diagnose_str));
end    
end
    
methods (Test)

function testFullConstructor(testCase)
    name = 'test';
    dim_outin = 2;
    lower_bound = -5.0;
    upper_bound = 10.0;
    horizon_period = [0, 1];
    delta_sltv = DeltaSltv(name,...
                           dim_outin,...
                           lower_bound,...
                           upper_bound,...
                           horizon_period);
    verifyEqual(testCase, delta_sltv.name, name)
    verifyEqual(testCase, delta_sltv.dim_in, dim_outin)
    verifyEqual(testCase, delta_sltv.dim_out, dim_outin)
    verifyEqual(testCase, delta_sltv.dim_in, delta_sltv.dim_out)
    verifyEqual(testCase, delta_sltv.lower_bound, lower_bound)
    verifyEqual(testCase, delta_sltv.upper_bound, upper_bound)
    verifyEqual(testCase, delta_sltv.horizon_period, horizon_period)
end

function testFullConstructorDifferentHorizonPeriod(testCase)
    name = 'test';
    dim_outin = 2;
    lower_bound = -5.0;
    upper_bound = [10.0, 7, 8, 9, 9];
    horizon_period = [3, 2];
    total_time = sum(horizon_period);
    delta_sltv = DeltaSltv(name,...
                           dim_outin,...
                           lower_bound,...
                           upper_bound,...
                           horizon_period);
    verifyEqual(testCase, delta_sltv.name, name)
    verifyEqual(testCase,...
                delta_sltv.dim_in,...
                dim_outin * ones(1, total_time))
    verifyEqual(testCase,...
                delta_sltv.dim_out,...
                dim_outin * ones(1, total_time))
    verifyEqual(testCase, delta_sltv.dim_in, delta_sltv.dim_out)
    verifyEqual(testCase,...
                delta_sltv.lower_bound,...
                lower_bound * ones(1, total_time))
    verifyEqual(testCase, delta_sltv.upper_bound, upper_bound)
    verifyEqual(testCase, delta_sltv.horizon_period, horizon_period)
end

function testOneArgConstructor(testCase)
    name = 'test';
    delta_sltv = DeltaSltv(name);
    verifyEqual(testCase, delta_sltv.name, name)
    verifyEqual(testCase, delta_sltv.dim_in, 1)
    verifyEqual(testCase, delta_sltv.dim_out, 1)
    verifyEqual(testCase, delta_sltv.dim_in, delta_sltv.dim_out)
    verifyEqual(testCase, delta_sltv.lower_bound, -1.0)
    verifyEqual(testCase, delta_sltv.upper_bound, 1.0)
    verifyEqual(testCase, delta_sltv.horizon_period, [0, 1])
end

function testTwoArgConstructor(testCase)
    name = 'test';
    dim_outin = 7;
    delta_sltv = DeltaSltv(name, dim_outin);
    verifyEqual(testCase, delta_sltv.name, name)
    verifyEqual(testCase, delta_sltv.dim_in, dim_outin)
    verifyEqual(testCase, delta_sltv.dim_out, dim_outin)
    verifyEqual(testCase, delta_sltv.dim_in, delta_sltv.dim_out)
    verifyEqual(testCase, delta_sltv.lower_bound, -1.0)
    verifyEqual(testCase, delta_sltv.upper_bound, 1.0)
    verifyEqual(testCase, delta_sltv.horizon_period, [0, 1])
end

function testFourArgConstructor(testCase)
    name = 'test';
    dim_outin = 7;
    upper_bound = 4.5;
    lower_bound = -0.1;
    delta_sltv = DeltaSltv(name, dim_outin, lower_bound, upper_bound);
    verifyEqual(testCase, delta_sltv.name, name)
    verifyEqual(testCase, delta_sltv.dim_in, dim_outin)
    verifyEqual(testCase, delta_sltv.dim_out, dim_outin)
    verifyEqual(testCase, delta_sltv.dim_in, delta_sltv.dim_out)
    verifyEqual(testCase, delta_sltv.lower_bound, lower_bound)
    verifyEqual(testCase, delta_sltv.upper_bound, upper_bound)
    verifyEqual(testCase, delta_sltv.horizon_period, [0, 1])
end

function testDeltaToMultiplier(testCase)
    name = 'test';
    delta_sltv = DeltaSltv(name);
    default_mult = deltaToMultiplier(delta_sltv);
    verifyEqual(testCase, default_mult.quad_time_varying, true)
    
    quad_time_varying = false;
    mult = deltaToMultiplier(delta_sltv,...
                             'quad_time_varying',...
                             quad_time_varying);
    verifyEqual(testCase, mult.quad_time_varying, quad_time_varying)
end

function testFailedName(testCase)
    verifyError(testCase, @() DeltaSltv(), ?MException)
    verifyError(testCase, @() DeltaSltv(1), ?MException)
end

function testFailedDimension(testCase)
    verifyError(testCase, @() DeltaSltv('test', -2), ?MException)
    verifyError(testCase, @() DeltaSltv('test', 2.2), ?MException)
end

function testFailedBounds(testCase)
    hp = [2, 2];
    total_time = sum(hp);
    lb = -1.2 * ones(1,total_time);
    ub = linspace(4, 5, total_time);
    lb(total_time) = ub(total_time) + 1;

    verifyError(testCase, @() DeltaSltv('test', 1, lb, ub, hp), ?MException)

    lb(total_time) = -inf;
    verifyError(testCase, @() DeltaSltv('test', 1, lb, ub, hp), ?MException)

    lb(total_time) = ub(total_time) - 1;
    ub(1) = nan;
    verifyError(testCase, @() DeltaSltv('test', 1, lb, ub, hp), ?MException)
end

function testHorizonPeriod(testCase)
   name = 'test';
   delta_sltv = DeltaSltv(name);
   assertEqual(testCase, delta_sltv.horizon_period, [0, 1])
   
   % Checking horizon_period and making sure it fits for all properties
   horizon_period2 = [4, 3];
   total_time2 = sum(horizon_period2);
   lb = -1;
   ub = linspace(2, 3, total_time2);
   dim_outin2 = ones(1, total_time2);
   delta_sltv = DeltaSltv(name, dim_outin2, lb, ub, horizon_period2);
   delta_sltv = matchHorizonPeriod(delta_sltv);
   assertEqual(testCase, delta_sltv.horizon_period, horizon_period2)
   assertEqual(testCase, delta_sltv.lower_bound, lb * ones(1, total_time2))
   assertEqual(testCase, delta_sltv.upper_bound, ub)
   assertEqual(testCase, delta_sltv.dim_in, dim_outin2)
   assertEqual(testCase, delta_sltv.dim_out, dim_outin2)

   % Resetting horizon_period and making sure it fits for all properties
   horizon_period3 = [total_time2, horizon_period2(2) * 2];
   total_time3 = sum(horizon_period3);
   dim_outin3 = ones(1, total_time3);
   delta_sltv = matchHorizonPeriod(delta_sltv, horizon_period3);
   verifyEqual(testCase, delta_sltv.horizon_period, horizon_period3)
   verifyEqual(testCase, delta_sltv.lower_bound, lb * ones(1, total_time3));
   verifyEqual(testCase,...
               delta_sltv.upper_bound,...
               [ub,...
                ub((1:(horizon_period3(2) / 2)) + horizon_period2(1)),...
                ub((1:(horizon_period3(2) / 2)) + horizon_period2(1))])
   verifyEqual(testCase, delta_sltv.dim_in,  dim_outin3)
   verifyEqual(testCase, delta_sltv.dim_out, dim_outin3)
end

function testFailedHorizonPeriod(testCase)
   name = 'test';
   delta_sltv = DeltaSltv(name);
   assertEqual(testCase, delta_sltv.horizon_period, [0, 1])
   
   % Resetting horizon_period and making sure it fits for all properties
   horizon_period2 = [4, 7];
   delta_sltv.horizon_period = horizon_period2;
   delta_sltv = matchHorizonPeriod(delta_sltv);
   assertEqual(testCase, delta_sltv.horizon_period, horizon_period2)
   
   % Resetting horizon_period and incorrectly trying to force a fit with
   % other properties
   horizon_period3 = [5, 3];
   horizon_period3 = commonHorizonPeriod([horizon_period2; horizon_period3]);
   delta_sltv.horizon_period = horizon_period3;
   verifyError(testCase, @() matchHorizonPeriod(delta_sltv), ?MException)
end

function testNormalization(testCase)
    % 10 randomly generated LFTs
    for i = 1:10
        dim_outin = randi([1, 10]);
        hp = [randi([0, 10]), randi([1, 10])];
        total_time = sum(hp);
        upper_bound = 10 * rand(1, total_time);
        lower_bound = upper_bound - 20 * rand(1, total_time);
        width = upper_bound - lower_bound;
        del = DeltaSltv('test', dim_outin, lower_bound, upper_bound, hp);
        lft = Ulft.random('num_deltas', 1,...
                          'req_deltas', {del},...
                          'horizon_period', hp);
        lft_n = normalizeLft(lft);
        
        % Check correctness
        % delta
        expected_del = DeltaSltv('test', dim_outin, -1, 1, hp);
        verifyEqual(testCase, lft_n.delta.deltas{1}, expected_del)
        % a, b, c, d (via sampling)
        for j = 1 :3
            kern = rand;
            for k = 1:total_time
                del_samp = (lower_bound(k) + width(k) * kern) * eye(dim_outin);
                del_n_samp = (-1 + 2 * kern) * eye(dim_outin);

                lft_samp = lft.d{k} + lft.c{k} * del_samp * ...
                           inv(eye(dim_outin) - lft.a{k} * del_samp) * ...
                           lft.b{k};
                lft_n_samp = lft_n.d{k} + lft_n.c{k} * del_n_samp * ...
                             inv(eye(dim_outin) - lft_n.a{k} * del_n_samp) * ...
                             lft_n.b{k};
                verifyLessThan(testCase,...
                               max(abs(lft_samp - lft_n_samp)),...
                               1e-5)
            end
        end
    end
end

end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0)
% Aug. 26, 2021 (v.0.5.0): Initial release - Micah Fry (micah.fry@ll.mit.edu)