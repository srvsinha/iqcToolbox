%% Requirements:
% 1. DeltaSubclass.validateSample shall take a Ulft object and sampling-time as 
%     an input and throw an error if the input is a invalid element of the set 
%     represented by DeltaSubclass
% 2. DeltaSublcass.validateSample shall not throw an error if the input is a 
%     valid element of the set represented by DeltaSubclass

%% Requirements:
% 1. DeltaSublcass.sample shall take a sampling-time and output a valid element
%     of the set represented by DeltaSubclass
% 2. If sampling-time is empty, the returned element will be memoryless. If
%     sampling-time is 0, the returned element will be a continuous-time system.
%     If the sampling-time is non-zero, the returned element will be a discrete-
%     time system.

%%
%  Copyright (c) 2021 Massachusetts Institute of Technology 
%  SPDX-License-Identifier: GPL-2.0
%%

%% Test class for sampling methods.
classdef testSamplingMethods < matlab.unittest.TestCase
    
methods (TestMethodSetup)
function seedAndReportRng(testCase)
    seed = floor(posixtime(datetime('now')));
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
function testDeltaSltiSampleAndValidate(testCase)
    % Check that sample and validateSample are consistent with each other
    for i = 1:5
        lft = Ulft.random('req_deltas', {'DeltaSlti'}, 'num_delta', 1);
        del = lft.delta.deltas{1};
        timestep = {-1, 0, 1, []};
        samp = sample(del, timestep{randi([1,4])});
        validateSample(del, samp);
    end
    
    % Check that validateSample rejects known bad samples
    % Out of bounds
    del = DeltaSlti('d');
    bad_samp = toLft(2);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample')
    % Bad dimensions
    bad_samp = toLft(zeros(2));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample')
    % Bad horizon_period
    hp = [1, 1];
    del = matchHorizonPeriod(del, hp);
    bad_samp = toLft(0.5);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample')
    % Not time-invariant
    hp = [1, 1];
    del = matchHorizonPeriod(del, hp);
    a = {[], []};
    b = {zeros(0, 1), zeros(0, 1)};
    c = {zeros(1, 0), zeros(1, 0)};
    d = {.5, -.5};
    bad_samp = Ulft(a, b, c, d, SequenceDelta(), 'horizon_period', hp); 
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample')
    % Not diagonal
    del = DeltaSlti('d', 2, 0, 1);
    bad_samp = toLft(ones(2));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample')
    % No bad Deltas
    bad_samp = sample(del, timestep) + DeltaSlti('b', 2);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample')   
    bad_samp = removeUncertainty(bad_samp, 1) + DeltaDelayZ(2);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSlti:validateSample') 
end

function testDeltaSltvSampleAndValidate(testCase)
    for i = 1:5
        lft = Ulft.random('req_deltas', {'DeltaSltv'}, 'num_delta', 1);
        del = lft.delta.deltas{1};
        timestep = {-1, 0, 1, []};
        samp = sample(del, timestep{randi([1,4])});
        validateSample(del, samp);
    end
    
    % Check that validateSample rejects known bad samples
    % Out of bounds
    del = DeltaSltv('d');
    bad_samp = toLft(2);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')
    % Bad dimensions
    bad_samp = toLft(zeros(2));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')
    % Bad horizon_period
    hp = [1, 1];
    del = matchHorizonPeriod(del, hp);
    bad_samp = toLft(0.5);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')
    % Not diagonal
    del = DeltaSltv('d', 2, 0, 1);
    bad_samp = toLft(ones(2));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')
    % No bad Deltas
    bad_samp = sample(del, timestep) + DeltaSlti('b', 2);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')   
    bad_samp = removeUncertainty(bad_samp, 1) + DeltaDelayZ(2);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')  
end

function testDeltaSltvRateBndSampleAndValidate(testCase)
    for i = 1:5
        lft = Ulft.random('req_deltas', {'DeltaSltvRateBnd'}, 'num_delta', 1);
        del = lft.delta.deltas{1};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% LINES BELOW NEED TO BE REMOVED TO TEST FULL FUNCTIONALITY %%%%%%
        mag_diff = del.upper_bound(1) - del.lower_bound(1);
        del.upper_rate = mag_diff * ones(1, length(del.upper_bound));
        del.lower_rate = -del.upper_rate;
        %%%%%% LINES ABOVE NEED TO BE REMOVED TO TEST FULL FUNCTIONALITY %%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        timestep_choices = {-1, 0, 1, []};
        timestep = timestep_choices{randi([1,4])};
        samp = sample(del, timestep);
        validateSample(del, samp, timestep);
    end    
    % Check that validateSample rejects known bad samples
    % Out of bounds
    del = DeltaSltvRateBnd('d');
    bad_samp = toLft(2);
    timestep = [];
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')
    % Bad dimensions
    bad_samp = toLft(zeros(2));
    timestep = -1;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')
    % Bad horizon_period
    hp = [1, 1];
    del = matchHorizonPeriod(del, hp);
    bad_samp = toLft(0.5);
    timestep = 0;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')
    % Not diagonal
    del = DeltaSltvRateBnd('d', 2, 0, 1);
    bad_samp = toLft(ones(2));
    timestep = 1;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')
    % Out of rate
    dim_outin = 1;
    lbnd = -1;
    ubnd = 1;
    lr = lbnd;
    ur = ubnd;
    hp = [1, 3];
    del = DeltaSltvRateBnd('d', dim_outin, lbnd, ubnd, lr, ur, hp);
    total_time = sum(hp);
    a = repmat({[]}, 1, total_time);
    b = repmat({zeros(0, 1)}, 1, total_time);
    c = repmat({zeros(1, 0)}, 1, total_time);
    % Rate-of-change is unacceptable for timestep = 1 -> 2
    d = {-1, 1, 0, 0};
    bad_samp = Ulft(a, b, c, d, SequenceDelta(), 'horizon_period', hp); 
    timestep = -1;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')
    % Rate-of-change is acceptable for timestep = 1:4, 
    % but transitioning from timesteps 4 -> 5 (or 4 -> 2) is out of rate-bounds.
    d = {-1, -1, 0, 1};
    bad_samp = Ulft(a, b, c, d, SequenceDelta(), 'horizon_period', hp); 
    timestep = -1;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')
    % No bad Deltas
    bad_samp = sample(del, timestep) + DeltaSlti('b');
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')   
    bad_samp = removeUncertainty(bad_samp, 1) + DeltaDelayZ();
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRateBnd:validateSample')   
end

function testDeltaSltvRepeatedSampleAndValidate(testCase)
    for i = 1:5
        num_dels = 3;
        names = {'a', 'b', 'c'};
        horizon_period = [randi([0, 10]), randi([1, 10])];
        total_time = sum(horizon_period);
        dim_dels = randi([1, 5], num_dels, 1) * ones(1, total_time);
        upper_bound = rand(num_dels, total_time);
        lower_bound = -upper_bound - 1;
        region = cellfun(@(lb, ub) [lb, ub],...
                         num2cell(lower_bound, 1),...
                         num2cell(upper_bound, 1),...
                         'UniformOutput', false);
        del = DeltaSltvRepeated(names, dim_dels, 'box', region, horizon_period);
        timestep = {-1, 0, 1, []};
        samp = sample(del, timestep{randi([1,4])});
        validateSample(del, samp);
    end
    
    % Check that validateSample rejects known bad samples
    % Out of bounds
    del = DeltaSltvRepeated({'d1', 'd2'}, [1; 1]);
    bad_samp = toLft(diag([0, 2]));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')
    % Bad dimensions
    bad_samp = toLft(zeros(3));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltvRepeated:validateSample')
    % Bad horizon_period
    hp = [1, 1];
    del = matchHorizonPeriod(del, hp);
    bad_samp = toLft(0.5 * eye(2));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltvRepeated:validateSample')
    % Not diagonal
    del = DeltaSltvRepeated({'d1', 'd2'}, [1; 2]);
    bad_samp = toLft(blkdiag(1, ones(2)));
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')
    % No bad Deltas
    bad_samp = sample(del, timestep) + DeltaDelayZ(3);
    verifyError(testCase,...
                @() validateSample(del, bad_samp),...
                'DeltaSltv:validateSample')   
    bad_samp = removeUncertainty(bad_samp, 1) + DeltaSlti('b', 3);
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaSltvRepeated:validateSample')
end

function testDeltaDltiSampleAndValidate(testCase)
    for i = 1:5
        lft = Ulft.random('req_deltas', {'DeltaDlti'},...
                          'num_delta', 1,...
                          'horizon_period', [0,1]);
        del = lft.delta.deltas{1};
        timestep_choices = {-1, 0, 1, []};
        timestep = timestep_choices{randi([1,4])};
        % What if timestep indicates no dynamics?
        if isempty(timestep)
            verifyError(testCase,...
                        @() sample(del, timestep),...
                        'DeltaDlti:sample')
        else
            samp = sample(del, timestep);
            % What if timestep indicates discrete-time?
            if timestep 
                verifyClass(testCase, samp.delta.deltas{1}, 'DeltaDelayZ')
            else
                verifyClass(testCase, samp.delta.deltas{1}, 'DeltaIntegrator')
            end
            validateSample(del, samp, timestep); 
        end
    end
    
    % Sampling time-invariant with non-default horizon periods
    del = DeltaDlti('d');
    hp = [2, 3];
    del = matchHorizonPeriod(del, hp);
    samp = sample(del, 1);
    verifyEqual(testCase, samp.horizon_period, del.horizon_period)
    verifyError(testCase, @() sample(del, false), 'DeltaDlti:sample')
    
    % Check that validateSample rejects known bad samples
    % Out of bounds
    del = DeltaDlti('d');
    upper_bound = 2;
    dim_state = 3;
    bad_samp = rss(dim_state);
    bad_samp.a = bad_samp.a - 0.1 * eye(dim_state);
    bad_samp = bad_samp * upper_bound / norm(bad_samp, 'inf');
    bad_samp = toLft(bad_samp);
    discrete = false;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, discrete),...
                'DeltaDlti:validateSample')
    % Bad dimensions
    bad_samp = toLft(zeros(2));
    discrete = [];
    verifyError(testCase,...
                @() validateSample(del, bad_samp, discrete),...
                'DeltaDlti:validateSample') 
    % Bad horizon_period
    hp = [1, 1];
    bad_samp = drss(dim_state);
    bad_samp.a = bad_samp.a * 0.9;
    bad_samp = 0.99 * bad_samp / norm(bad_samp, 'inf');
    bad_samp = toLft(bad_samp).matchHorizonPeriod(hp);
    discrete = true;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, discrete),...
                'DeltaDlti:validateSample') % Convert to time-invariant
    % Not time-invariant
    bad_samp.a{1} = zeros(dim_state);
    bad_samp.a{2} = 0.25 * ones(dim_state);
    verifyError(testCase,...
                @() validateSample(del, bad_samp, discrete),...
                'DeltaDlti:validateSample')
    % No bad Deltas
    timestep = timestep_choices{randi([1,3])};
    bad_samp = sample(del, timestep) + DeltaSlti('a');
    verifyError(testCase,...
                @() validateSample(del, bad_samp, discrete),...
                'DeltaDlti:validateSample')    
end

function testDeltaBoundedSampleAndValidate(testCase)
    dynamic_system = false;
    while ~dynamic_system
        timestep_choices = {-1, 0, 1, []};
        timestep = timestep_choices{randi([1,4])};
        if timestep == 0
            lft = Ulft.random('req_deltas', {'DeltaBounded'},...
                              'num_deltas', 1,...
                              'horizon_period', [0,1]);
        else
            lft = Ulft.random('req_deltas', {'DeltaBounded'}, 'num_deltas', 1);
        end
        del = lft.delta.deltas{1};
        samp = sample(del, timestep);
        validateSample(del, samp);
        dynamic_system = ~isempty(samp.delta.deltas);
    end
    
    % Check that sample won't make a bad horizon_period continuous-time LFT
    del = DeltaBounded('d').matchHorizonPeriod([1, 1]);
    timestep = 0;
    verifyError(testCase, @() del.sample(timestep), 'DeltaBounded:sample')
    
    % Check that validateSample rejects known bad samples
    % Out of bounds
    del = DeltaBounded('d');
    out_of_bound = 2;
    dim_state = 3;
    bad_samp = drss(dim_state);
    bad_samp.a = bad_samp.a - 0.1 * eye(dim_state);
    bad_samp = bad_samp * out_of_bound / norm(bad_samp, 'inf');
    bad_samp = toLft(bad_samp);
    timestep = 1;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaBounded:validateSample')
    bad_samp = toLft(out_of_bound);
    timestep = [];
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaBounded:validateSample')
    
    % Bad dimensions
    bad_samp = toLft(zeros(2));
    timestep = [];
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaBounded:validateSample') 
    % Bad horizon_period
    hp = [1, 1];
    bad_samp = drss(dim_state);
    bad_samp.a = bad_samp.a * 0.9;
    bad_samp = 0.99 * bad_samp / norm(bad_samp, 'inf');
    bad_samp = toLft(bad_samp).matchHorizonPeriod(hp);
    timestep = -1;
    verifyError(testCase,...
                @() validateSample(del, bad_samp, timestep),...
                'DeltaBounded:validateSample') 
end

function testDeltaSampleAndValidate(testCase)
    del = DeltaSltvRateBnd('d');
    [~, ~, ~, del] = recastMatricesAndDelta(del);
    verifyError(testCase, @() sample(del), 'Delta:sample')
    verifyError(testCase, @() validateSample(del), 'Delta:validateSample')
end

function testDeltaDelayZSampleAndValidate(testCase)
    del = DeltaDelayZ();
    verifyError(testCase, @() sample(del), 'DeltaDelayZ:sample')
    verifyError(testCase, @() validateSample(del), 'DeltaDelayZ:validateSample')
end

function testDeltaIntegratorSampleAndValidate(testCase)
    d = DeltaIntegrator();
    verifyError(testCase, @() sample(d), 'DeltaIntegrator:sample')
    verifyError(testCase,...
                @() validateSample(d),...
                'DeltaIntegrator:validateSample')
end
end
end

%%  CHANGELOG
% Sep. 28, 2021 (v0.6.0): Added after v0.5.0 - Micah Fry (micah.fry@ll.mit.edu)