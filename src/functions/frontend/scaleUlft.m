function this_lft = scaleUlft(this_lft, eps, delta)
%% SCALEULFT method for scaling uncertainty bounds
%
%     lft_out = removeUncertainty(this_lft, 0.5)
%     lft_out = removeUncertainty(this_lft, 0.5, 'DX')
%     lft_out = removeUncertainty(this_lft, [0.5,0.2], {'DX', 'DY'})
%     lft_out = removeUncertainty(this_lft, 0.5, {'DX', 'DY'})
%     lft_out = removeUncertainty(this_lft, [0.5,0.2,0.1], [1,5,6])
%
%     Variables:
%     ---------
%       Input:
%         this_lft : Ulft object :: the lft whose uncertainty bounds are to be scaled
%         delta : string/cell of strings/array of natural numbers :: names/indices of the 
%                 uncertainties that are to be scaled
%         eps : scalar/array of real numbers :: scales for each/all
%               uncertainties in delta (or all uncertainties in this_lft)
%       Output:
%         this_lft : Ulft object :: the lft whose uncertainty bounds are scaled
%
%     See also Ulft, removeDisturbance, removePerformance.


%% Check and process inputs

switch nargin        
    case 1
        error('Must provide the scaling factor')     
    case 2
        delta = {};                
end

validateattributes(eps, {'numeric'}, {'nonempty', 'nonnegative',...
                                      'real', 'finite', 'nonnan'})

if isempty(delta)   
    % scale all uncertainties with the same factor
    
    validateattributes(eps, {'numeric'}, {'scalar'})  
    delta_ind = 1:length(this_lft.delta.deltas);
    delay_ind = find(strcmp(this_lft.delta.types,'DeltaDelayZ'));
    delta_ind(delay_ind) = [];
    
    if eps == 0
        % nominal lft (by removing all uncertainties)
        this_lft = removeUncertainty(this_lft,delta_ind);
        return;
    end
    
else
    % scale the uncertainties whose names are provides   
    
    validateattributes(delta, {'numeric', 'char', 'cell'}, {'nonempty'})

    isChar = @(del) isa(del, 'char');
    % Convert all input types to be indices
    if isChar(delta)
        delta_ind = find(strcmp(delta, this_lft.delta.names), 1);
        if isempty(delta_ind)
            warning('Ulft:removeUncertainty',...
                    ['The named uncertainty, ',...
                     delta,...
                     ', does not appear in the lft'])
        end
    elseif isa(delta, 'cell') && all(all(cellfun(isChar, delta)))
        delta_ind = [];
        for i = 1:length(delta)
            delta_i = delta{i};
            delta_ind = [delta_ind, find(strcmp(delta_i, this_lft.delta.names))];
            if isempty(find(strcmp(delta_i, this_lft.delta.names), 1))
                warning('Ulft:removeUncertainty',...
                        ['The named uncertainty, ',...
                         delta_i,...
                         ', does not appear in the lft'])
            end
        end
    elseif isnumeric(delta)
        validateattributes(delta, {'numeric'}, {'positive',...
                                              'integer',...
                                              'finite',...
                                              'nonnan'})
            bad_inds = delta(delta > length(this_lft.delta.names));
            if ~isempty(bad_inds)
                warning('Ulft:removeUncertainty',...
                        ['The numbered uncertainty, ',...
                         num2str(bad_inds),...
                         ', does not appear in the lft'])
            end
        delta_ind  = delta;
    end
       
    nunc = numel(delta_ind); % # of uncertainties to scale
    if numel(eps) ~= 1 && numel(eps) ~= nunc
        error(['sacling factor must be a scalar or a vector with ', num2str(nunc), ' elements']);
    end   
end


if numel(eps) == 1
    eps = eps*ones(1,numel(delta_ind));
end


%% Scaling the uncertainties
epsind = 1;
for i = delta_ind
    
    scale = eps(epsind);
    
    if scale == 0
        error('Use removeUncertainty function to remove specific uncertainties')
    end
    
        del = this_lft.delta.deltas{1,i};
        prop = properties(del);

        if ~isempty(find(strcmp(prop,'upper_bound')))
            this_lft.delta.deltas{1,i}.upper_bound = scale*this_lft.delta.deltas{1,i}.upper_bound;
        end

        if ~isempty(find(strcmp(prop,'lower_bound')))
            this_lft.delta.deltas{1,i}.lower_bound = scale*this_lft.delta.deltas{1,i}.lower_bound;
        end

        if ~isempty(find(strcmp(prop,'upper_rate')))
            this_lft.delta.deltas{1,i}.upper_rate = scale*this_lft.delta.deltas{1,i}.upper_rate;
        end

        if ~isempty(find(strcmp(prop,'lower_rate')))
            this_lft.delta.deltas{1,i}.lower_rate = scale*this_lft.delta.deltas{1,i}.lower_rate;
        end
    
    epsind = epsind + 1;
end


  
end

%%  CHANGELOG
% Oct. 28, 2022 - Sourav Sinha (srvsinha@vt.edu)