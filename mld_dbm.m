% de Boyer Montegut et al. [2004] mld algorithm
% NAP-deuce

% input in-situ temperature, salinity, and cast depths/pressures

% returns mixed layer depth, isothermal layer depth, and sorted profile of
% potential density
function [mld_out,ild_out,sig_theta] = mld_dbm(temp_cast,salin_cast,depth_cast,theta_yes)

% sort cast by increasing depth - depth and pressure treated as equivalent
% from this point on though this could be later modified using sw_pres if a
% latitude argument is added above
    [dep_sort,ix] = sort(depth_cast);
    temp = temp_cast(ix);
    salin = salin_cast(ix);

% compute potential temperature and potential density
    if theta_yes
        theta = temp_cast(ix);
    else
        theta = sw_ptmp(salin,temp,dep_sort,0);
    end
    sig_theta = sw_dens0(salin,theta);

% find reference depth
    ref_ind = find(dep_sort > 9,1,'first'); % unevenly spaced data
    % ref_ind = find(dep_sort == 10,1,'first');

% reference depth is too deep - no shallow depths in cast
    if dep_sort(ref_ind) > 25
        mld_out = NaN;
        ild_out = NaN;
        return
    end

% density increase from 10m reference using delta-theta of 0.2 deg C
    ref_sig_theta = sw_dens0(salin(ref_ind),theta(ref_ind)-0.2);

% search for MLD
    if sum(~isnan(sig_theta)) > 1 % not a one-point or all-NaN cast
    % find mixed layer depth
        not_found = 1;
        start_ind = ref_ind;
        iter_count = 1;
        while not_found            
    % begin search at reference (10 m) index
    % finds next point below reference depth that exceeds criterion

            ml_ind = find(sig_theta(start_ind:end) > ref_sig_theta,1,'first');
            if length(sig_theta) >= ml_ind+start_ind % ml_ind is within the interior of the cast
                if sig_theta(ml_ind+start_ind) > ref_sig_theta % next point also meets criterion, therefore likely not a spike
                    not_found = 0;
                    ml_ind = ml_ind + start_ind - 1; % final index
                end
            else % last point in cast
                not_found = 0;
                ml_ind = ml_ind + start_ind - 1;
            end
            % if a spike, start search again at first point after spike
            start_ind = ml_ind+start_ind;
            iter_count = iter_count+1;
            % break loop if cast is all spikes/no MLD found
            if iter_count > length(sig_theta)
                break
            end
        end
        % if an MLD is found, interpolate to find depth at which density = ref_sig_theta
        if ~isempty(ml_ind)
            if ~not_found && sum(isnan(sig_theta(ml_ind-1:ml_ind))) == 0
                mld_out = interp1(sig_theta(ml_ind-1:ml_ind),dep_sort(ml_ind-1:ml_ind),ref_sig_theta);
            else
                mld_out = NaN;
            end
        else
            mld_out = NaN;
        end

        % find isothermal layer depth
        not_found = 1;
        start_ind = ref_ind;
        iter_count = 1;
        while not_found
            il_ind = find(theta(start_ind:end) < theta(ref_ind)-0.2,1,'first');
            if length(theta) >= il_ind+start_ind
                if theta(il_ind+start_ind) < theta(ref_ind)-0.2 % next point also meets criterion, therefore likely not a spike
                    not_found = 0;
                    il_ind = il_ind + start_ind - 1; % final index
                end
            else
                il_ind = il_ind + start_ind - 1;
                not_found = 0;
            end
            start_ind = il_ind+start_ind;
            iter_count = iter_count+1;
            if iter_count > length(theta)
                break
            end
        end
        if ~isempty(il_ind)
            if ~not_found && sum(isnan(theta(il_ind-1:il_ind))) == 0
                ild_out = interp1(theta(il_ind-1:il_ind),dep_sort(il_ind-1:il_ind),theta(ref_ind)-0.2);
            else
                ild_out = NaN;
            end
        else
            ild_out = NaN;
        end
    else
        mld_out = NaN;
        ild_out = NaN;
    end