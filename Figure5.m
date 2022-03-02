%% load data
load('contrast_thresholds.mat', 'active', 'passive', 'active_fix', 'passive_fix');

markersize = 6;
capsize = 8;

%% get normalized sensitivity values from thresholds
sUse = 1; % only use first spatial bin (0-30arcmin)
tUse = 1:8;

for ddi = 1:2
    if ddi == 1
        dt = active;
        dt_fix = active_fix;
    else
        dt = passive;
        dt_fix = passive_fix;
    end
    
    % compute sensitivity and normalized sensitivity
    dt.sensitivityBS = 1 ./ dt.threshBS;
    dt_fix.sensitivityBS = 1 ./ dt_fix.threshBS;
    
    dt_fix.sensitivity_mean = nanmean(dt_fix.sensitivityBS, 4);
    
    dt.nsensitivityBS = dt.sensitivityBS ./ dt_fix.sensitivity_mean;
    
    dt.nsensitivity_mean = nanmean(dt.nsensitivityBS, 4);
    dt.nsensitivity_se = nanstd(dt.sensitivityBS, [], 4);
    
    
    if ddi == 1
        active = dt;
        active_fix = dt_fix;
    else
        passive = dt;
        stimulated_fix = dt_fix;
    end
end

%% two-sided non-parametric bootstrap test at each space and time bin
for subjIdx = 1:2
    for sv = sUse
        for tv = tUse
            dt1 = squeeze(active.nsensitivityBS(sv, tv, subjIdx, :)); 
            dt2 = squeeze(passive.nsensitivityBS(sv, tv, subjIdx, :));
            
            p = mean(dt1 < dt2);
            if p > 0.5
                p = 1 - p;
            end
            pvals(sv, tv, subjIdx) = 2 * p; % times 2 for two-sided test
        end
    end
end

%% plot
cols = [...
    0    0.4470    0.7410; ...
    0.9290    0.6940    0.1250];
cols_fill = [cols(1, :); ones(1, 3)];
symbs = 'ds';

tm = (active.TmpVec(1:end-1) + active.TmpVec(2:end))/2;
tshifts = [-3, 3];

ttls = {'Subject A', 'Subject B'};

figure();
for subjIdx = 1:2
    subplot(1, 2, subjIdx); hold on;
    plot([-200, 200], [1, 1], 'k--');
    plot([0, 0], [0, 3], 'k--');
    
    for ddi = 1:2
        if ddi == 1
            dt = active;
        else
            dt = passive;
        end
        
        for sv = sUse
            h_leg(ddi) = errorbar(tm + tshifts(ddi), dt.nsensitivity_mean(sv, :, subjIdx), dt.nsensitivity_se(sv, :, subjIdx),...
                [symbs(ddi), '-'], 'markersize', markersize, 'capsize', capsize, 'Color', cols(ddi, :),...
                'MarkerFaceColor', cols_fill(ddi, :), 'linewidth', 1);
        
            for tv = 1:size(pvals, 2)
                if pvals(sv, tv, subjIdx) < .05
                    plot(tm(tv), 1.1, 'k*');
                end
            end
        end
    end
    title(ttls{subjIdx});
    xlabel('time from saccade or motion onset (ms)');
    ylabel('normalized sensitivity');
    xlim(active.TmpVec([1, end]));
    xticks(active.TmpVec);
    set(gca, 'FontSize', 12);
    ylim([0.6, 1.2]);
end
legend(h_leg, {'Active', 'Passive'}, 'Location', 'southwest');