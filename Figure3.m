%%% Plot raw and normalized contrast sensitivity aligned to saccade offset
%%% and saccade peak speed. Show the corresponding eye speed profiles.
%%% author: Janis Intoy

load('contrast_thresholds.mat', 'sacc_offset', 'sacc_pkspeed', 'fixation');

%% plotting parameters
markersize = 6;
capsize = 8;
col = [1 0 0;0 0 1; 0 0.5 0];
star_height = 2.9;
bar_height = [1.27, 1.24, 1.21];

speed_tlim = 200;

for ddi = 1:2 % loop through saccade offset or peak speed alignment
    if ddi == 1
        dt = sacc_offset;
    else
        dt = sacc_pkspeed;
    end
    
    %% Figure 3A/B - raw sensitivity values and speed profile
    SpVec = dt.SpVec;
    tm = dt.times_average;
    
    %%%%%% compute sensitivity
    dt.sensitivity = 1 ./ dt.thresh;
    
    dt.sensitivity_mean = mean(dt.sensitivity, 3);
    dt.sensitivity_se = std(dt.sensitivity, [], 3) / sqrt(size(dt.sensitivity, 3));
    
    fixation.sensitivity = 1 ./ fixation.thresh;
    
    fixation.sensitivity_mean = mean(fixation.sensitivity, 3);
    fixation.sensitivity_se = std(fixation.sensitivity, [], 3) / sqrt(size(fixation.sensitivity, 3));
    
    %%%%%% one-way ANOVAs
    
    % At each time, does sensitivity vary with eccentricity?
    p_ecc = nan(size(dt.thresh, 2), 1);
    for tv = 1:length(p_ecc)
        p_ecc(tv) = anova1(squeeze(dt.sensitivity(:, tv, :))', [], 'off');
    end
    
    % At each eccentricity, when does sensitivity differ from fixation?
    p_compareToFixation = nan(size(dt.sensitivity_mean));
    for sv = 1:size(p_compareToFixation, 1)
        tmp = [squeeze(fixation.sensitivity(sv, 1, :))'; squeeze(dt.sensitivity(sv, :, :))];
        [~, ~, stats] = anova1(tmp', 1:9, 'off'); %
        c = multcompare(stats, 'Display', 'off');
        
        p_compareToFixation(sv, :) = c(1:length(dt.TmpVec)-1, end);
    end
    
    %%%%%% plot
    pp = nan(length(SpVec)-1, 1);
    ss = cell(length(SpVec)-1, 1);
    
    figure((ddi - 1) * 2 + 1); clf; 
    axes('Position', [0.2, .35 .75, .6]); hold on;
    for sv = 1:length(SpVec)-1
        pp(sv) = plot(tm(sv, :), dt.sensitivity_mean(sv, :),...
            '-d', 'color', col(sv,:),...
            'MarkerSize', markersize , 'MarkerFaceColor', col(sv,:),...
            'markerEdgeColor', 'none');
        errorbar(tm(sv, :), dt.sensitivity_mean(sv, :), dt.sensitivity_se(sv,:),...
            'color', col(sv,:), 'capsize', capsize);
        ss{sv} = sprintf('%i-%i', SpVec(sv), SpVec(sv+1));
    end
    
    % plot significance
    for tv = 1:length(p_ecc)
        if p_ecc(tv) < 0.05
            plot(mean(tm(:, tv)), star_height, 'k*');
        end
    end
    
    for sv = 1:size(p_compareToFixation, 1)
        for tv = 1:size(p_compareToFixation, 2)
            if p_compareToFixation(sv, tv) < 0.05
                plot(dt.TmpVec(tv:tv+1), bar_height(sv)*[1, 1], '-', ...
                    'Color', col(sv, :), 'LineWidth', 2);
            end
        end
    end
    
    ylim([1.2, 3]);
    set(gca,'YTick',1:0.5:3, 'YScale', 'log');
    plot([0, 0], [1, 3], 'k--');
    xlim([-250, dt.TmpVec(end)]);
    set(gca, 'XTick', dt.TmpVec, 'XTickLabels', []);
    lb = legend(pp, ss,'Location', 'southwest');
    
    ylabel('contrast sensitivity');
    
    yl = [0, 40];
    axes('Position', [0.2, .17 .75, .15]); hold on;
    [hl, hp] = boundedline(-speed_tlim:speed_tlim, dt.eye_speed(:, 1)/60, dt.eye_speed(:, 2)/60, 'k');
    ylim(yl);
    
    plot([0, 0], yl, 'k--');
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    ylabel({'speed', '(deg/s)'},'FontName', 'Arial', 'FontSize', 12);
    
    if ddi == 1
        xlabel('time from saccade offset (ms)');
    else
        xlabel('time from saccade peak speed (ms)');
    end
    
    %% Figure 3C/D: contrast sensitivity normalized to fixation
    %%%%%% compute sensitivity ratio
    dt.sensitivity_norm = bsxfun(@rdivide, ...
        bsxfun(@minus, dt.sensitivity, dt.sensitivity(:, 1, :)),...
        dt.sensitivity(:, 1, :));
    dt.sensitivity_norm_mean = mean(dt.sensitivity_norm, 3);
    dt.sensitivity_norm_se = std(dt.sensitivity_norm, [], 3) / sqrt(size(dt.sensitivity_norm, 3));
    
    %%%%%% When does the change in sensitivity vary across eccentricities?
    p_change_ecc = nan(size(tm, 2),1 );
    for tv = 1:size(tm, 2)
        tmp = squeeze(dt.sensitivity_norm(:, tv, :));
        p_change_ecc(tv) = anova1(tmp', ss, 'off');
    end
    
    %%%%%% plot
    figure((ddi - 1) * 2 + 2); clf; hold on;
    plot(dt.TmpVec([1, end]), [0, 0], 'k--');
    for sv = 1:length(SpVec)-1
        pp(sv) = plot(tm(sv, :), dt.sensitivity_norm_mean(sv, :),...
            '-d', 'color', col(sv,:),...
            'MarkerSize', markersize , 'MarkerFaceColor', col(sv,:),...
            'markerEdgeColor', 'none');
        errorbar(tm(sv, :), dt.sensitivity_norm_mean(sv, :), dt.sensitivity_norm_se(sv,:),...
            'color', col(sv,:), 'capsize', capsize);
        ss{sv} = sprintf('%i-%i', SpVec(sv), SpVec(sv+1));
    end
    
    for tv = 1:length(p_change_ecc)
        if p_change_ecc(tv) < 0.05
            plot(mean(tm(:, tv)), 0.25, 'k*');
        end
    end
    
    ylim([-.6, 0.3]);
    set(gca,'YTick', -.6:.2:0.3);
    plot([0, 0], [-.6, .3], 'k--');
    xlim(dt.TmpVec([1, end]));
    set(gca, 'XTick', dt.TmpVec);
    
    ylabel('change in sensitivity');
    if ddi == 1
        xlabel('time from saccade offset (ms)');
    else
        xlabel('time from saccade peak speed (ms)');
    end
end