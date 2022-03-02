%%% Compare contrast sensitivity before and after saccades. Show
%%% differences in pre-saccadic sensitivity for foveated versus
%%% non-foveated probes.
%%% author: Janis Intoy

load('contrast_thresholds.mat', 'presaccadic', 'postsaccadic', 'foveated', 'nonfoveated');

markersize = 6;
capsize = 8;

%% Figure 4A
for ddi = 1:2
    if ddi == 1
        dt = presaccadic;
    else
        dt = postsaccadic;
    end
    
    % compute sensitivity
    dt.sensitivity = 1 ./ dt.thresh;
    dt.sensitivity_mean = mean(dt.sensitivity, 3);
    nnn = sum(~isnan(dt.thresh), 3);
    dt.sensitivity_se = std(dt.sensitivity, [], 3) ./ sqrt(nnn);
    
    if ddi == 1
        presaccadic = dt;
    else
        postsaccadic = dt;
    end
end

% plot
symbols = {'s-', 'o--'};

coli = hsv(6);
col = zeros(2, 3);

sens1 = presaccadic.sensitivity;
sens2 = postsaccadic.sensitivity;

m1 = presaccadic.sensitivity_mean;
s1 = presaccadic.sensitivity_se;

m2 = postsaccadic.sensitivity_mean;
s2 = postsaccadic.sensitivity_se;

h_legs = nan(2, 1);

figure(); hold on;
for sv = 1:size(m1, 1)
    hh = plot((1:2) + .1 * sign(sv - 1.5), [squeeze(sens1(sv, :, :)), squeeze(sens2(sv, :, :))],...
        symbols{sv}(1),...
        'markersize', 5);
    for hi = 1:length(hh)
        set(hh(hi), 'MarkerFaceColor', (coli(hi, :) + [2, 2, 2]) / 3); 
        set(hh(hi), 'MarkerEdgeColor', 'none');
    end
    
    [~, pval] = ttest(squeeze(sens1(sv, :, :)), squeeze(sens2(sv, :, :)));
    if pval < 0.05
        plot(1.5, (m1(sv, :) + m2(sv, :)) / 2 + .05, '*',...
            'Color',  col(sv, :), 'MarkerSize', 5);
    else
        text(1.5, (m1(sv, :) + m2(sv, :)) / 2 + .05, 'n.s.',...
            'Color',  col(sv, :));
    end
    
    h_legs(sv) = errorbar((1:2) + .05 * sign(sv - 1.5), [m1(sv, :), m2(sv, :)], [s1(sv, :), s2(sv, :)],...
        symbols{sv}, 'Color', col(sv, :),...
        'MarkerFaceColor', col(sv, :),...
        'MarkerEdgeColor', 'none',...
        'capsize', capsize,...
        'markersize', 7);
    
end
xlim([0.5, 2.5]);
ylim([1.8, 2.5]);
set(gca, 'YScale', 'log');
xticks(1:2);
xticklabels({'before saccade', 'after saccade'});
set(gca, 'FontSize', 10, 'FontName', 'Arial');
ylabel('contrast sensitivity', 'FontSize', 12);
legend(h_legs, {'0''-20''', '20''-60'''}, 'Location', 'south');

%% Figure 4B-C
sUse = 2:3; % only use peripheral spatial bins
tUse = 1:4; % only compare pre-saccadic time bins

for ddi = 1:2
    if ddi == 1
        dt = foveated;
    else
        dt = nonfoveated;
    end
    
    % compute sensitivity and CI from threshold values
    dt.sensitivityBS = 1 ./ dt.threshBS;
    dt.sensitivity_mean = nanmean(dt.sensitivityBS, 4);
    for sv = sUse
        for tv = tUse
            dt.sensitivity_95(sv, tv, :) = cat(3,...
                quantile(dt.sensitivityBS(sv, tv, :, :), .025, 4),...
                quantile(dt.sensitivityBS(sv, tv, :, :), .975, 4));
        end
    end
    
    if ddi == 1
        foveated = dt;
    else
        nonfoveated = dt;
    end
end

% non-parametric one-sided bootstrap tests
pvals = nan(length(foveated.SpVec)-1, length(foveated.TmpVec)-1);
for sv = sUse
    for tv = tUse
        dt1 = squeeze(foveated.sensitivityBS(sv, tv, :, :));
        dt2 = squeeze(nonfoveated.sensitivityBS(sv, tv, :, :));
        
        use = ~isnan(dt1) & ~isnan(dt2);
        n = sum(use);
        pvals(sv, tv) = nanmean(dt1(use) < dt2(use)); % one-sided
    end
end

% plot
tm = (foveated.TmpVec(1:end-1) + foveated.TmpVec(2:end))/2;

col = cool(2);
for sv = sUse
    figure(); hold on;
    for ddi = 1:2
        if ddi == 1
            dt = foveated;
        else
            dt = nonfoveated;
        end
        sens = dt.sensitivity_mean(sv, :, :);
        senslow = dt.sensitivity_95(sv, tUse, 1);
        senshigh = dt.sensitivity_95(sv, tUse, 2);
        errorbar(tm(tUse), sens(tUse), ...
            sens(tUse) - senslow(tUse), ...
            senshigh(tUse) - sens(tUse),...
            '-d', ...
            'linewidth', 1, 'MarkerSize', markersize, ...
            'capsize', capsize,...
            'Color', col(ddi, :), ...
            'MarkerFaceColor', col(ddi, :), 'MarkerEdgeColor', col(ddi, :));
    end
    xlim([foveated.TmpVec(1), 0]);
    set(gca,'XTick',foveated.TmpVec,'XTickLabel',foveated.TmpVec(1:end));
    set(gca, 'yscale', 'log','FontSize', 12)
    ylim([1 3]);
    set(gca,'YTick',1:0.5:3);
    set(gca, 'FontSize', 10, 'FontName', 'Arial');
    ylabel('contrast sensitivity','FontSize',12);
    xlabel('time from saccade onset (ms)','FontSize',12);
    
    legend({'Foveated', 'Non-foveated'}, 'Location', 'southeast');
    
    for tv = tUse
        if pvals(sv, tv) < 0.05
            text( tm(tv), 2.7, '*', 'FontSize', 12);
        end
    end
    
    text (-170, 1.1,...
        sprintf('%i''-%i''', foveated.SpVec(sv), foveated.SpVec(sv+1)),'FontSize',14);
end
