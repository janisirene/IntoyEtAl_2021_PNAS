%%% Plot raw and normalized contrast sensitivity aligned to saccade onset 
%%% in three spatial bins and eight temporal bins. Compute the strength of
%%% peri-microsaccadic suppression.
%%% author: Janis Intoy

load('contrast_thresholds.mat', 'sacc_onset', 'fixation');

%% plotting parameters
markersize = 6;
capsize = 8;
col = [1 0 0;0 0 1; 0 0.5 0];
star_height = 2.7;
bar_height = [1.27, 1.24, 1.21];

SpVec = sacc_onset.SpVec;
tm = sacc_onset.times_average;

%% Figure 2C: sensitivity aligned to saccade onset
%%%%%% compute sensitivity
sacc_onset.sensitivity = 1 ./ sacc_onset.thresh;

sacc_onset.sensitivity_mean = mean(sacc_onset.sensitivity, 3);
sacc_onset.sensitivity_se = std(sacc_onset.sensitivity, [], 3) / sqrt(size(sacc_onset.sensitivity, 3));

fixation.sensitivity = 1 ./ fixation.thresh;

fixation.sensitivity_mean = mean(fixation.sensitivity, 3);
fixation.sensitivity_se = std(fixation.sensitivity, [], 3) / sqrt(size(fixation.sensitivity, 3));

%%%%%% one-way ANOVAs

% Does fixational sensitivity vary with eccentricity?
[p_fixation, tbl_fixation] = anova1(squeeze(fixation.sensitivity)', [], 'off');

% At each time, does sensitivity vary with eccentricity?
p_ecc = nan(size(sacc_onset.thresh, 2), 1);
for tv = 1:length(p_ecc)
     p_ecc(tv) = anova1(squeeze(sacc_onset.sensitivity(:, tv, :))', [], 'off');
end

% At each eccentricity, when does sensitivity differ from fixation?
p_compareToFixation = nan(size(sacc_onset.sensitivity_mean));
for sv = 1:size(p_compareToFixation, 1)
    tmp = [squeeze(fixation.sensitivity(sv, 1, :))'; squeeze(sacc_onset.sensitivity(sv, :, :))];
    [~, ~, stats] = anova1(tmp', 1:9, 'off'); %
    c = multcompare(stats, 'Display', 'off');
    
    p_compareToFixation(sv, :) = c(1:length(sacc_onset.TmpVec)-1, end);
end

%%%%%% plot
pp = nan(length(SpVec)-1, 1);
ss = cell(length(SpVec)-1, 1);

tm_fixation = -220; % for plotting purposes

figure(1); clf; hold on;
for sv = 1:length(SpVec)-1
    pp(sv) = plot(tm(sv, :), sacc_onset.sensitivity_mean(sv, :),...
        '-d', 'color', col(sv,:),...
        'MarkerSize', markersize , 'MarkerFaceColor', col(sv,:),...
        'markerEdgeColor', 'none');
    errorbar(tm(sv, :), sacc_onset.sensitivity_mean(sv, :), sacc_onset.sensitivity_se(sv,:),...
        'color', col(sv,:), 'capsize', capsize);
    ss{sv} = sprintf('%i-%i', SpVec(sv), SpVec(sv+1));
    
    
    errorbar(tm_fixation,fixation.sensitivity_mean(sv,:),fixation.sensitivity_se(sv,:),...
        'd', 'markerSize', markersize , 'MarkerFaceColor', col(sv, :),...
        'markerEdgeColor', 'none',...
        'color',col(sv,:), 'capsize', capsize);
end

% plot significance
if p_fixation < 0.05
    plot(tm_fixation, star_height, 'k*'); 
end

for tv = 1:length(p_ecc)
    if p_ecc(tv) < 0.05
        plot(mean(tm(:, tv)), star_height, 'k*');
    end
end

for sv = 1:size(p_compareToFixation, 1)
    for tv = 1:size(p_compareToFixation, 2)
        if p_compareToFixation(sv, tv) < 0.05
            plot(sacc_onset.TmpVec(tv:tv+1), bar_height(sv)*[1, 1], '-', ...
                'Color', col(sv, :), 'LineWidth', 2);
        end
    end
end

ylim([1.2, 3]);
set(gca,'YTick',1:0.5:3, 'YScale', 'log');
plot([0, 0], [1, 3], 'k--');
xlim([-250, sacc_onset.TmpVec(end)]);
set(gca, 'XTick', sacc_onset.TmpVec);
lb = legend(pp, ss,'Location', 'southwest');

ylabel('contrast sensitivity');
xlabel('time from saccade onset (ms)');

%% Figure 2D: contrast sensitivity normalized to fixation
%%%%%% compute sensitivity ratio
sacc_onset.sensitivity_norm = bsxfun(@rdivide, sacc_onset.sensitivity, fixation.sensitivity);
sacc_onset.sensitivity_norm_mean = mean(sacc_onset.sensitivity_norm, 3);
sacc_onset.sensitivity_norm_se = std(sacc_onset.sensitivity_norm, [], 3) / sqrt(size(sacc_onset.sensitivity_norm, 3));

%%%%%% plot
figure(2); clf; hold on;
plot(sacc_onset.TmpVec([1, end]), [1, 1], 'k--');
for sv = 1:length(SpVec)-1
    pp(sv) = plot(tm(sv, :), sacc_onset.sensitivity_norm_mean(sv, :),...
        '-d', 'color', col(sv,:),...
        'MarkerSize', markersize , 'MarkerFaceColor', col(sv,:),...
        'markerEdgeColor', 'none');
    errorbar(tm(sv, :), sacc_onset.sensitivity_norm_mean(sv, :), sacc_onset.sensitivity_norm_se(sv,:),...
        'color', col(sv,:), 'capsize', capsize);
    ss{sv} = sprintf('%i-%i', SpVec(sv), SpVec(sv+1));
end

ylim([0.5 1.1]);
set(gca,'YTick',[0.5, 0.75, 1]);
plot([0, 0], [0.5 1.1], 'k--');
xlim(sacc_onset.TmpVec([1, end]));
set(gca, 'XTick', sacc_onset.TmpVec);

ylabel('normalized sensitivity');
xlabel('time from saccade onset (ms)');

%% Figure 2E: strength of suppression
tm = (sacc_onset.TmpVec(1:end-1) + sacc_onset.TmpVec(2:end)) / 2;
perisaccadic = tm > -50 & tm < 50;

suppression_strength = bsxfun(@rdivide,...
    bsxfun(@minus, fixation.sensitivity, sacc_onset.sensitivity), ...
    fixation.sensitivity);

% average over perisaccadic times
suppression_strength = squeeze( mean(suppression_strength(:, perisaccadic, :), 2) );

mean_suppression = mean(suppression_strength, 2);
std_suppression = std(suppression_strength, [], 2);

% Does perisaccadic suppression differ between eccentricities?
[p_strength, tbl, stats] = anova1(suppression_strength', [], 'off');
c = multcompare(stats, 'Display', 'off');

figure(3); clf; 
axes('Position', [0.35, .1, .3, .8]); hold on;
bar(1:length(mean_suppression), 100 * mean_suppression, 'FaceColor', [0.5 0.5 0.5]);
h_e = errorbar(1:length(mean_suppression), 100 * mean_suppression, 100 * std_suppression,...
    'k');
set(h_e, 'LineStyle', 'none');
ylabel('suppression strength');

yticks(0:10:50); yticklabels({'0%', '', '', '', '40%', ''});
ylim([0, 50]);
xlim([0.5, 3.5]);

xticks(1:length(ss));
xticklabels(ss);

% significant different between 1st and 3rd spatial bins
if c(2, end) < 0.05
    plot([1, 3], 45*[1 1], 'k-'); 
    text(2, 45, '*', 'FontSize', 20);
end