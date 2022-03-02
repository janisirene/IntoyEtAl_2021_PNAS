%%% Script to fit psychometric curves and estimate contrast thresholds
%%% using example data as shown in Figure 2B.
%%% author: Janis Intoy

%%% This script requires the psignifit toolbox for psychometric fitting 
%%% implemented by the Wichmann lab (Wichmann and Hill, 2001). The toolbox 
%%% is available for download from github:
%%% https://github.com/wichmann-lab/psignifit (Aug 11, 2021)

%addpath('psignifit-master');
load('psychometric_data.mat');

nBoots = 50; % 500 was used in manuscript figures

%% psychometric fit parameters
Threshlev = 0.25;

options             = struct;   % initialize as an empty struct

options.sigmoidName = 'logn';
options.expType     = 'YesNo';   
                       
options.threshPC = Threshlev;
options.stimulusRange = [0.25, .7]; % used to set the prior on the threshold

options.fixedPars = NaN(5,1); %(threshold, width, lapse, guess and eta)
options.fixedPars(3:end) = 0; % only fit mean and sigma

%%
cols = 'rbg';
lvl_plot = logspace(-0.6021, -.1549, 100);

for pp = 1:length(psychometric_data)
    lvl_all = psychometric_data(pp).contrast;
    hits_all = psychometric_data(pp).hits;
    
    p = nan(nBoots, length(lvl_plot)); % pschometric curves
    thr = nan(nBoots, 1); % contrast thresholds
    
    for bi = 1:nBoots+1
        rs = randsample(length(lvl_all), length(lvl_all), true);
        
        if bi <= nBoots
            lvl = lvl_all(rs);
            hits = hits_all(rs);
        else % last iteration uses all data
            lvl = lvl_all;
            hits = hits_all;
        end
        
        % reformat data for psignifit
        x = double( unique( lvl(lvl>0) ) );
        k = zeros(size(x));
        n = zeros(size(x));
        for ii=1:size(x,2)
            k(ii) = sum(  hits ( (x(ii)==lvl) )  );
            n(ii) = sum(  (x(ii)==lvl)  );
        end
        xkn = [x(:), k(:), n(:)];
        
        % fit psychometric curve
        result = psignifit(xkn, options);
        
        % estimate full psychometric curve with these parameters
        fitValues = (1-result.Fit(3)-result.Fit(4))*arrayfun(@(x) result.options.sigmoidHandle(x,result.Fit(1),result.Fit(2)),lvl_plot)+result.Fit(4);
        
        if bi <= nBoots
            thr(bi) = exp(result.Fit(1));
            p(bi, :) = fitValues;
        end
    end

    p_mean = mean(p, 1);
    p_se = std(p, [], 1);
    
    thr_mean = mean(thr);
    thr_se = std(thr);
    
    n_size = (n - min(n)) / (max(n) - min(n)) * 30 + 20;
    
    figure(pp); clf; hold on;
    scatter(x, k ./ n, n_size, [0.3, .3 .3], 'filled');
    [hl, hp] = boundedline(lvl_plot, p_mean, p_se);
    plot(thr_mean + [-thr_se, thr_se], [0.25, 0.25], 'linewidth', 2, 'Color', 'k');
    plot([lvl_plot(1), thr_mean], [Threshlev, Threshlev], 'k--');
    plot(thr_mean*[1,1], [0, 0.25], 'k--');
    set(hl, 'Color', cols(pp), 'LineWidth', 2);
    set(hp, 'FaceColor', cols(pp), 'FaceAlpha', .3);
    yticks(0:.25:1);
    xlim(lvl_plot([1, end]));
    ylim([-.05, 1.05]);
    
    xlabel('michelson contrast');
    ylabel('proportion correct');
    
    sv = psychometric_data(pp).spatialIndex;
    tv = psychometric_data(pp).temporalIndex;
    title(sprintf('%i''-%i'', %i - %i ms', SpVec(sv), SpVec(sv+1), TmpVec(tv), TmpVec(tv+1)));
end