%% Plot results

data_dirs = {<Insert a list of paths to datafiles produced by the *_experiment.m files.>};

styles = {'-', '-.', '--', '-.', '--', ':', ':'};
CM = colormap('lines');
colors = [CM(1:5, :); CM(1:2, :)];
figure(1), clf
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultTextFontName', 'Times New Roman')
for j = 1:length(data_dirs)
    % Stack all accs.
    data_dir = data_dirs{j};
    accs = [];
    n_reps = 0;
    files = dir(strcat(data_dir, '*.mat'));
    for file = files'
        ws = load(strcat(data_dir, file.name));
        accs_part = ws.accs;
        if length(size(accs_part)) == 3
            n_reps = n_reps + size(accs_part, 3);
        else
            n_reps = n_reps + 1;
        end
        if length(accs) == 0
            accs = accs_part;
        else
            accs = cat(3, accs, accs_part);
        end
    end
    mean_accs = mean(accs, 3);
    fprintf('%d repetitions.\n', n_reps);
    fprintf('Solver: %d\n', ws.solver);
    fprintf('k=%d\n', ws.k);

    subplot(1, length(data_dirs), j);
    legends = {};
    for i = 1:length(ws.method_names)-1
        style = styles{mod(i-1, length(styles)) + 1};
        method = ws.method_names{i};
        if iscell(ws.query_counts)
            query_counts = ws.query_counts{i};
        else
            query_counts = ws.query_counts;
        end
        plot(query_counts, mean_accs(i,1:length(query_counts)), ...
            style, 'color', colors(i,:), 'Linewidth',1.5), hold on
        if strcmp(method, 'TopMatchings10Batch')
            method = 'TopMatchingsBatch';
        end
        legends{end+1} = method;
    end
    if j == 1
        legend(legends, 'Location', 'SouthEast', 'FontSize', 9)
    end
    set(gca, 'FontSize', 13)
    xlabel('# of queried nodes', 'fontsize', 15);
    ylabel('Accuracy (unqueried nodes)', 'fontsize', 15);
    ylim([0.5, 1]);
    if ws.solver == 1
        title('NetAlignMP++', 'fontsize', 12);
    else
        title('Natalie', 'fontsize', 12);
    end
end
%fname = 'social_networks.eps';
%fname = 'family_trees.eps';
fname = 'synthetic.eps';
