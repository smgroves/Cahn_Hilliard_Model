% Define paths to the datasets and corresponding t0 values
% /Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r "plot_inflection_pt();quit;"

alpha = "0.0";
indir =sprintf("/Users/smgroves/Documents/GitHub/Cahn_Hilliard_Model/Cahn_Hilliard_solvers/julia_multigrid/manuscript_output/critical_radius/alpha_%s", alpha);
% epsilon ="0.015009";

% for epsilon = ['0.060037'] 
    plot_inflection('0.0075047', indir, alpha);
% end

function [] = plot_inflection(epsilon, indir, alpha)
    path = sprintf("%s/radius_0.5_level_set_epsilon_%s_alpha_%s_nx_256.txt", indir, epsilon, alpha);
    all_data = readtable(path);
    % Predefine colors for each dataset for consistency in plotting
    colors = {'b', 'r', 'g', 'c', 'm', 'y', 'k',"#d142f5","#f542b6","#f59e42","#42f5d7","#42b3f5"};

    % kth degree polynomial fit
    k = 10;

    % Initialize figure
    figure;
    hold on;
    R0s = unique(all_data.R0)
    for i = 1:length(R0s)
        try
            R0 = R0s(i)
            % do stuff with i, item
            % Load the dataset
            data = all_data(all_data.R0==R0,:);
            tt = data.time;
            rr = data.radius;
            

            % Process the dataset
            % Step 1: Find rows with NaN in 'radius'
            nan_rows = isnan(data.radius);

            % Step 2: Find index of the first NaN
            first_nan_index = find(nan_rows, 1);

            % Step 3: Retrieve value of 'AnotherColumn' at that index
            t0 = data.time(first_nan_index)
            % t0 = datasets{i}.t0;

            idx = find(tt > .9*t0, 1);

            % Ensure tt and rr are column vectors and select data up to idx
            tt = tt(1:idx)';
            rr = rr(1:idx)';

            % Fit a polynomial to the dataset
            p = polyfit(tt, rr, k); % kth degree polynomial fit
            rr_fit = polyval(p, tt);

            % Calculate the first derivative of the fitted polynomial
            p_deriv1 = polyder(p);
            rr_deriv1 = polyval(p_deriv1, tt);

            % Calculate the second derivative of the fitted polynomial
            p_deriv2 = polyder(p_deriv1);
            rr_deriv2 = polyval(p_deriv2, tt);

            % Find the inflection points (where the second derivative changes sign)
            inflection_idx = find(diff(sign(rr_deriv2)) ~= 0); % Indices where second derivative changes sign

            % Get the x and y coordinates of the inflection points
            inflection_tt = tt(inflection_idx);
            inflection_rr = rr_fit(inflection_idx);

            % Plot the data and the fitted curve
            % [~, name, ~] = fileparts(datasets{i}.path); % Extract the filename
            % name = strrep(name, '_', '='); % Replace all underscores with equal signs
            name=sprintf("R0 = %.5g",R0);
            plot(tt, rr, '.', 'Color', colors{i}, 'DisplayName', sprintf('Data %s', name));
            plot(tt, rr_fit, '-', 'Color', colors{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('Fitted Curve %s', name));

            % Highlight inflection points
            plot(inflection_tt, inflection_rr, 'o', 'Color', colors{i}, 'MarkerFaceColor', colors{i}, 'DisplayName', sprintf('Inflection Points %s', name));
            
            % Label the inflection points
            if ~isempty(inflection_tt)
                [sorted_tt, sort_idx] = sort(inflection_tt);
                sorted_rr = inflection_rr(sort_idx);

                % Highlight the first inflection point
                text(sorted_tt(1), sorted_rr(1), sprintf(' (%.2g, %.4g)', sorted_tt(1), sorted_rr(1)), ...
                    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

                % Highlight the more inflection points if it exists
                if length(sorted_tt) == 2
                    text(sorted_tt(2), sorted_rr(2), sprintf(' (%.2g, %.4g)', sorted_tt(2), sorted_rr(2)), ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
                elseif length(sorted_tt) == 3
                    text(sorted_tt(3), sorted_rr(3), sprintf(' (%.2g, %.4g)', sorted_tt(3), sorted_rr(3)), ...
                        'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
                end
            end
        catch
            warning('No inflection point. Skipping this R0.');
        end
    end

    xlabel('Time');
    ylabel('Radius');
    title(sprintf('Curve and Inflection Points for Multiple Datasets, epsilon = %s', epsilon));
    legend show;
    grid on;
    hold off;

    % saveas(gcf,sprintf('%s/%s_inflection_points.pdf', indir, epsilon));
    % exportgraphics(gcf,sprintf('%s/%s_inflection_points.pdf', indir, epsilon),'Resolution',300) 
    set(gcf, 'PaperSize', [11, 20])
    orient(gcf,'landscape')

    print(gcf,sprintf('%s/%s_inflection_points_nx_256.pdf', indir, epsilon),"-dpdf",'-fillpage')
end
