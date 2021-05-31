function outlier_ids = find_univariate_outliers(data, std_threshold)
    % Find univariate outliers given data vector `data` and the
    % `std_threshold`. Returns a boolean array (True if observation is an
    % outlier).
    outlier_ids = find(data < (mean(data) - std_threshold * std(data)) | ...
                       data > (mean(data) + std_threshold * std(data)));           
end