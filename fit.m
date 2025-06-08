function FIT = fit(y_true, y_model)
    % This function calculates the fitting percentage (FIT).
    % Input:
    %   y_true  - Actual data (vector)
    %   y_model - Predicted model data (vector)
    % Output:
    %   FIT - Fitting percentage
    
    % Number of data points
    N = length(y_true);
    
    % Mean of true data
    y_true_mean = mean(y_true);
    
    % Calculate FIT
    FIT = 100 * (1 - sum((y_true - y_model).^2) / sum((y_true - y_true_mean).^2));
    
    % Display FIT value
    disp(['Fitting Percentage (FIT): ', num2str(FIT), '%']);
end
