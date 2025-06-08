function PEC = pec(y_true, y_model)
    % This function calculates the prediction error cost (PEC).
    % Input:
    %   y_true  - Actual data (vector)
    %   y_model - Predicted model data (vector)
    % Output:
    %   PEC - Prediction Error Cost
    
    % Number of data points
    N = length(y_true);
    
    % Calculate PEC
    PEC = sum((y_true - y_model).^2) / N;
    
    % Display PEC value
    disp(['Prediction Error Cost (PEC): ', num2str(PEC)]);
end
