function formattedStr = floatToString(num)
    % Convert float to string with up to 3 significant digits
    str = sprintf('%.5g', num)
    % Remove trailing zeros after the decimal point
    if contains(str, '.') % Check if there is a decimal point
        idx = strfind(str, '.');
        if ~isempty(idx)
            % Find the last non-zero digit after the decimal point
            endIdx = idx + numel(str) - 1;
            while endIdx > idx && str(endIdx) == '0'
                endIdx = endIdx - 1;
            end
            if str(endIdx) == '.'
                formattedStr = str(1:endIdx-1); % Remove the decimal point if it's the last character
            else
                formattedStr = str(1:endIdx);
            end
        else
            formattedStr = str;
        end
    else
        formattedStr = str;
    end
end