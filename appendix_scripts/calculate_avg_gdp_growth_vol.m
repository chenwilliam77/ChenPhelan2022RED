% Calculate average GDP growth over 1950-2019
%
% Written by William Chen, Jan. 2022

% Grab US data and define auxiliary functions
input_data_fp = '../save/input_data/';
pwt = readtable([input_data_fp, 'pwt100.xlsx'], 'Sheet', 'Data');
us_pwt_data = pwt(strcmp(pwt.countrycode, 'USA'), :);
growth_fnct = @(x) (x - lagmatrix(x, 1)) ./ lagmatrix(x, 1);

% GDP growth statistics
gdp_growth = growth_fnct(us_pwt_data.rgdpna);
disp(['Mean GDP Growth = ', ...
    num2str(100 * mean(gdp_growth, 'omitnan')), '%']);
disp(['Vol. GDP Growth = ', ...
    num2str(100 * std(gdp_growth, 'omitnan')), '%']);
