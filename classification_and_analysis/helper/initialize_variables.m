%% Script to initialize a set of variables from a list of string variable names
% Intended for use with checkboxList.m to initialize a list of variables
%
% INPUTS: variables_to_create, target_values

%% Check for bad input
if length(variables_to_create) ~= length(target_values)
    error("Variables and targets list must be same length.")
end

%% Create variables
disp("Initializing Variables:")
for i=1:length(variables_to_create)
    if target_values(i) %Fuck matlab doesn't have ternary operators
        value = "true";
    else
        value = "false";
    end
    
    eval(strcat(variables_to_create(i), " = ", value, ";"));
    disp(strcat(" - ", variables_to_create(i), " : ", value))
end

%% Clear extra variables
clear value