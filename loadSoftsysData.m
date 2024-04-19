function [data, dt] = loadSoftsysData(filename)
%LOADSOFTSYSDATA Summary of this function goes here
%   Detailed explanation goes here

    % columns 1-3 accel [m/s^2]
    % 4-6 gyro [rad/s]
    % 7-9 mag [uT]
    data = readmatrix(filename);
    dt = 0.0078; % from hw description
end

