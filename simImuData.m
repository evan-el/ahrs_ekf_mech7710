function [data, orientationNED, mag_field_vec] = simImuData(dt)
%SIMIMUDATA Summary of this function goes here
%   Detailed explanation goes here
    fs = 1/dt;
    firstLoopNumSamples = fs*4;
    secondLoopNumSamples = fs*2;
    totalNumSamples = firstLoopNumSamples + secondLoopNumSamples;
    
    traj = kinematicTrajectory('SampleRate',fs);
    
    accBody = zeros(totalNumSamples,3);
    angVelBody = zeros(totalNumSamples,3);
    angVelBody(1:firstLoopNumSamples,3) = (2*pi)/4;
    angVelBody(firstLoopNumSamples+1:end,3) = (2*pi)/2;
    
    [~,orientationNED,~,accNED,angVelNED] = traj(accBody,angVelBody);
    
    IMU = imuSensor('accel-gyro-mag','SampleRate',fs);
    
    [accelReadings, gyroReadings, magReadings] = IMU(accNED,angVelNED,orientationNED);
    
    data = [accelReadings, gyroReadings, magReadings];
    orientationNED = rad2deg(quat2eul(orientationNED, "XYZ"));

end

