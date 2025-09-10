% 
% Created by Rhet Hailey - 4/28/24
%
% hard coded as shir

clear, clc

addpath '/home/rhet/School/Spring24/MECH7710/final_project/Data'

static = rosbag("2024-04-25-05-06-19_static1.bag");
dynamic = rosbag("dynamic1_2024-04-25-05-21-27.bag");
magCal = rosbag("mag_cal_4_25_2024-04-25-05-10-44.bag");


list = {'static','dynamic','magCal'};
listBag = {static,dynamic,magCal};

data = struct();

for i = 1:length(list)
    
    % from what bag
    imu = select((listBag{i}),'Topic','/razor/imu');
    mag = select((listBag{i}),'Topic','/razor/mag');

    imuStruct = readMessages(imu,"DataFormat","struct");
    magStruct = readMessages(mag,"DataFormat","struct"); % maybe later
    
    % time
    time = imu.MessageList.Time; time = time - time(i);
    data.((list{i})).time = time;

    type = "Orientation";
    vars = {"X","Y","Z","W"};
    data.((list{i})).orient = getXYZFromPoint(imuStruct,type,vars);

    type = "OrientationCovariance";
    data.((list{i})).orientCov = getArrFromPoint(imuStruct,type);

    type = "AngularVelocity";
    vars = {"X","Y","Z"};
    data.((list{i})).omega = getXYZFromPoint(imuStruct,type,vars);

    type = "AngularVelocityCovariance";
    data.((list{i})).omegaCov = getArrFromPoint(imuStruct,type);

    type = "LinearAcceleration";
    vars = {"X","Y","Z"};
    data.((list{i})).acc = getXYZFromPoint(imuStruct,type,vars);

    type = "LinearAccelerationCovariance";
    data.((list{i})).accCov = getArrFromPoint(imuStruct,type);

    type = "MagneticField_";
    vars = {"X","Y","Z"};
    data.((list{i})).mag = getXYZFromPoint(magStruct,type,vars);

    type = "MagneticFieldCovariance";
    data.((list{i})).magCov = getArrFromPoint(magStruct,type);

end

% save data to never open this file again
% save("Dog_RosBag","data")

%% functions I stole from mathworks and modified
function [data] = getXYZFromPoint(msg,type,vars)
% Extract internal message first for efficiency
% (matters more for more nested messages)
data = struct();


for i = 1:length(msg)
    vec = msg{i}.(type);
    for j = 1:length(vars)
        data.(vars{j})(i) = vec.(vars{j});
    end
end

end

function [data] = getArrFromPoint(msg,type)
% Extract internal message first for efficiency
% (matters more for more nested messages)
% data = struct();


for i = 1:length(msg)
    vec = msg{i}.(type);
    data(:,i) = vec;
end

end