clear all
close all
clc

% columns 1-3 accel [m/s^2]
% 4-6 gyro [rad/s]
% 7-9 mag [uT]
% [data, dt] = loadSoftsysData("hw3_data.txt");

[data, dt, gps_course, gps_course_time] = loadNavData("data/normal_vecNav.mat");

% dt = 0.01;
% [data, true_orientation] = simImuData(dt);

accel = data(:, 1:3);
% gyro = deg2rad(data(:, 4:6));
gyro = data(:, 4:6);
mag = data(:, 7:9);

% mag calibration pre-processing step
[cal_data, ~, ~, ~] = loadNavData("data/6_28 Data/6-28_Mag_Cal_vecNav.mat");
mag_cal_data = data(:, 7:9);
[A, b, expmfs] = magcal(mag_cal_data);

% produce corrected mag data
mag = (mag - b)*A;

% initial rotation matrix, normalize columns
C_init = [cross(cross(accel(1,:)', mag(1,:)'), accel(1,:)') cross(accel(1,:)',mag(1,:)') accel(1,:)'];
C_init(:, 1) = C_init(:,1)./norm(C_init(:,1));
C_init(:, 2) = C_init(:,2)./norm(C_init(:,2));
C_init(:, 3) = C_init(:,3)./norm(C_init(:,3));

gyro_int(:,1) = [deg2rad(-90); deg2rad(50); deg2rad(-50)];

q_hat(:, 1) = [(1/2)*sqrt(C_init(1,1)+C_init(2,2)+C_init(3,3)+1);
               (1/2)*sign(C_init(3,2)-C_init(2,3))*sqrt(C_init(1,1)-C_init(2,2)-C_init(3,3)+1);
               (1/2)*sign(C_init(1,3)-C_init(3,1))*sqrt(C_init(2,2)-C_init(3,3)-C_init(1,1)+1);
               (1/2)*sign(C_init(2,1)-C_init(1,2))*sqrt(C_init(3,3)-C_init(1,1)-C_init(2,2)+1)];
q_hat_gyro_only(:, 1) = q_hat(:,1);
eul_gyro_only(:, 1) = [0; 0; 0];
eul(:, 1) = [0; 0; 0];
t(1) = 0;

P = 10*eye(4);
sigma_gyro = 0.075; % assuming equal noise on each axis
Sigma_gyro = (sigma_gyro^2).*eye(3);

% Q1 = (sigma_gyro^2).*eye(4);
% Q = Q1;

% assuming equal noise on each axis
sigma_acc = 0.5; 
sigma_mag = 1.1;
R = [(sigma_acc^2).*eye(3) zeros(3);
     zeros(3) (sigma_mag^2).*eye(3)];

magnetic_dip_angle = deg2rad(-14.78); % deg, for auburn, al, from magnetic-declination.com
magnetic_field_vec = [cos(magnetic_dip_angle) 0 sin(magnetic_dip_angle)]'./sqrt(cos(magnetic_dip_angle)^2 + sin(magnetic_dip_angle)^2); % ned
grav_vec = [0 0 -1]'; % ned

% magnetic_field_vec = [0 cos(magnetic_dip_angle) -sin(magnetic_dip_angle)]'./sqrt(cos(magnetic_dip_angle)^2 + sin(magnetic_dip_angle)^2); % enu
% grav_vec = [0 0 1]'; % enu

for i = 1:length(accel(:, 1))-1
    t(i+1) = dt*i;

    % --- time update----
    Omega = [0  -gyro(i+1,:);
             gyro(i+1,:)' -vec3ToSkewSym(gyro(i+1,:))];

    % from stack overflow - integrating quaternion using angular rates
    q_w_a = cos(norm(gyro(i+1, :))*dt/2);
    q_w_v = sin(norm(gyro(i+1, :))*dt/2).*gyro(i+1,:)'/(norm(gyro(i+1,:)));
    q_w = [q_w_a; q_w_v];

    % state_prop = cos(norm(gyro(i+1, :))*dt/2)*eye(4) + (2/(norm(gyro(i+1,:))))...
    %                 *sin(norm(gyro(i+1,:))*dt/2)*Omega;

    gyro_int(:,i+1) = gyro_int(:, i) + gyro(i+1, :)'*dt;

    % q_hat(:, i+1) = state_prop*q_hat(:,i);
    q_hat(:, i+1) = quatmult(q_hat(:,i), q_w);

    % normalize quaternion
    q_hat(:, i+1) = q_hat(:, i+1)./norm(q_hat(:, i+1));

    % q_hat_gyro_only(:, i+1) = state_prop*q_hat_gyro_only(:, i);
    q_hat_gyro_only(:, i+1) = quatmult(q_hat_gyro_only(:,i), q_w);


    % normalize quaternion
    q_hat_gyro_only(:, i+1) = q_hat_gyro_only(:, i+1)./norm(q_hat_gyro_only(:, i+1));    
    eul_gyro_only(:, i+1) = rad2deg(unwrap(quat2eul(q_hat_gyro_only(:, i+1)')));


    A_d = [1 (-dt/2)*gyro(i+1,1) (-dt/2)*gyro(i+1,2) (-dt/2)*gyro(i+1,3);
           (dt/2)*gyro(i+1,1) 1 (dt/2)*gyro(i+1,3) (-dt/2)*gyro(i+1,2);
           (dt/2)*gyro(i+1,2) (-dt/2)*gyro(i+1,3) 1 (dt/2)*gyro(i+1,1);
           (dt/2)*gyro(i+1,3) (dt/2)*gyro(i+1,2) (-dt/2)*gyro(i+1,1) 1;];

    % Jacobian of f w.r.t. angular rates (omega)
    W = (dt/2).*[-q_hat(2, i+1)  -q_hat(3, i+1)  -q_hat(4, i+1);
                  q_hat(1, i+1)  -q_hat(4, i+1)   q_hat(3, i+1);
                  q_hat(4, i+1)   q_hat(1, i+1)  -q_hat(2, i+1);
                 -q_hat(3, i+1)   q_hat(2, i+1)  -q_hat(1, i+1)];

    Q = W*Sigma_gyro*W'; % Q changes over time depending on angular vel
    
    % covariance propagation
    P = A_d*P*A_d' + Q;

    % --- measurement correction ---

    % normalize the accel and mag so we can compare them with our gravity
    % and magnetic field unit vectors
    accel_normalized = accel(i+1,:)./(norm(accel(i+1,:)));
    mag_normalized = mag(i+1,:)./(norm(mag(i+1,:)));

    % rotate gravity and mag vectors into sensor's frame (this is the
    % nonlinear measurement model, h(x))
    % rot_mat = quatToRotMat(q_hat(:, i+1));
    rot_mat = quat2rotm(q_hat(:,i+1)');
    accel_hat = rot_mat'*grav_vec;
    mag_hat = rot_mat'*magnetic_field_vec;
    Y_hat = [accel_hat; mag_hat];

    % form jacobian of measurement model (H)
    u_g = cross(grav_vec, q_hat(2:4, i+1));
    u_r = cross(magnetic_field_vec, q_hat(2:4, i+1));
    skew_sym_first_row = vec3ToSkewSym(u_g+q_hat(1,i+1)*grav_vec);
    skew_sym_second_row = vec3ToSkewSym(u_r+q_hat(1,i+1)*magnetic_field_vec);
    H = 2.*[u_g skew_sym_first_row+dot(q_hat(2:4,i+1),grav_vec)*eye(3)-grav_vec*q_hat(2:4,i+1)';
            u_r skew_sym_second_row+dot(q_hat(2:4,i+1),magnetic_field_vec)*eye(3)-magnetic_field_vec*q_hat(2:4,i+1)';
            ];

    Y = [accel_normalized'; mag_normalized'];

    % compute Kalman gain
    v = Y - Y_hat;
    S = H*P*H' + R;
    K = P*H'*inv(S);

    % apply measurement correction to state and covariance
    corrected_state = q_hat(:, i+1) + K*v;
    q_hat(:, i+1) = corrected_state(1:4);
    P = (eye(4) - K*H)*P';

    % normalize quaternion
    q_hat(:, i+1) = q_hat(:, i+1)./norm(q_hat(:, i+1));

    eul(:, i+1) = rad2deg(unwrap(quat2eul(q_hat(:, i+1)')));
end


% figure
% plot(t, eul(1, :));
% hold on
% plot(t, eul_gyro_only(1,:));
% % plot(t(1:end-1), true_orientation(:,1));
% plot(t, rad2deg(gyro_int(3,:)));
% legend("EKF", "Gyro Only", "Gyro Only (Euler)");
% title("Rotation about North Axis");
% 
% figure
% plot(t, eul(2, :));
% hold on
% plot(t, eul_gyro_only(2,:));
% % plot(t(1:end-1), true_orientation(:,2));
% plot(t, rad2deg(gyro_int(2,:)));
% legend("EKF", "Gyro Only", "Gyro Only (Euler)");
% title("Rotation about East Axis");

figure
plot(t, eul(3, :));
hold on
% plot(t, eul_gyro_only(1,:));
% plot(t(1:end-1), true_orientation(:,3));
plot(gps_course_time, rad2deg(gps_course));
% plot(t, rad2deg(gyro_int(1,:)));
legend("EKF", "GPS Course");
title("Heading");
