function q1_q2 = quatmult(q1,q2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    q1_q2_w = q1(1)*q2(1) - q1(2)*q2(2) - q1(3)*q2(3) - q1(4)*q2(4);
    q1_q2_x = q1(1)*q2(2) + q1(2)*q2(1) + q1(3)*q2(4) - q1(4)*q2(3);
    q1_q2_y = q1(1)*q2(3) - q1(2)*q2(4) + q1(3)*q2(1) + q1(4)*q2(2);
    q1_q2_z = q1(1)*q2(4) + q1(2)*q2(3) - q1(3)*q2(2) + q1(4)*q2(2);
    q1_q2 = [q1_q2_w; q1_q2_x; q1_q2_y; q1_q2_z];
end

