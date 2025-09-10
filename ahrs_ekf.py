import numpy as np

class AhrsEKF:

    mag_cal_A = 0 # TODO
    mag_cal_b = 0

    def __init__(self, A, b):
        self.mag_cal_A = A
        self.mag_cal_b = b

    def initialize(self, accel_init, mag_init):
        # create initial rotation matrix
        C1 = np.cross(np.cross(accel_init, mag_init), accel_init)
        C1 = C1/np.linalg.norm(C1)
        C2 = np.cross(accel_init, mag_init)
        C2 = C2/np.linalg.norm(C2)
        C3 = accel_init
        C3 = C3/np.linalg.norm(C3)
        self.C_init = np.column_stack((C1, C2, C3))

        q_hat1 = (1/2)*np.sqrt(self.C_init[0][0]+self.C_init[1][1]+self.C_init[2,2]+1)
        q_hat2 = (1/2)*np.sign(self.C_init[2][1]-self.C_init[1][2])*np.sqrt(self.C_init[0][0]-self.C_init[1][1]-self.C_init[2][2]+1)
        q_hat3 = (1/2)*np.sign(self.C_init[0][2]-self.C_init[2][0])*np.sqrt(self.C_init[1][1]-self.C_init[2][2]-self.C_init[0][0]+1)
        q_hat4 = (1/2)*np.sign(self.C_init[1][0]-self.C_init[0][1])*np.sqrt(self.C_init[2][2] - self.C_init[0][0]-self.C_init[1][1]+1)
        self.q_hat = np.vstack((q_hat1, q_hat2, q_hat3, q_hat4))
        
        self.P = 10*np.identity(4)

        # set tunings
        sigma_gyro = 0.075
        self.Sigma_gyro = (sigma_gyro**2)*np.identity(3)

        sigma_acc = 0.5 # TODO: change tunings
        sigma_mag = 1.1
        self.R = np.identity(6)
        self.R[0][0] = (sigma_acc**2)
        self.R[1][1] = (sigma_acc**2)
        self.R[2][2] = (sigma_acc**2)
        self.R[3][3] = (sigma_mag**2)
        self.R[4][4] = (sigma_mag**2)
        self.R[5][5] = (sigma_mag**2)

        magnetic_dip_angle = np.deg2rad(-14.78)
        self.magnetic_field_vec = np.array([np.cos(magnetic_dip_angle), 0, np.sin(magnetic_dip_angle)]).T/np.sqrt(np.cos(magnetic_dip_angle)**2 + np.sin(magnetic_dip_angle)**2)
        self.grav_vec = np.array([0, 0, -1]).T

        print("AHRS EKF Initialized")

    def time_update(self, gyro, dt):
        omega_skew_sym = vec3ToSkewSym(gyro)
        Omega_row1 = np.hstack(([0], -gyro.T))
        Omega_block2 = np.column_stack((gyro.T, vec3ToSkewSym(gyro)))
        Omega = np.vstack((Omega_row1, Omega_block2))

        q_w_a = np.cos(np.linalg.norm(gyro)*dt/2)
        q_w_v = np.sin(np.linalg.norm(gyro)*dt/2)*gyro/np.linalg.norm(gyro)
        q_w = np.hstack((q_w_a, q_w_v)).T
        self.q_hat = quatmult(self.q_hat, q_w)
        self.q_hat = self.q_hat[0]

        self.A_d = np.array([[1, (-dt/2)*gyro[0], (-dt/2)*gyro[1], (-dt/2)*gyro[2]],
                             [(dt/2)*gyro[0], 1, (dt/2)*gyro[2], (-dt/2)*gyro[1]],
                             [(dt/2)*gyro[1], (-dt/2)*gyro[2], 1, (dt/2)*gyro[0]],
                             [(dt/2)*gyro[2], (dt/2)*gyro[1], (-dt/2)*gyro[0], 1]])
        # Jacobian of f w.r.t. angular rates (omega)
        W = (dt/2)*np.array([[-self.q_hat[1],  -self.q_hat[2],  -self.q_hat[3]],
                             [self.q_hat[0],   -self.q_hat[3],   self.q_hat[2]],
                             [self.q_hat[3],    self.q_hat[0],  -self.q_hat[1]],
                             [-self.q_hat[2],   self.q_hat[1],  -self.q_hat[0]]])
        Q = np.matmul(np.matmul(W,self.Sigma_gyro),W.T)

        P = np.matmul(np.matmul(self.A_d, self.P), self.A_d.T) + Q

    def measurement_update(self, accel, mag):
        accel_normalized = accel/np.linalg.norm(accel)
        mag_normalized = mag/np.linalg.norm(mag)

        rot_mat = quat2rot(self.q_hat)
        accel_hat = np.matmul(rot_mat.T,self.grav_vec.T)
        mag_hat = np.matmul(rot_mat.T,self.magnetic_field_vec.T)
        Y_hat = np.hstack((accel_hat, mag_hat))

        u_g = np.cross(self.grav_vec, self.q_hat[1:])
        u_r = np.cross(self.magnetic_field_vec, self.q_hat[1:])
        skew_sym_first_row = vec3ToSkewSym(u_g + self.q_hat[0]*self.grav_vec)
        skew_sym_second_row = vec3ToSkewSym(u_r + self.q_hat[0]*self.magnetic_field_vec)
        H_12_block = skew_sym_first_row + np.dot(self.q_hat[1:], self.grav_vec)*np.identity(3) - np.matmul(self.grav_vec, self.q_hat[1:].T)
        H_22_block = skew_sym_second_row +np.dot(self.q_hat[1:], self.magnetic_field_vec)*np.identity(3) - np.matmul(self.magnetic_field_vec, self.q_hat[1:].T)
        H_row1_block = np.hstack((u_g[:, None], H_12_block))
        H_row2_block = np.hstack((u_r[:, None], H_22_block))
        H = 2*np.vstack((H_row1_block, H_row2_block))

        Y = np.hstack((accel_normalized, mag_normalized))

        # compute Kalman gain
        v = Y - Y_hat
        S = np.matmul(np.matmul(H, self.P), H.T) + self.R
        K = np.matmul(np.matmul(self.P, H.T), np.linalg.inv(S))

        self.q_hat = self.q_hat + np.matmul(K, v)
        self.P = np.matmul((np.identity(4) - np.matmul(K, H)), self.P.T)

        print(self.q_hat)
        print(self.P)

def vec3ToSkewSym(q):
    skew_sym_mat = np.array([[0, -q[2], q[1]],
                             [q[2], 0, -q[0]],
                             [-q[1], q[0], 0]])
    return skew_sym_mat

def quatmult(q0,q1):
   # quaternion multiplication
    
   q0_q1_w = q0[0]*q1[0] - q0[1]*q1[1] - q0[2]*q1[2] - q0[3]*q1[3]
   q0_q1_x = q0[0]*q1[1] + q0[1]*q1[0] + q0[2]*q1[3] - q0[3]*q1[2]
   q0_q1_y = q0[0]*q1[2] - q0[1]*q1[3] + q0[2]*q1[0] + q0[3]*q1[1]
   q0_q1_z = q0[0]*q1[3] + q0[1]*q1[2] - q0[2]*q1[1] + q0[3]*q1[1]
   
   q0_q1 = np.array([q0_q1_w, q0_q1_x, q0_q1_y, q0_q1_z]).T;    
   
   return q0_q1


# converted from Howard Chen's matlab function
def quat2rot(q):
   #  quatToDCM converts from quaternion vector to DCM
   #  input: quaternion vector [nx3] first element is real
   #  output: Rotation matrix [2x2xn]
   #
   #  Author: not Howard Chen


   R = np.zeros((3,3))


   R[0, 0] = q[0]**2+q[1]**2-q[2]**2-q[3]**2
   R[0, 1] = -2*q[0]*q[3]+2*q[1]*q[2]
   R[0, 2] = 2*q[0]*q[2]+2*q[1]*q[3]
   R[1, 0] = 2*q[0]*q[3]+2*q[1]*q[2]
   R[1, 1] = q[0]**2-q[1]**2+q[2]**2-q[3]**2
   R[1, 2] = -2*q[0]*q[1]+2*q[2]*q[3]
   R[2, 0] = -2*q[0]*q[2]+2*q[1]*q[3]
   R[2, 1] = 2*q[0]*q[1]+2*q[2]*q[3]
   R[2, 2] = q[0]**2-q[1]**2-q[2]**2+q[3]**2

   return R


# A = np.identity(3)
# b = np.array([-0.2061, -0.0178, -0.0520])
# a = AhrsEKF(A, b)
# mag_init = np.array([0.0542322857859514, -0.215706835372906, -0.0152498115500106])
# accel_init = np.array([9.31001853942871,	0.113502115011215,	-2.58657431602478])
# gyro_test = np.array([-0.0567961037158966,	0.0142987715080380,	0.0285795908421278])
# dt = 0.1
# a.initialize(accel_init, mag_init)
# a.time_update(gyro_test, dt)
# a.measurement_update(accel_init, mag_init)