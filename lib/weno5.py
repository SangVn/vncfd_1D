#!/usr/bin/env python
# coding: utf-8
# Copyright (C) 2018  Nguyen Ngoc Sang, <https://github.com/SangVn> 

from numpy import array

#WENO schemes
epsilon = 1e-6

# 1. WENO5 Interpolation
#các hệ số trong công thức chỉ thị độ trơn 
a_b0 = [4., -19., 25., 11., -31., 10.]
a_b1 = [4., -13., 13., 5., -13., 4.]
a_b2 = [10., -31., 25., 11., -19., 4.]
a_bItp  = [a_b0, a_b1, a_b2]

#tái cấu trúc u_{i+1/2}
#trọng lượng tuyến tính 
gamma_Itp_p05 = array([1./16, 5./8, 5./16])

#các hệ số trong công thức xấp xỉ bậc ba 
a_u0 = [3./8, -5./4, 15./8]
a_u1 = [-1./8, 3./4, 3./8]
a_u2 = [3./8, 3./4, -1./8]
a_uItp_p05  = [a_u0, a_u1, a_u2]

#tái cấu trúc u_{i-1/2}
gamma_Itp_m05 = array([5./16, 5./8, 1./16])

#các hệ số trong công thức xấp xỉ bậc ba 
a_u0 = [-1./8, 3./4, 3./8]
a_u1 = [3./8, 3./4, -1./8]
a_u2 = [15./8, -5./4, 3./8]
a_uItp_m05  = [a_u0, a_u1, a_u2]


# 2. WENO5 Reconstruction
a_b0 = [1., -4., 3.]
a_b1 = [1., 0., -1.]
a_b2 = [3., -4., 1.]
a_bRct = [a_b0, a_b1, a_b2]

#reconstruction u_{i+1/2}
gamma_Rct_p05 = array([1./10, 3./5, 3./10])

a_u0 = [1./3, -7./6, 11./6]
a_u1 = [-1./6, 5./6, 1./3]
a_u2 = [1./3, 5./6, -1./6]
a_uRct_p05  = [a_u0, a_u1, a_u2]

#reconstruction u_{i-1/2}
gamma_Rct_m05 = array([3./10, 3./5, 1./10])

a_u0 = [-1./6, 5./6, 1./3]
a_u1 = [1./3, 5./6, -1./6]
a_u2 = [11./6, -7./6, 1./3]
a_uRct_m05  = [a_u0, a_u1, a_u2]

#WENO5 Interpolation: smothness indicators
def beta_Itp(a, u):
    b = 1./3*(a[0]*u[0]**2 + a[1]*u[0]*u[1] + a[2]*u[1]**2 + a[3]*u[0]*u[2] + a[4]*u[1]*u[2] + a[5]*u[2]**2)
    return b

#WENO5 Reconstruction: smothness indicators
def beta_Rct(a, u):
    b = 13./12*(u[0] - 2*u[1] + u[2])**2 + 0.25*(a[0]*u[0] + a[1]*u[1] + a[2]*u[2])**2
    return b

#Hàm tìm trọng lượng phi tuyến 
def omega(gamma, beta):
    omgs = gamma/(epsilon + beta)**2
    omgs_sum = sum(omgs)
    omg = omgs/omgs_sum
    return omg

def psum(a, P3):
    return a[0]*P3[0] + a[1]*P3[1] + a[2]*P3[2]

#Hàm tái cấu trúc nghiệm
def reconstr(Ps, beta, a_b, gamma, a_P, Ps_05):
    for i in range(2, len(Ps)-2): #for each cell in cells
        #bước 1. tính chỉ thị độ trơn trong từng khuôn S_j
        b0 = beta(a_b[0], Ps[i-2:i+1])
        b1 = beta(a_b[1], Ps[i-1:i+2])
        b2 = beta(a_b[2], Ps[i  :i+3])
        
        #bước 2. tìm trọng lượng phi tuyến từng khuôn
        omg_r = omega(gamma, array([b0[0], b1[0], b2[0]]))
        omg_u = omega(gamma, array([b0[1], b1[1], b2[1]]))
        omg_p = omega(gamma, array([b0[2], b1[2], b2[2]]))
    
        #bước 3. tìm xấp xỉ bậc 3
        P0 = psum(a_P[0], Ps[i-2:i+1])
        P1 = psum(a_P[1], Ps[i-1:i+2])
        P2 = psum(a_P[2], Ps[i  :i+3])
        
        #bước 4. tìm xấp xỉ bậc 5
        Ps_05[i, 0] = sum(omg_r*[P0[0], P1[0], P2[0]])
        Ps_05[i, 1] = sum(omg_u*[P0[1], P1[1], P2[1]])
        Ps_05[i, 2] = sum(omg_p*[P0[2], P1[2], P2[2]])
        
        #Điều kiện biên: u_05[0], u_05[1], u_05[-2], u_05[-1]
    return


#sơ đồ WENO5
def weno5_reconstr(Ps):
    P_left = Ps.copy()
    P_right = Ps.copy()
    reconstr(Ps, beta_Rct, a_bRct, gamma_Rct_p05, a_uRct_p05, P_left)
    reconstr(Ps, beta_Rct, a_bRct, gamma_Rct_m05, a_uRct_m05, P_right)
    
    return P_left[:-1], P_right[1:]

#sơ đồ WENO5
def weno5_interpolation_reconstr(Ps):
    P_left = Ps.copy()
    P_right = Ps.copy()
    reconstr(Ps, beta_Itp, a_bItp, gamma_Itp_p05, a_uItp_p05, P_left)
    reconstr(Ps, beta_Itp, a_bItp, gamma_Itp_m05, a_uItp_m05, P_right)
    
    return P_left[:-1], P_right[1:]
