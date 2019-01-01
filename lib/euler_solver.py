#!/usr/bin/env python
# coding: utf-8
# Copyright (C) 2018  Nguyen Ngoc Sang, <https://github.com/SangVn> 

# Giải hệ phương trình Euler 1D sơ đồ Godunov

import numpy as np
import matplotlib.pyplot as plt
from decay import*
from constants import*
from convert_variables import*
from runge_kutta import*
from muscl import*
from weno5 import*

#tìm nghiệm chính xác bài toán Riemann 
def riemann_exact_solution(Pl, Pr, x, xstar, time_target):
    #mặc định tìm áp suất P bằng phương pháp lặp Godunov    
    pressure_solver = pressure_classic_godunov
    
    r1, u1, p1 = Pl[0], Pl[1], Pl[2]
    r2, u2, p2 = Pr[0], Pr[1], Pr[2]
    
    P, c1, c2, a1, a2 = pressure_solver(Pl, Pr)
    U = (a1*u1 + a2*u2 + p1 - p2)/(a1 + a2)
    #xét gián đoạn bên trái 
    if P > p1: #nếu là sóng xung kích 
        D1 = u1 - a1/r1
        D1star = D1
        R1 = r1*a1/(a1-r1*(u1-U))
    else: #nếu là sóng giãn
        D1 = u1 - c1
        c1star = c1 + gm1d2*(u1-U)
        D1star = U - c1star
        R1 = g*P/c1star**2
    #tương tự cho gián đoạn bên phải 
    if P > p2: #nếu là sóng xung kích 
        D2 = u2 + a2/r2
        D2star = D2
        R2 = r2*a2/(a2 + r2*(u2-U))
    else: #nếu là sóng giãn
        D2 = u2 + c2
        c2star = c2 - gm1d2*(u2-U)
        D2star = U + c2star
        R2 = g*P/c2star**2
        
    #Căn cứ vào vị trí của các gián đoạn tại thời điểm t và vị trí điểm x_i 
    w = (x - xstar)/time_target
    P_out = np.zeros((3, len(x)))
    for i, wi in enumerate(w):
        if wi<=D1:   #nằm bên trái sóng trái
            P_out[0, i] = r1
            P_out[1, i] = u1
            P_out[2, i] = p1
        elif wi>=D2: #nằm bên phải sóng phải
            P_out[0, i] = r2
            P_out[1, i] = u2
            P_out[2, i] = p2
        elif D1star<= wi <= U: #nằm giữa hai sóng 
            P_out[0, i] = R1
            P_out[1, i] = U
            P_out[2, i] = P
        elif U<= wi <= D2star: #nằm giữa hai sóng 
            P_out[0, i] = R2
            P_out[1, i] = U
            P_out[2, i] = P
        elif D1< wi < D1star: #nằm trong sóng giãn trái
            cstar = gm1dgp1*(u1 - wi) + c1/gp1d2
            P_out[1, i] = wi + cstar
            P_out[2, i] = p1*(cstar/c1)**g2dgm1           
            P_out[0, i] = g*P_out[2, i]/cstar**2
        elif D2star< wi < D2: #nằm trong sóng giãn phải 
            cstar = gm1dgp1*(wi - u2) + c2/gp1d2
            P_out[1, i] = wi - cstar
            P_out[2, i] = p2*(cstar/c2)**g2dgm1          
            P_out[0, i] = g*P_out[2, i]/cstar**2
        else:
            print 'Decay Error!'

    return P_out

#tính bước thời gian 
def time_step(P, CFL, dx):
    r, u, p = P[:, 0], P[:, 1], P[:, 2]
    c = np.sqrt(g*p/r) #vận tốc âm thanh 
    u_max = max(abs(u) + c)        
    dt = CFL*dx/u_max
    return dt


def euler_solver(Ps, reconstr, x, time_target, CFL):
    #xác định phương pháp runge-kutta 
    rk = np.zeros((3,3))
    runge_kutta_order = 1.
    if reconstr == godunov_reconstr:
        runge_kutta_order = 1
        runge_kutta_1(rk)
    elif reconstr == muscl_reconstr:
        runge_kutta_order = 2
        runge_kutta_2(rk)
    elif reconstr == weno5_reconstr:
        runge_kutta_order = 3
        runge_kutta_3(rk)
    else:
        print ('Error: reconstruction!')
        
    #mặc định tìm áp suất P bằng phương pháp lặp Godunov    
    pressure_solver = pressure_classic_godunov
    
    
    nx = len(x)
    dx = x[1] - x[0] #xét lưới đều 
    Us = Ps.copy()
    P2U(Ps, Us)
    Pstar = Ps[:-1].copy()
    Fstar = Ps[:-1].copy()
    time = 0.0
    
    while(time < time_target):
        Un = Us.copy()
        #tìm dt
        dt = time_step(Ps, CFL, dx) 
        if(time+dt > time_target): dt = time + dt - time_target
        time += dt
        #print time
        for stage in range(runge_kutta_order):
            #bước 1: tái cấu trúc - reconstruction 
            P_left, P_right = reconstr(Ps)
            #bước 2: tìm nghiệm phân rã gián đoạn, tính hàm dòng
            for i in range(nx-1): Pstar[i] = decay_godunov(P_left[i], P_right[i], pressure_solver)
            P2F(Pstar, Fstar)
                
            #bước 3: tích phân theo thời gian
            Us[1:-1] = rk[stage][0]*Un[1:-1] + rk[stage][1]*Us[1:-1] - rk[stage][2]*dt/dx*(Fstar[1:] - Fstar[:-1])
            
            U2P(Us, Ps) #tìm biến nguyên thủy
            
            #điều kiện biên: outflow
            Ps[0] = Ps[1]
            Ps[-1] = Ps[-2]        
    return  Ps.T

#hàm vẽ đồ thị so sánh 4 nghiệm: chính xác, godunov, muscl và weno5 
def plot_4P(x, Ps1, Ps2, Ps3, Ps4, img_name):
    e1 = Ps1[2]/(gm1*Ps1[0])
    e2 = Ps2[2]/(gm1*Ps2[0])
    e3 = Ps3[2]/(gm1*Ps3[0])
    e4 = Ps4[2]/(gm1*Ps4[0])
    
    f, axarr = plt.subplots(2, 2, figsize=(18,12))
    axarr[0, 0].plot(x, Ps1[0], x, Ps2[0], x, Ps3[0], x, Ps4[0])
    axarr[0, 0].set_title('density')
    axarr[0, 1].plot(x, Ps1[1], x, Ps2[1], x, Ps3[1], x, Ps4[1])
    axarr[0, 1].set_title('velocity')
    axarr[1, 0].plot(x, Ps1[2], x, Ps2[2], x, Ps3[2], x, Ps4[2])
    axarr[1, 0].set_title('pressure')
    axarr[1, 1].plot(x, e1, x, e2, x, e3, x, e4)
    axarr[1, 1].set_title('energy')
    f.subplots_adjust(hspace=0.2)
    plt.legend(['exact solution', 'godunov reconstr', 'muscl reconstr', 'weno5 reconstr'])
    plt.savefig(img_name+'.png')
    plt.show()
    
