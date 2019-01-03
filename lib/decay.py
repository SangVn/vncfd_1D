#!/usr/bin/env python
# coding: utf-8
# Copyright (C) 2018  Nguyen Ngoc Sang, <https://github.com/SangVn> 

from scipy.optimize import fsolve
from constants import *

#phương pháp lặp tìm P
def pressure_classic_godunov(Pl, Pr):
    r1, u1, p1 = Pl[0], Pl[1], Pl[2]
    r2, u2, p2 = Pr[0], Pr[1], Pr[2]
    
    #vận tốc âm thanh 
    c1 = (g*p1/r1)**0.5
    c2 = (g*p2/r2)**0.5
    
    #phương pháp lặp 
    P0 = (p1*r2*c2 + p2*r1*c1 + (u1-u2)*r1*c1*r2*c2)/(r1*c1+r2*c2)
    if P0 < eps: P0 = eps
    iterations = 50 # max_iteration 
    while (True):
        P = P0 #áp suất P^{n-1}
        if P >= p1: a1 = (r1*(gp1d2*P + gm1d2*p1))**0.5
        else:
            pp = max(eps, P/p1)
            op = 1. - pp**gm1d2g
            if op>=eps: a1 = gm1d2g*r1*c1*(1. - pp)/op
            else: a1 = r1*c1
        if P >= p2: a2 = (r2*(gp1d2*P + gm1d2*p2))**0.5
        else:
            pp = max(eps, P/p2)
            op = 1. - pp**gm1d2g
            if op>=eps: a2 = gm1d2g*r2*c2*(1. - pp)/op
            else: a2 = r2*c2

        z = P/(p1+p2)
        alpha = gm1/(3*g)*(1. - z)/(z**gp1d2g)/(1. - z**gm1d2g) - 1.
        if alpha < 0.: alpha = 0.
        phi = (a2*p1 + a1*p2 + a1*a2*(u1 - u2))/(a1+a2)

        P0 = (alpha*P + phi)/(1. + alpha)#tính P^n
        iterations -= 1
        if (abs(P0 - P) < eps) or (not iterations): break
    #kết thúc vòng lặp! 
  
    return P, c1, c2, a1, a2

#giải phương trình f(P) = 0 , tìm P
def pressure_root_finding(Pl, Pr):
    r1, u1, p1 = Pl[0], Pl[1], Pl[2]
    r2, u2, p2 = Pr[0], Pr[1], Pr[2]
    
    #vận tốc âm thanh 
    c1 = (g*p1/r1)**0.5
    c2 = (g*p2/r2)**0.5

    fl_shock = lambda P: (P - p1)/(r1*c1*(gp1d2g*(P/p1) + gm1d2g)**0.5)
    fl_expansion = lambda P: c1/gm1d2*((P/p1)**gm1d2g - 1.)
    fr_shock = lambda P: (P - p2)/(r2*c2*(gp1d2g*(P/p2) + gm1d2g)**0.5)
    fr_expansion = lambda P: c2/gm1d2*((P/p2)**gm1d2g - 1.)
    
    def fl(P):
        if P >= p1: return fl_shock(P)
        else: return fl_expansion(P)
    def fr(P):
        if P >= p2: return fr_shock(P)
        else: return fr_expansion(P)
    
    f = lambda P: fl(P) + fr(P) + u2 - u1
    
    P,info, ier, msg = fsolve(f, (p1+p2)/2.,full_output=True,xtol=1.e-14)
    #trường hợp sóng giãn mạnh hay chân không 
    if ier!=1:
        P,info, ier, msg = fsolve(f, (p1+p2)/2.,full_output=True,factor=0.1,xtol=1.e-10)
        # trường hợp nghiệm không hội tụ 
        if ier!=1: 
            print ('Warning: fsolve did not converge.')
            print (msg)
   
    if P >= p1-eps: a1 = (r1*(gp1d2*P + gm1d2*p1))**0.5
    else: a1 =  gm1d2g*r1*c1*(1. - P/p1)/(1. - (P/p1)**gm1d2g)
    if P >= p2-eps: a2 = (r2*(gp1d2*P + gm1d2*p2))**0.5
    else: a2 =  gm1d2g*r2*c2*(1. - P/p2)/(1. - (P/p2)**gm1d2g)
    
    return P, c1, c2, a1, a2

#Giải bài toán phân rã gián đoạn trong trên từng bề mặt thể tích hữu hạn
#phân rã gián đoạn Godunov
def decay_godunov(Pl, Pr, pressure_solver):
    r1, u1, p1 = Pl[0], Pl[1], Pl[2]
    r2, u2, p2 = Pr[0], Pr[1], Pr[2]
    
    #tìm P, U, R
    P, c1, c2, a1, a2 = pressure_solver(Pl, Pr)
    U = (a1*u1 + a2*u2 + p1 - p2)/(a1 + a2)
    #xét gián đoạn bên trái 
    if P > p1: #nếu là sóng xung kích 
        D1 = u1 - a1/r1
        R1 = r1*a1/(a1-r1*(u1-U))
    else: #nếu là sóng giãn
        D1 = u1 - c1
        c1star = c1 + gm1d2*(u1-U)
        D1star = U - c1star
        R1 = g*P/c1star**2
    #tương tự cho gián đoạn bên phải 
    if P > p2: #nếu là sóng xung kích 
        D2 = u2 + a2/r2
        R2 = r2*a2/(a2 + r2*(u2-U))
    else: #nếu là sóng giãn
        D2 = u2 + c2
        c2star = c2 - gm1d2*(u2-U)
        D2star = U + c2star
        R2 = g*P/c2star**2
        
    #xét cấu hình phân rã xác định nghiệm PStar = (Rstar, Ustar, Pstar)
    #tùy theo vị trí biên i+1/2 nằm trong vùng nào (xem bài 12)

    if D1>0 and D2>0:   #nằm bên trái sóng trái
        Rstar = r1
        Ustar = u1
        Pstar = p1
    elif D1<0 and D2<0: #nằm bên phải sóng phải
        Rstar = r2
        Ustar = u2
        Pstar = p2
    elif D1<0 and D2>0: #nằm giữa hai sóng 
        if U>=0: Rstar = R1 #nằm bên trái gián đoạn tiếp xúc
        else:    Rstar = R2 #nằm bên phải gián đoạn tiếp xúc
        Ustar = U
        Pstar = P
    elif D1<0 and D1star>0: #nằm trong sóng giãn trái
        Ustar = gm1dgp1*u1 + c1/gp1d2
        Pstar = p1*(Ustar/c1)**g2dgm1
        Rstar = g*p1/Ustar**2   
    elif D2>0 and D2star<0: #nằm trong sóng giãn phải 
        Ustar = gm1dgp1*u2 - c2/gp1d2
        Pstar = p2*(Ustar/c2)**g2dgm1
        Rstar = g*p2/Ustar**2
    else:
        print ('Godunov _ Decay Error!')
        
    #vector biến gốc PStar
    PStar = [Rstar, Ustar, Pstar]
    return PStar
