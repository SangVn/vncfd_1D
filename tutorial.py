#!/usr/bin/env python
# coding: utf-8
# Copyright (C) 2018  Nguyen Ngoc Sang, <https://github.com/SangVn> 

# import thư viện 
import lib.muscl
from lib.euler_solver import*

#chia lưới
xl = 0    #giới hạn biên trái
xr = 1    #giới hạn biên phải
nx = 81   #số điểm lưới 
x, dx = np.linspace(xl, xr, nx, retstep=True)  #chia lưới

#điều kiện ban đầu 
Pl = [1.0, 0.0, 1000.0]
Pr = [1.0, 0.0, 0.01]
xstar = 0.5
t_max = 0.01

#tạo mảng chứa nghiệm 
Ps = np.zeros((nx, 3))
Ps[:int(xstar/dx)] = Pl
Ps[int(xstar/dx):] = Pr



P_exact = riemann_exact_solution(Pl, Pr, x, xstar, t_max)

P_godunov = euler_solver(Ps.copy(), godunov_reconstr, x, t_max, CFL=0.45)

lib.muscl.limiter = ospre
P_muscl = euler_solver(Ps.copy(), muscl_reconstr, x, t_max, CFL=0.45)

P_weno5 = euler_solver(Ps.copy(), weno5_reconstr, x, t_max, CFL=0.45)

plot_4P(x, P_exact, P_godunov, P_muscl, P_weno5, 'riemann')
