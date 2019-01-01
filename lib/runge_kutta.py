#!/usr/bin/env python
# coding: utf-8
# Copyright (C) 2018  Nguyen Ngoc Sang, <https://github.com/SangVn> 

#các hệ số trong runge-kutta
#runge-kutta bậc 1 (euler)
def runge_kutta_1(rk):
    rk[0][0], rk[0][1], rk[0][2] = 1.0, 0.0, 1.0
    
#runge-kutta bậc 2 
def runge_kutta_2(rk):
    rk[0][0], rk[0][1], rk[0][2] = 1.0, 0.0, 1.0
    rk[1][0], rk[1][1], rk[1][2] = 0.5, 0.5, 0.5

#runge-kutta bậc 3 
def runge_kutta_3(rk):
    rk[0][0], rk[0][1], rk[0][2] = 1.0, 0.0, 1.0
    rk[1][0], rk[1][1], rk[1][2] = 0.75, 0.25, 0.25
    rk[2][0], rk[2][1], rk[2][2] = 1./3, 2./3, 2./3
