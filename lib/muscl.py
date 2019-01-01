#!/usr/bin/env python
# coding: utf-8
# Copyright (C) 2018  Nguyen Ngoc Sang, <https://github.com/SangVn> 

from numpy import array

#MUSCL schemes

#P=(r, u, p): amax(a, P) = (max(a, r), max(a, u), max(a, p))
def amax(a, P):
    return array([max(a, pi) for pi in P])

def amin(a, P):
    return array([min(a, pi) for pi in P])

#slope limiter: phi(r), rs - mảng r
def minmod(rs):
    phi = array([amax(0, amin(1, r)) for r in rs])
    return phi

def koren(rs):
    phi = array([amax(0, amin(2*r, amin((1+2*r)/3., 2))) for r in rs])
    return phi

def superbee(rs):
    phi = array([amax(0, amin(2*r, 1), amin(r,2)) for r in rs])
    return phi

def vanleer74(rs):
    phi = array([(r+abs(r))/(1.+abs(r)) for r in rs])
    return phi

def vanleer77(rs):
    phi = array([amax(0, amin(2*r, 0.5*(1+r),2)) for r in rs])
    return phi

def ospre(rs):
    phi = array([1.5*(r+r*r)/(1+r+r*r) for r in rs])
    return phi

limiter = minmod

#muscl reconstruction
def muscl_reconstr(Ps):
    P_left = Ps[:-1].copy()
    P_right = Ps[1:].copy()
    rs = (Ps[1:-1] - Ps[:-2])/(Ps[2:] - Ps[1:-1] + 1e-10)
    
    P_left[1:] = Ps[1:-1] + 0.5*limiter(rs[:])*(Ps[2:] - Ps[1:-1])
    P_left[0] = Ps[0] #điều kiện biên 
    
    P_right[:-1] = Ps[1:-1] - 0.5*limiter(rs[:])*(Ps[2:] - Ps[1:-1])
    P_right[-1] = Ps[-1] #điều kiện biên
    
    return P_left, P_right

#tái cấu trúc nghiệm hằng số từng mảnh - Godunov
def godunov_reconstr(Ps):
    P_left = Ps[:-1].copy()
    P_right = Ps[1:].copy()
    return P_left, P_right
