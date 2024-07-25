import numpy as np
import math

def gamma_pi0_nu(G_F, f, m_pi0, m_HNL):
	if m_HNL < m_pi0: return 0
	mass_ratio = np.power(m_pi0/m_HNL,2)
	factor = np.power(G_F,2)*np.power(f,2)*np.power(m_HNL,3)/(32*math.pi)
	return factor*np.power(1-mass_ratio,2)

def gamma_H_l(G_F, V, f, m_H, m_l, m_HNL):
	if m_HNL < (m_H+m_l): return 0
	mass_ratio1 = np.power(m_l/m_HNL,2)
	mass_ratio2 = np.power(m_H/m_HNL,2)
	mass_diff_minus = np.power((m_H-m_l)/m_HNL,2)
	mass_diff_plus = np.power((m_H+m_l)/m_HNL,2)
	factor = np.power(G_F*V*f,2)*np.power(m_HNL,3)/(16*math.pi)
	return factor*(np.power(1-mass_ratio1,2)-(mass_ratio2*(1+mass_ratio1)))*np.sqrt((1-mass_diff_minus)*(1-mass_diff_minus))

def gamma_rho_l(g_p, G_F, V, m_rho, m_l, m_HNL):
	if m_HNL < (m_rho+m_l): return 0
	mass_ratio1 = np.power(m_l/m_HNL,2)
	mass_ratio2 = np.power(m_rho/m_HNL,2)
	mass_diff_minus = np.power((m_rho-m_l)/m_HNL,2)
	mass_diff_plus = np.power((m_rho+m_l)/m_HNL,2)
	mass_diff_minus_2 = (np.power(m_l,2)-2*np.power(m_rho,2))/np.power(m_HNL,2)
	factor = np.power(g_p/m_rho,2)*np.power(G_F,2)*np.power(V,2)*np.power(m_HNL,3)/(8*math.pi)
	return factor*(np.power(1-mass_ratio1,2)-mass_ratio2*(1+mass_diff_minus_2))*np.sqrt((1-mass_diff_minus)*(1-mass_diff_minus))

def gamma_rho0_nu(g_p, G_F, m_rho0, m_HNL):
	if m_HNL < m_rho0: return 0
	mass_ratio = np.power(m_rho0/m_HNL,2)
	factor = np.power(g_p/m_rho0,2)*np.power(G_F,2)*np.power(m_HNL,3)/(16*math.pi)
	return factor*(1+2*mass_ratio)*np.power(1-mass_ratio,2)

def gamma_nu1_l1_l2(G_F, m_l1, m_l2, m_HNL):
	if m_HNL < (m_l1+m_l2): return 0
	factor = np.power(G_F,2)*np.power(m_HNL,5)/(192*np.power(math.pi,3))
	mass_ratio = max(m_l1,m_l2)/m_HNL
	x_factor = 1-(8*np.power(mass_ratio,2))+(8*np.power(mass_ratio,6))-np.power(mass_ratio,8)-(12*np.power(mass_ratio,4)*np.log(np.power(mass_ratio,2)))
	return factor*x_factor

def gamma_nu_ll(G_F, sin2thetaW, m_l, m_HNL):
	if m_HNL < 2*m_l: return 0
	x = np.power(m_l/m_HNL,2)
	L = 0
	num = 0
	den = 0
	if x <= 0.25:
		num = 1-(3*x)-(1-x)*np.sqrt(1-4*x)
		den = x*(1+np.sqrt(1-4*x))
	if den != 0 and num/den > 0: L = np.log10(num/den)
	C1 = 0.25*(1-4*sin2thetaW+8*np.power(sin2thetaW,2))
	C2 = 0.5*sin2thetaW*(2*sin2thetaW-1)
	C3 = 0.25*(1+4*sin2thetaW+8*np.power(sin2thetaW,2))
	C4 = 0.5*sin2thetaW*(2*sin2thetaW+1)
	factor = np.power(G_F,2)*np.power(m_HNL,5)/(192*np.power(math.pi,3))
	gamma_same_flavor = (C3*((1-14*x-2*np.power(x,2)-12*np.power(x,3))*np.sqrt(1-4*x)+12*np.power(x,2)*(np.power(x,2)-1)*L)+(4*C4*(x*(2+10*x-12*np.power(x,2))*np.sqrt(1-4*x)+6*np.power(x,2)*(1-2*x+2*np.power(x,2))*L)))
	gamma_opposite_flavor = (C1*((1-14*x-2*np.power(x,2)-12*np.power(x,3))*np.sqrt(1-4*x)+12*np.power(x,2)*(np.power(x,2)-1)*L)+(4*C2*(x*(2+10*x-12*np.power(x,2))*np.sqrt(1-4*x)+6*np.power(x,2)*(1-2*x+2*np.power(x,2))*L)))
	return factor*gamma_opposite_flavor

def FullWidth(m_HNL):
    ## Declare global physical constants
    G_F = 1.17E-05 ## Fermi constant in GeV^-2
    sin2thetaW = 0.23122 ## sin^2(weak angle)

    return gamma_nu_ll(G_F,sin2thetaW,511E-06,m_HNL)+gamma_pi0_nu(G_F,0.130,0.135,m_HNL)+gamma_nu_ll(G_F,sin2thetaW,0.1056,m_HNL)+gamma_pi0_nu(G_F,0.156,0.548,m_HNL)+gamma_rho0_nu(0.102,G_F,0.775,m_HNL)+gamma_pi0_nu(G_F,-0.0585,0.958,m_HNL)+gamma_nu1_l1_l2(G_F,511E-06,1.776,m_HNL)+gamma_nu1_l1_l2(G_F,0.1056,1.776,m_HNL)+gamma_H_l(G_F,0.97427,0.130,0.140,1.776,m_HNL)+gamma_H_l(G_F,0.224,0.130,0.494,1.776,m_HNL)+gamma_rho_l(0.102,G_F,9.74E-01,0.775,1.776,m_HNL)
