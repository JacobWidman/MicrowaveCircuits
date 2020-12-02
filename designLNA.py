import math as m
import cmath as cm
import numpy as np

##algorimic steps (1)
def is_s12_unilateral(s11a,s12a,s21a,s22a): ##Unilateral calculations
    u = (s11a*s12a*s21a*s22a)/((1-m.pow(s11a,2))*(1-m.pow(s22a,2)))
    print(f'U is {u}')
    u_lower = (1/(m.pow((1+u),2)))
    u_lowerLog = 10*m.log10(u_lower)
    u_upper = (1/(m.pow((1-u),2)))
    u_upperLog = 10*m.log10(u_upper)
    print(f'Unilateral Lower {u_lower} Unilateral upper {u_upper}')
    print(f'Unilateral Lower dB {u_lowerLog} Unilateral upper dB {u_upperLog}')
    if abs(u_upperLog - u_lowerLog)> 1:
        return False
    else:
        return True

def get_stability(s11,s12,s21,s22,s11a,s12a,s21a,s22a): ## Checks stability
    
    delta_d = (s11*s22) - (s21*s12)
    d_abs = abs(delta_d)
    print(f'deltad {delta_d}, dabs {d_abs}')
    k = (1+m.pow(d_abs,2)-pow(s11a,2)-pow(s22a,2))/(2*s21a*s12a)
    print(f'k is {k}')
    if k > 1 and abs(delta_d) <1:
        print("unconditionally stable from K") ##Unconditionally stable, no special smith chart fun
        return delta_d, d_abs
    
    elif k < 1 and abs(delta_d) <1:
        print("not unconditionally stable from K") ##breaking out smith chart to plot if the s11 and s22 are possible for this circuit
        center_load, radius_load, center_source, radius_source = get_con_stablity_circles(s11,s12,s21,s22,s11a,s12a,s21a,s22a, d_abs, delta_d)
        center_sourcem = abs(center_source)
        center_loadm = abs(center_load)
        if center_loadm > 1:
            total = abs(radius_load - center_loadm)
            print(f'total load {total}')
        if center_sourcem > 1:
            center_sourcem = abs(center_source)
            total = abs(radius_source - center_sourcem)
            print(f'total source {total}')
        return delta_d, d_abs, center_load, radius_load, center_source, radius_source
    
    elif k > 1 and abs(delta_d) >1:
        print("not unconditionally stable from delta")
        center_load, radius_load, center_source, radius_source = get_con_stablity_circles(s11,s12,s21,s22,s11a,s12a,s21a,s22a, d_abs, delta_d)
        if center_load > 1:
            center_loadm = abs(center_load)
            total = radius_load - center_loadm
            print(f'total load {total}')
        elif center_source > 1:
            center_sourcem = abs(center_source)
            total = radius_source - center_source
            print(f'total source {total}')
        return delta_d, d_abs
    
    else:
        print("I had hoped that we would never go to this place. Something went wrong in the stability")
    return None

def get_con_stablity_circles(s11,s12,s21,s22,s11a,s12a,s21a,s22a, d_abs, delta_d): ## draw the center and the radius for the smith chart
    center_load = np.conjugate(s22-(delta_d*np.conjugate(s11)))/(m.pow(s22a,2)-m.pow(d_abs,2))
    radius_load = abs(abs(s12*s21)/(m.pow(s22a,2)-m.pow(d_abs,2)))

    center_source = np.conjugate(s11-(delta_d*np.conjugate(s22)))/(m.pow(s11a,2)-m.pow(d_abs,2))
    radius_source = abs(abs(s12*s21)/(m.pow(s11a,2)-m.pow(d_abs,2)))

    center_load1 = cm.polar(center_load)
    center_source1 = cm.polar(center_source)
    c_l_degree = np.degrees(center_load1[1])
    c_s_degree = np.degrees(center_source1[1])
    print(f"center load {center_load1}, degrees {c_l_degree}, \nradius load {radius_load}, \ncenter source {center_source1}, degrees {c_s_degree}, \nradius source {radius_source}")
    return center_load, radius_load, center_source, radius_source

def get_max_transducer_gain(s11,s12,s21,s22,s11a,s12a,s21a,s22a, delta_d, d_abs): ###get max transducer gain and calculate the gamma s and gamma l
    print(f'\nMAX GAIN\n')
    b1 = 1 + m.pow(s11a,2) - m.pow(s22a,2) - m.pow(d_abs,2)
    b2 = 1 + m.pow(s22a,2) - m.pow(s11a,2) - m.pow(d_abs,2)
    c1 = s11 - (delta_d*np.conjugate(s22))
    c2 = s22 - (delta_d*np.conjugate(s11))
    
    print(f'b1 {b1}, b2 {b2}, c1 {c1}, c2 {c2}')
    
    discr_s = m.pow(b1,2)-(4*m.pow(abs(c1),2))
    discr_l = m.pow(b2,2)-(4*m.pow(abs(c2),2))
    try:
        if discr_s > 0:
            gam_s_sqrt = m.sqrt(discr_s)
        elif discr_s < 0:
            gam_s_sqrt = m.sqrt(-m.pow(b1,2)+(4*m.pow(abs(c1),2)))
        else:
            gam_s_sqrt = 0
            
        if discr_l > 0:
            gam_l_sqrt = m.sqrt(discr_l)
        elif discr_s < 0:
            gam_l_sqrt = m.sqrt(-m.pow(b2,2)+(4*m.pow(abs(c2),2)))
        else:
            gam_l_sqrt = 0
    except ValueError:
        print("Value Error, discriminate is negative")
        
    print(f'dis s{discr_s}, dis l {discr_l}')
    gamma_s_lower = (b1+gam_s_sqrt)/(2*c1)
    
    gamma_s_higher = (b1-gam_s_sqrt)/(2*c1)
    if abs(gamma_s_lower) < 1 and abs(gamma_s_higher) > 1:
        gamma_s = gamma_s_lower
    elif abs(gamma_s_higher) < 1 and abs(gamma_s_lower) > 1:
        gamma_s = gamma_s_higher
    else:
        print("both potential magnitudes of gamma s are less than one. Try something else")

    gamma_l_lower = (b2+gam_l_sqrt)/(2*c2)
    gamma_l_higher =(b2-gam_l_sqrt)/(2*c2)

    if abs(gamma_l_lower) < 1 and abs(gamma_l_higher) > 1:
        gamma_l = gamma_l_lower
    elif abs(gamma_l_higher) < 1 and abs(gamma_l_lower) > 1:
        gamma_l = gamma_l_higher
    else:
        print("both potential magnitudes of gamma s are less than one. Try something else")

    gamma_lp = np.degrees(cm.phase(gamma_l))
    gamma_lm = abs(gamma_l)

    gamma_sp = np.degrees(cm.phase(gamma_s))
    gamma_sm = abs(gamma_s)


    print(f" Gamma S {gamma_sm} {gamma_sp} \n Gamma L {gamma_lm} {gamma_lp}")

    gamma_in = np.conjugate(gamma_s)
    gamma_out = np.conjugate(gamma_l)

    print(f'Gamma S {gamma_s}, gamma L {gamma_l}, gamma in {gamma_in}, gamma out {gamma_out}')
          
    gs = (1-m.pow(abs(gamma_s),2))/(m.pow(abs(1-(gamma_in*gamma_s)),2))
    g0 = m.pow(s21a,2)
    gl = (1-m.pow(abs(gamma_l),2))/(m.pow(abs(1-(s22*gamma_l)),2))
    print(f'GS {gs}, g0 {g0}, gl {gl}')
    gt = gs*g0*gl
    print(f'{gt} GT less')
    gt_full = 10*m.log10(gt)
    
    return gt_full, gamma_l, gamma_s



def get_gt_min_noise(s11,s12,s21,s22,gamma_opt): ## find the minimum possible noise this circuit can produce, will effect gain
    print(f'\nMIN NOISE\n')

    gamma_s = gamma_opt
    gamma_out = s22 + (s12*s21*gamma_s)/(1-(s11*gamma_s))
    gamma_l = np.conjugate(gamma_out)
    gamma_in = s11 + (s12*s21*gamma_l)/(1-(s22*gamma_l))

    print(f'{gamma_s} Gamma S, {gamma_l} gamma L, {gamma_in} gamma in, {gamma_out} gamma out')

    gs = (1-m.pow(abs(gamma_s),2))/(m.pow(abs(1-(gamma_in*gamma_s)),2))
    g0 = m.pow(s21a,2)
    gl = (1-m.pow(abs(gamma_l),2))/(m.pow(abs(1-(s22*gamma_l)),2))
    print(f'GS {gs}, g0 {g0}, gl {gl}')
    gt = gs*g0*gl
    print(f'{gt} GT less')
    gt_full = 10*m.log10(gt)
    return gt_full

def noise_calc(n, gamma_opt, r_n, f_min, z_0): ## calculates noise circle
    f          = f_min + ((4*r_n)/z_0)*n*(1/m.pow(abs(1+gamma_opt),2))
    f_1        = 10*m.log(f,10)
    center_f   = gamma_opt/(n+1)
    center_f1  = cm.polar(center_f)
    c_f_degree = np.degrees(center_f1[1])
    center_f   = abs(center_f)
    radius_f   = m.sqrt(n*(n+1-m.pow(abs(gamma_opt),2)))/(n+1)
    
    print(f"center f {center_f}, degrees {c_f_degree} \nf {f}, f dB {f_1}, radius {radius_f}")
    
    return center_f, radius_f, f_1

def gain_calc(s11, s11a, gs): ## Calculates gain circle
    center_s = (gs*np.conjugate(s11))/(1-(1-gs)*m.pow(s11a, 2))
    center_s1  = cm.polar(center_s)
    c_s_degree = np.degrees(center_s1[1])
    center_s = abs(center_s)
    print(f"center s {center_s}, degrees {c_s_degree}")
    radius_s = abs((m.sqrt(1-gs)*(1-m.pow(s11a,2)))/(1-(1-gs)*m.pow(s11a, 2)))
    return center_s, radius_s

def uni_case(s11,s12,s21,s22, gamma_opt, s11a,s12a,s21a,s22a, delta_d, d_abs, r_n): ## Unilateral case where S12 is close enough to zero, makes case so much simpler
    print(f'\n Unilateral Case\n')
    G_s = 1/m.pow((1-s11a),2)
    G_0 = m.pow(s21a,2)
    G_l = 1/m.pow((1-s22a),2)
    print(f'Gs {G_s}, G0 {G_0}, Gl {G_l}')
    
    gain_a_max = G_s * G_0 * G_l
    G_a_maxdB = 10*m.log10(gain_a_max)

    print(f'G A Max {gain_a_max} \nG A Max dB {G_a_maxdB}')
    gamma_in  = s11
    gamma_out = s22
    gamma_s   = np.conjugate(s11)
    gamma_l   = np.conjugate(s22)
    
    center_g, radius_g = gain_calc(s11, s11a, gs)
    center_f, radius_f, f_1 = noise_calc(n, gamma_opt, r_n, f_min, z_0)
    g_in_out  = cm.polar(gamma_in)
    g_out_out = cm.polar(gamma_out)
    g_s_out   = cm.polar(gamma_s)
    g_l_out   = cm.polar(gamma_l)
    print (f'Gamma in {g_in_out} Gamma out {g_out_out} \nGamma s {g_s_out} Gamma l {g_l_out}')

    return None

def conditionally_stable_vswr(vswr_in, vswr_out, gamma_s_guess): ## In case she asks for vswr for some reason
    print(f'\nVSWR Calculations\n')
    gamma_in_a   = abs((vswr_in-1)/(vswr_in+1))
    gamma_out_a   = abs((vswr_out-1)/(vswr_out+1))

    gamma_l      = np.conjugate(s22) + (np.conjugate(s12)*np.conjugate(s21)*np.conjugate(gamma_s_guess))/(1-(s11*np.conjugate(gamma_s_guess)))
    gamma_in     = s11 + (s12*s21*gamma_l)/(1-(s22*gamma_l))
    gamma_out    = s22 + (s12*s21*gamma_s_guess)/(1-(s11*gamma_s_guess)) ## Only going to guess one variable
    g_in_out  = cm.polar(gamma_in)
    g_out_out = cm.polar(gamma_out)
    g_s_out   = cm.polar(gamma_s_guess)
    g_l_out   = cm.polar(gamma_l)
    print (f'Gamma in {g_in_out} Gamma out {g_out_out} \nGamma s {g_s_out} Gamma l {g_l_out}')
    
    center_v_in    = (np.conjugate(gamma_in)*(1-m.pow(gamma_in_a,2)))/(1-m.pow(abs(gamma_in_a*gamma_in),2))
    center_v_inpo  = cm.polar(center_v_in)
    radius_v_in    = (abs(gamma_in_a)*(1-abs(m.pow(gamma_in,2))))/(1-m.pow(abs(gamma_in_a*gamma_in),2))
    print(f'Center In mag {center_v_inpo[0]} Center In degrees{center_v_inpo[1]} \nRadius in {radius_v_in}')

## This section is not working due to asking to change to a float as a complex.  
##    center_v_out   = (np.conjugate(gamma_out)*(1-m.pow(gamma_out_a,2)))/(1-m.pow(abs(gamma_out_a*gamma_out),2))
##    center_v_outpo = cm.polar(center_v_out)
##    radius_v_out   = (abs(gamma_out_a)*(1-abs(m.pow(gamma_out,2))))/(1-m.pow(abs(gamma_out_a*gamma_out),2))
##    print(f'Center Out mag {center_v_outpo[0]} Center Out degrees{center_v_outpo[1]} \nRadius out {radius_v_out}')
    return None

##def gp_gain_circle(s11,s12,s21,s22,s11a,s12a,s21a,s22a, delta_d, Gp, k): ##get the Gp circles for some reason, also doesn't work. 
##    
##    gp              = Gp/m.pow(s21a,2)
##    center_p        = (gp* np.conjugate(s22-(delta_d*np.conjugate(s11))))/(1+gp(m.pow(s22a,2)-m.pow(abs(delta_d),2)))
##    center_p_polar  = cm.polar(center_p)
##    radius_p        = m.sqrt(1-2*k(abs(s12*s21))*gp+(m.pow(abs(s12*s21),2)*m.pow(gp,2))) /abs((1+gp(m.pow(s22a,2)-m.pow(abs(delta_d),2))))
##    print(f'Center Gain P mag {center_p_polar[0]} Center Gain P degrees{center_p_polar[1]} \nRadius Gain P {radius_p}')
##    return None

###############
##MAIN Will get it organized someday
###############

##def main(s11,s12,s21,s22, f_min, gamma_opt, r_ohm):

## HARD CODED VALUES

s11 = cm.rect(0.6, np.radians(146))
s12 = cm.rect(0.017, np.radians(62))
s21 = cm.rect(1.97, np.radians(32))
s22 = cm.rect(0.52, np.radians(-63))
f_min = 2.5
gamma_opt = cm.rect(0.475, np.radians(166))

Rn = 9
z_0 = 50
k = 1
##GUESSING Values
gs = 0.8
n = .115
Gp = 0.8

s11a = abs(s11)
s12a = abs(s12)
s21a = abs(s21)
s22a = abs(s22)

vswr_in = 1.5
vswr_out = 1.5
print(f'{s11a} s11, {s12a}, {s21a}, {s22a}')

delta_d, d_abs = get_stability(s11,s12,s21,s22,s11a,s12a,s21a,s22a)
boolie = is_s12_unilateral(s11a, s12a, s21a, s22a)

print(f'Unilateral is {boolie}')

##if boolie:
##    uni_case(s11,s12,s21,s22, gamma_opt, s11a,s12a,s21a,s22a, delta_d, d_abs, Rn)
##else:
##    
##    gt, gamma_gt_s, gamma_gt_l = get_max_transducer_gain(s11,s12,s21,s22,s11a,s12a,s21a,s22a, delta_d, d_abs)
##    print(f'{gt} GT max transducer')
##
##    gt_min_noise = get_gt_min_noise(s11,s12,s21,s22, gamma_opt)
##    print(f'{gt_min_noise} GT min noise')
##
##    ##Needs N
##    center_f, radius_f, f = noise_calc(n, gamma_opt, Rn, f_min, z_0)
##
##    ##
##    center_g, radius_g = gain_calc(s11, s11a, gs)
##    print(f'{center_g} center, {radius_g} radius')
##    g_t_final = 10*m.log(gs,10) + gt
##    print(f"Gt in dB {g_t_final}")
##    gamma_s_final = cm.rect(0.04, np.radians(-127)) ## Must change as values change
##    
####    conditionally_stable_vswr(vswr_in, vswr_out, gamma_s_final)
####    gp_gain_circle(s11,s12,s21,s22,s11a,s12a,s21a,s22a, delta_d, Gp, k)
##    
##    gamma_s_finalp = np.degrees(cm.phase(gamma_s_final))
##    gamma_s_finalm = abs(gamma_s_final)
##    print(f'{gamma_s_finalm} gamma s m, {gamma_s_finalp} phase')
##    gamma_out_final= s22+(s12*s21*gamma_s_final)/(1-s11*gamma_s_final)
##    gamma_out_finalp = np.degrees(cm.phase(gamma_out_final))
##    gamma_out_finalm = abs(gamma_out_final)
##    print(f'{gamma_out_finalm} gamma out m, {gamma_out_finalp} phase')
##
##    gamma_l_final = np.conjugate(gamma_out_final)
##    gamma_l_finalp = np.degrees(cm.phase(gamma_l_final))
##    gamma_l_finalm = abs(gamma_l_final)
##    print(f'{gamma_l_finalm} gamma l m, {gamma_l_finalp} phase')
##    gamma_in_final= s11+(s12*s21*gamma_l_final)/(1-s22*gamma_l_final)
##    gamma_in_finalp = np.degrees(cm.phase(gamma_in_final))
##    gamma_in_finalm = abs(gamma_in_final)
##    print(f'{gamma_in_finalm} gamma in m, {gamma_in_finalp} phase')

gt, gamma_gt_s, gamma_gt_l = get_max_transducer_gain(s11,s12,s21,s22,s11a,s12a,s21a,s22a, delta_d, d_abs)
print(f'{gt} GT max transducer')

gt_min_noise = get_gt_min_noise(s11,s12,s21,s22, gamma_opt)
print(f'{gt_min_noise} GT min noise')

##Needs N
center_f, radius_f, f = noise_calc(n, gamma_opt, Rn, f_min, z_0)

##
center_g, radius_g = gain_calc(s11, s11a, gs)
print(f'{center_g} center, {radius_g} radius')
g_t_final = 10*m.log(gs,10) + gt
print(f"Gt in dB {g_t_final}")
gamma_s_final = cm.rect(0.04, np.radians(-127)) ## Must change as values change

####conditionally_stable_vswr(vswr_in, vswr_out, gamma_s_final)
####gp_gain_circle(s11,s12,s21,s22,s11a,s12a,s21a,s22a, delta_d, Gp, k)

gamma_s_finalp = np.degrees(cm.phase(gamma_s_final))
gamma_s_finalm = abs(gamma_s_final)
print(f'{gamma_s_finalm} gamma s m, {gamma_s_finalp} phase')
gamma_out_final= s22+(s12*s21*gamma_s_final)/(1-s11*gamma_s_final)
gamma_out_finalp = np.degrees(cm.phase(gamma_out_final))
gamma_out_finalm = abs(gamma_out_final)
print(f'{gamma_out_finalm} gamma out m, {gamma_out_finalp} phase')

gamma_l_final = np.conjugate(gamma_out_final)
gamma_l_finalp = np.degrees(cm.phase(gamma_l_final))
gamma_l_finalm = abs(gamma_l_final)
print(f'{gamma_l_finalm} gamma l m, {gamma_l_finalp} phase')
gamma_in_final= s11+(s12*s21*gamma_l_final)/(1-s22*gamma_l_final)
gamma_in_finalp = np.degrees(cm.phase(gamma_in_final))
gamma_in_finalm = abs(gamma_in_final)
print(f'{gamma_in_finalm} gamma in m, {gamma_in_finalp} phase')

##gamma_l_final = np.conjugate(gamma_out_final)
##gamma_l_finalp = np.degrees(cm.phase(gamma_l_final))
##gamma_l_finalm = abs(gamma_l_final)
##print(f'{gamma_l_finalm} gamma l m, {gamma_l_finalp} phase')
"""
After finding the gammas
1. Use smith chart to match the gamma_s and gamma_l to 0 using
    transmission lines and lumped elements ( The transimmsion lines act as
    the resistive part and the lumped elements to negate the reactive componet)
    i.e. The final impedance is r+jx so the lumped element is -jx    
2. Change the lumped elements to stubs (Series or shunt, maybe focus on shunt)
    have t
"""






