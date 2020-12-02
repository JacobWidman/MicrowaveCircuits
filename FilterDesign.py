"""

For making ideal filters with microstrip transistions

ideal filters invoving prototyping a low pass filter and converting to bandpass

To start problem:
    - Look at the charts for the right one
    (omega/ omega_cutoff ) -1  find N and plug and chug here

(Voltage) N = 10
(Current) N = 5
"""

import math
from scipy import constants as spc
## HARDCODED


def voltage_source(var_list, initial_impedance, res_freq, freq, wave_len):
    print(f'VOLTAGE')
##capacitors
    cap1 = var_list[1] / (initial_impedance * res_freq) ##C1`

    cap3 = var_list[3] / (initial_impedance * res_freq) ##C3`

    cap5 = var_list[5] / (initial_impedance * res_freq) ##C5

    cap7 = var_list[7] / (initial_impedance * res_freq) ##C7

    cap9 = var_list[9] / (initial_impedance * res_freq) ##C9

    print(f'Cap1 {cap1} Cap3 {cap3} \nCap5 {cap5} Cap7 {cap7} \nCap9 {cap9}  ')

    ##inductors
    ind2 = var_list[2] * initial_impedance / res_freq ##I2
    ind4 = var_list[4] * initial_impedance / res_freq ##I4
    ind6 = var_list[6] * initial_impedance  / res_freq##I6
    ind8 = var_list[8] * initial_impedance / res_freq ##I8
    ind10 = var_list[10] * initial_impedance / res_freq ##I10


    print(f'Ind2 {ind2} Ind4 {ind4} \nInd6 {ind6} Ind8 {ind8}\nInd10 {ind10}  ')

    ##microstrip code for lowpass
    ##wider line higher cap
    ## thinner line higher ind

    ##capacitive section

    ##c1 = (1/Z_width)* (length_line/(freq*wave_len))

    c1_ratio = Z_width * freq * cap1 ##need this to double check

    c1_len = c1_ratio * wave_len

    print(f'cap 1 ratio {c1_ratio} ')

    ##inducitve section

    i2_ratio = (ind2*freq)/Z_height
    i2_len = i2_ratio * wave_len
    print(f'ind 2 ratio {i2_ratio} ')

    return None

def current_source(var_list, initial_impedance, res_freq, freq, wave_len):
    print(f'CURRENT')
####resistors the same R L C L C L R for current loads
##
####Inductors
    ##Aint broke don't fix, but could use a for loop for assigning things
    ind1 = var_list[1] * initial_impedance / res_freq ##I1 == I5` Not always

    ind3 = var_list[3] * initial_impedance / res_freq ##I3

####capaciters
    cap2 = var_list[2] / (initial_impedance * res_freq) ##I2 == I4 not always

    ##microstrip code for lowpass
    ##capacitive section
    ##c1 = (1/Z_width)* (length_line/(freq*wave_len))

    c2_ratio = Z_width * freq * cap2 ##need this to double check

    c2_len = c2_ratio * wave_len

    print(f'cap 2 ratio {c2_ratio} ')

    ##inducitve section

    i1_ratio = (ind1*freq)/Z_height
    i1_len = i1_ratio * wave_len
    print(f'ind 1 ratio {i1_ratio} ')

    return None

########################
## MAIN
########################




initial_impedance = 50 ##Z_0 Ohms
##bw_delta = 0.1 ##bandwidth delta
freq = 2.4 * math.pow(10, 9) ##Hz THIS IS CUTOFF FREQ
Z_height = 20 ##Ohms
Z_width = 120 ##Ohms


"""
It might help to know that Z_0 = math.sqrt(ind/cap)

Remember that this is for a voltage source,
for a current source the series inductor is the g1, g3 etc and cap is g2, g4

also g0 and the g[n+1] are resistors and should be 1

should be smaller than the 10th of the wavelength
"""

wave_len = spc.c/freq
res_freq =  freq * 2 * math.pi

##this var list is for a max flat filter of 10
##edit as necessary, as I am not going to put all the tables here
var_list = {
    0: 1,
    1: 0.6305,
    2: 0.3002,
    3: 0.2384,
    4: 0.2066,
    5: 0.1808,
    6: 0.1539,
    7: 0.1240,
    8: 0.0911,
    9: 0.0557,
    10: 0.0187,
    11:1,
    }

resistor1 = var_list[0] * initial_impedance ##R
resistor2 = var_list[11] * initial_impedance
print(f'RS {resistor1} RL {resistor2}')

voltage_source(var_list, initial_impedance, res_freq, freq, wave_len)
current_source(var_list, initial_impedance, res_freq, freq, wave_len)




####bandpass conversion
####resistors remain the same 
##
####conversion for Caps
####C1
##bp_ind1 = bw_delta / (res_freq * cap1)
##bp_cap1 = cap1 / (res_freq * bw_delta)
##print("BPIN 1 " + str(bp_ind1) + " BPCAP1 " + str(bp_cap1))
##
####C3
##bp_ind3 = bw_delta / (res_freq * cap3)
##bp_cap3 = cap3 / (res_freq * bw_delta)
##print("BPIN 3 " + str(bp_ind3) + " BPCAP3 " + str(bp_cap3))
##
####C5
##bp_ind5 = bw_delta / (res_freq * cap5)
##bp_cap5 = cap5 / (res_freq * bw_delta)
##print("BPIN 5 " + str(bp_ind5) + " BPCAP5 " + str(bp_cap5))
##
####C7
##bp_ind7 = bw_delta / (res_freq * cap7)
##bp_cap7 = cap7 / (res_freq * bw_delta)
##print("BPIN 7 " + str(bp_ind7) + " BPCAP7 " + str(bp_cap7))
##
####C9
##bp_ind9 = bw_delta / (res_freq * cap9)
##bp_cap9 = cap9 / (res_freq * bw_delta)
##print("BPIN 9 " + str(bp_ind9) + " BPCAP9 " + str(bp_cap9))
##
##
####Inductors
####I2
##bp_ind2 = ind2 / (res_freq * bw_delta)
##bp_cap2 = bw_delta / (res_freq * ind2)
##print("BPIN 2 " + str(bp_ind2) + " BPCAP2 " + str(bp_cap2))
####I4
##bp_ind4 = ind4 / (res_freq * bw_delta)
##bp_cap4 = bw_delta / (res_freq * ind4)
##print("BPIN4  " + str(bp_ind4) + " BPCAP4 " + str(bp_cap4))
####I6
##bp_ind6 = ind6 / (res_freq * bw_delta)
##bp_cap6 = bw_delta / (res_freq * ind6)
##print("BPIN 6 " + str(bp_ind6) + " BPCAP6 " + str(bp_cap6))
####I8
##bp_ind8 = ind8 / (res_freq * bw_delta)
##bp_cap8 = bw_delta / (res_freq * ind8)
##print("BPIN 8 " + str(bp_ind8) + " BPCAP8 " + str(bp_cap8))
####I10
##bp_ind10 = ind10 / (res_freq * bw_delta)
##bp_cap10 = bw_delta / (res_freq * ind10)
##print("BPIN 10 " + str(bp_ind10) + " BPCAP10 " + str(bp_cap10))
##



####microstrip line
##
####for inductors i2 = i4
##char_imp2 = (4 * res_freq * bp_ind2)/ math.pi
##dist2 = (3*math.pow(10,8))/ (4*freq)
##
##print("Z02, L " + str(char_imp2) + " " + str(dist2))
##
####for capacitors c1 = c5
##char_imp1 = res_freq * bp_ind1
##dist1 = ((3*math.pow(10,8))/ (4*freq)) + ((bp_cap1*char_imp1)/(2*math.pi)) 
##print("Z01, L " + str(char_imp1) + " " + str(dist1))
####c3
##char_imp3 = res_freq * bp_ind3
##dist3 = ((3*math.pow(10,8))/ (4*freq)) + ((bp_cap3*char_imp3)/(2*math.pi)) 
##print("Z03, L " + str(char_imp3) + " " + str(dist3))
##

"""
#########################################################################################

"""


##

##
####bandpass conversion
####resistors remain the same 
##
####conversion for C2 == C4
##bp_ind2 = bw_delta / (res_freq * cap2)
##bp_cap2 = cap2 / (res_freq * bw_delta)
##
####Ind1 == I5
##bp_ind1 = (ind1*bw_delta )/ (res_freq )
##bp_cap1 = bw_delta / (res_freq * ind1)
##
##
####I3
##bp_ind3 = (ind3*bw_delta )/ (res_freq )
##bp_cap3 = bw_delta / (res_freq * ind3)










