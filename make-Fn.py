#! /usr/bin/env python 

# Python 2.7

import sys


ABG = {
        'S0_odd' : {'alpha':132.5,'beta':163.4,'gamma':118.2},      # ground state geometry; odd/even are needed to create an alternant conformer
        'S0_even' : {'alpha':132.5,'beta':-163.4,'gamma':118.2},    # ground state geometry
        'S1_perfect' : {'alpha':120.0,'beta':180.0,'gamma':120.0},  # make the planes containing fluorenes parallel
        'S1_normal' : {'alpha':128.8,'beta':180.0,'gamma':120.0},   # arrange fluorenes as in F2S1
        'S1_neutral' : {'alpha':132.5,'beta':180.0,'gamma':120.0}   # use alpha as in F2S0
        }



#### STRUCTURAL CONSTANTS #####

S_FIRST_ME_BRIDGE = """C
C     1    1.546
Bq    2    0.968   1      120.0
H     1    1.09    2      110.7   3        0.0
H     1    1.09    2      110.7   3        120.0
H     1    1.09    2      110.7   3        -120.0"""

S_FLUORENE = """C     %i   1.178   %i     90.0    %i       -90.0
C     %i   1.402   %i     71.4    %i       180.0
C     %i   1.390   %i     120.0   %i       180.0
C     %i   1.390   %i     120.0   %i       0.0
C     %i   1.390   %i     120.0   %i       0.0
C     %i   1.390   %i     120.0   %i       0.0
C     %i   1.178   %i     90.0    %i       90.0
C     %i   1.390   %i     71.4    %i       180.0
C     %i   1.390   %i     120.0   %i       180.0
C     %i   1.390   %i     120.0   %i       0.0
C     %i   1.390   %i     120.0   %i       0.0
C     %i   1.390   %i     120.0   %i       0.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0
H     %i   1.08    %i     120.0   %i       180.0"""

S_BRIDGE = """C    %i   1.55    %i   _(alpha)f   %i   180.0
Bq   %i   0.664   %i   _(gamma)f   %i   _(beta)f
H    %i   0.867   %i   90.0        %i   -90.0
H    %i   0.867   %i   90.0        %i   90.0
C    %i   1.55    %i   _(gamma)f   %i   180.00
Bq   %i   0.968   %i   _(alpha)f   %i   _(beta)f"""

S_LAST_ME = """C   %i   1.546   %i  120.0   %i     180.0
H   %i   1.09    %i  110.7   %i     0.0
H   %i   1.09    %i  110.7   %i     120.0
H   %i   1.09    %i  110.7   %i     240.0"""


######### FUNC ############

def bridge(start_i,i_Bq,i_top,i_bridge):
    """
     CH2
    / \
   F   C
       |
       Bq

    New bridge unit will be constructed in the framework of CH2 (i_bridge), C (i_top), and Bq (i_bq)
    """

    ref_ati = (
          i_top,           i_Bq,            i_bridge,
          start_i,         i_top,           i_Bq,
          start_i+1,       start_i,         i_top,
          start_i+1,       start_i,         i_top,
          start_i,         start_i+1,       i_top,
          start_i+4,       start_i,         start_i+1
    )
    D = {'s' : S_BRIDGE % ref_ati,
         'i_bridge' : start_i,
         'i_top' : start_i+4,
         'i_Bq' : start_i+5
            }
    return D


def bridge_abg(s,abg):
    """
    Inserts values of alpha, beta, and gamma
    """
    s = s.replace('_','%')
    return s % abg


def fluorene(start_i,i_Bq,i_top,i_bridge):
    """
     CH2
    / \
   F   C
       |
       Bq

    Next fluorene units will be constructed in the framework of CH2 (i_bridge), C (i_top), and Bq (i_bq)
    """

    start_i -= 4 # calculate offset for the existing ref_ati
    i_Bq -= start_i
    i_top -= start_i
    i_bridge -= start_i

    ref_ati = (
          i_Bq,    i_top,    i_bridge,
          4,       i_Bq,     i_top,
          5,       4,        i_Bq,
          6,       5,        4,
          7,       6,        5,
          8,       7,        6,
          i_Bq,    i_top,    i_bridge,
          10,      i_Bq,     i_top,
          11,      10,       i_top,
          12,      11,       10,
          13,      12,       11,
          14,      13,       12,
          9,       8,        4,
          8,       9,        7,
          7,       8,        6,
          6,       7,        5,
          12,      13,       11,
          13,      14,       12,
          14,      15,       13,
          15,      10,       14
    )

    ati = [x + start_i for x in ref_ati]
    return S_FLUORENE % tuple(ati)




def closing_me(start_i,i_Bq,i_top,i_bridge):
    """
     CH2
    / \
       C
       |
       Bq

    Makes a close Me
    """

    ref_ati = (
        i_top,     i_Bq , i_bridge,
        start_i, i_top, i_Bq,
        start_i, i_top, i_Bq,
        start_i, i_top, i_Bq
    )

    return S_LAST_ME % ref_ati


def print_z(s):
    print s
    return len(s.split('\n'))



######################### MAIN ###################
##################################################

"""
FN = 6
Fstacked = 6
el_state = 'S1_normal'
"""

# COMMAND LINE ARGUMENTS
FN = int(sys.argv[1]) # FN, number of fluorenes
Fstacked = int(sys.argv[2]) # how many units of FN should be fully stacked?
el_state = sys.argv[3] # Determines the angle between adjacent stacked fluorene units; valid choices are the various 'S1_*' keys in ABG


start_i = 1

# Print first Me and bridge
start_i += print_z(S_FIRST_ME_BRIDGE)

# Print first fluorene
s = fluorene(start_i=start_i,i_Bq=3,i_top=2,i_bridge=1)
start_i += print_z(s)

state_suffix = ''
i_Bq = 3
i_top = 2
i_bridge = 1

for i in range(1,FN):
    if i >= Fstacked:
        el_state = 'S0'
        if ((i - Fstacked) % 2 == 0):
            state_suffix = '_even'
        else:
            state_suffix = '_odd'

    # Print next bridge
    bridge_params = bridge(start_i=start_i,i_Bq=i_Bq,i_top=i_top,i_bridge=i_bridge)
    s = bridge_abg(bridge_params['s'],ABG[el_state+state_suffix])
    start_i += print_z(s)


    # Print next fluorene
    i_Bq=bridge_params['i_Bq']
    i_top=bridge_params['i_top']
    i_bridge=bridge_params['i_bridge']
    s = fluorene(start_i=start_i,i_Bq=i_Bq,i_top=i_top,i_bridge=i_bridge)
    start_i += print_z(s)

# Print closing Me group
s = closing_me(start_i=start_i,i_Bq=bridge_params['i_Bq'],i_top=bridge_params['i_top'],i_bridge=bridge_params['i_bridge'])
start_i += print_z(s)

