#!/usr/bin/env python3
""" EvtGen .dec file parser

This module parses EvtGen decay model parameters from a .dec file and
converts them to helicity, transversity and cartesian conventions and
prints these.

"""

import sys
import os
import re
import math

__author__ = "Daniel Cervenkov"
__email__ = "cervenkov@ipnp.mff.cuni.cz"
__date__ = "2012/04/24"


def print_help():
    print("USAGE: " + sys.argv[0] + " FILE")


# Reads in all parameters of the relevant model
def get_raw_pars(dec_file):
    f = open(dec_file, 'r')
    text = f.read()
    f.close()

    text = re.sub('#.*\n', '', text)

    model_str = 'SVV_NONCPEIGEN'
    model_delimiter = ';'

    beg = text.find(model_str)
    end = text.find(model_delimiter, beg+1)

    if beg > -1 and end > beg:
        model = text[beg:end]
        pars = model.split()
        del pars[0]

        # If a parameter is a variable name, replace it by its value
        for index, par in enumerate(pars):
            if not par.replace('.', '0').replace('e', '0').replace('-', '0').isdigit():
                beg_const = text.find(par)
                end_const = text.find('\n', beg+1)
                if beg_const < beg and end_const > beg_const:
                    constant = text[beg_const:end_const].split()
                    pars[index] = constant[1]

        for index, par in enumerate(pars):
            if not par.replace('.', '0').replace('e', '0').replace('-', '0').isdigit():
                print('ERROR: Not all parameters are defined in the decay file.')
                return 0

        return pars

    else:
        print("'{0} ... {1}' not found".format(model_str, model_delimiter))
        return 0


# Calculates helicity strong phases and r_lambda from unphysical 'amplitudes'
def calcHelPars(pars):
    new_pars = pars[3:9]
    for i in range(3):
        new_pars.append(round(float(pars[9+2*i])/float(pars[3+2*i]), 3))
    for i in range(3):
        new_pars.append(round(float(pars[10+2*i])-float(pars[4+2*i]), 3))

    return new_pars


# Calculates transversity amplitudes, strong phases and r_m from unphysical 'amplitudes'
def calcTransPars(pars):
    new_pars = pars[3:9]

    hp = complex(float(pars[3])*math.cos(float(pars[4])),
                 float(pars[3])*math.sin(float(pars[4])))

    h0 = complex(float(pars[5])*math.cos(float(pars[6])),
                 float(pars[5])*math.sin(float(pars[6])))

    hm = complex(float(pars[7])*math.cos(float(pars[8])),
                 float(pars[7])*math.sin(float(pars[8])))

    hps = complex(float(pars[9])*math.cos(float(pars[10])),
                  float(pars[9])*math.sin(float(pars[10])))

    h0s = complex(float(pars[11])*math.cos(float(pars[12])),
                  float(pars[11])*math.sin(float(pars[12])))

    hms = complex(float(pars[13])*math.cos(float(pars[14])),
                  float(pars[13])*math.sin(float(pars[14])))

    ap = (hp + hm)/math.sqrt(2)
    a0 = h0
    at = (hp - hm)/math.sqrt(2)

    rhop = (hps + hms)/(hp + hm)
    rho0 = h0s/h0
    rhot = (hps - hms)/(hp - hm)

    new_pars[0] = round(abs(ap), 3)
    new_pars[1] = round(math.atan2(ap.imag, ap.real), 2)
    new_pars[2] = round(abs(a0), 3)
    new_pars[3] = round(math.atan2(a0.imag, a0.real), 2)
    new_pars[4] = round(abs(at), 3)
    new_pars[5] = round(math.atan2(at.imag, at.real), 2)

    new_pars.append(round(abs(rhop), 3))
    new_pars.append(round(abs(rho0), 3))
    new_pars.append(round(abs(rhot), 3))

    new_pars.append(round(math.atan2(rhop.imag, rhop.real), 3))
    new_pars.append(round(math.atan2(rho0.imag, rho0.real), 3))
    new_pars.append(round(math.atan2(rhot.imag, rhot.real), 3))

    return new_pars


# Calculates transversity amplitudes, strong phases and r_m from unphysical 'amplitudes'
def calcTransCartPars(pars):
    new_pars = pars[3:9]

    hp = complex(float(pars[3])*math.cos(float(pars[4])),
                 float(pars[3])*math.sin(float(pars[4])))

    h0 = complex(float(pars[5])*math.cos(float(pars[6])),
                 float(pars[5])*math.sin(float(pars[6])))

    hm = complex(float(pars[7])*math.cos(float(pars[8])),
                 float(pars[7])*math.sin(float(pars[8])))

    hps = complex(float(pars[9])*math.cos(float(pars[10])),
                  float(pars[9])*math.sin(float(pars[10])))

    h0s = complex(float(pars[11])*math.cos(float(pars[12])),
                  float(pars[11])*math.sin(float(pars[12])))

    hms = complex(float(pars[13])*math.cos(float(pars[14])),
                  float(pars[13])*math.sin(float(pars[14])))

    ap = (hp + hm)/math.sqrt(2)
    a0 = h0
    at = (hp - hm)/math.sqrt(2)

    phiw = 2*float(pars[1])+float(pars[2]) - 3.14159

    weak_part = complex(math.cos(phiw), -math.sin(phiw))
    weak_part_bar = complex(math.cos(phiw), math.sin(phiw))

    rhop = (hps + hms)/(hp + hm) * weak_part
    rho0 = h0s/h0 * weak_part
    rhot = (hps - hms)/(hp - hm) * weak_part

    new_pars[0] = round(abs(ap), 3)
    new_pars[1] = round(math.atan2(ap.imag, ap.real), 2)
    new_pars[2] = round(abs(a0), 3)
    new_pars[3] = round(math.atan2(a0.imag, a0.real), 2)
    new_pars[4] = round(abs(at), 3)
    new_pars[5] = round(math.atan2(at.imag, at.real), 2)

    rhop_bar = rhop * weak_part_bar * weak_part_bar
    rho0_bar = rho0 * weak_part_bar * weak_part_bar
    rhot_bar = rhot * weak_part_bar * weak_part_bar

    new_pars.append(round(rhop.real, 4))
    new_pars.append(round(rho0.real, 4))
    new_pars.append(round(rhot.real, 4))

    new_pars.append(round(rhop.imag, 4))
    new_pars.append(round(rho0.imag, 4))
    new_pars.append(round(rhot.imag, 4))

    new_pars.append(round(rhop_bar.real, 4))
    new_pars.append(round(rho0_bar.real, 4))
    new_pars.append(round(rhot_bar.real, 4))

    new_pars.append(round(rhop_bar.imag, 4))
    new_pars.append(round(rho0_bar.imag, 4))
    new_pars.append(round(rhot_bar.imag, 4))

    return new_pars


def print_pars(pars):
    dm = pars[0]
    phiw = 2*float(pars[1])+float(pars[2]) - 3.14159

    print('{0}: {1}'.format("dm".rjust(4), dm))
    print('{0}: {1}'.format("phiw".rjust(4), phiw))
    print('')

    hel_pars = calcHelPars(pars)
    trans_pars = calcTransPars(pars)
    trans_cart_pars = calcTransCartPars(pars)

    if hel_pars:
        labels = ['hp', 'hpa', 'h0', 'h0a', 'hm', 'hma',
                  'rp', 'r0',  'rm', 'sp',  's0', 'sm']
        for i, label in enumerate(labels):
            print('{0}: {1}'.format(label.rjust(4), hel_pars[i]))
        print('')

    if trans_pars:
        labels = ['ap', 'apa', 'a0', 'a0a', 'at', 'ata',
                  'rp', 'r0',  'rt', 'sp',  's0', 'st']
        for i, label in enumerate(labels):
            print('{0}: {1}'.format(label.rjust(4), trans_pars[i]))
        print('')

    if trans_cart_pars:
        labels = ['ap',  'apa', 'a0',  'a0a', 'at',  'ata',
                  'xp',  'x0',  'xt',  'yp',  'y0',  'yt',
                  'xpb', 'x0b', 'xtb', 'ypb', 'y0b', 'ytb']
        for i, label in enumerate(labels):
            print('{0}: {1}'.format(label.rjust(4), trans_cart_pars[i]))


def main():
    if len(sys.argv) != 2:
        print("ERROR: You must supply exactly 1 argument")
        print_help()
        sys.exit(1)

    dec_file = sys.argv[1]
    if os.path.exists(dec_file):
        pars = get_raw_pars(dec_file)
        if pars:
            print_pars(pars)
    else:
        print("ERROR:", dec_file, "doesn't exist.")


if __name__ == "__main__":
    main()
