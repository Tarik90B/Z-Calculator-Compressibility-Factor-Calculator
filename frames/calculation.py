"""
Z – Compressibility Factor Calculator
Copyright (c) 2026 Tarik

Independent implementation of AGA8 according to ISO 12213.
Licensed under the MIT License.
"""

import tkinter as tk
import numpy as np
import datetime
import csv
from collections import defaultdict
from utils import resource_path
from pathlib import Path

archive_file = resource_path("archive/z_factor_archive.csv", writable=True)
b1_file = resource_path("data/Equation_of_state_parameters.csv", writable=False)
b2_file = resource_path("data/Characterization_parameters.csv", writable=False)
b3_file = resource_path("data/Binary_interaction_parameter_values.csv", writable=False)

def save_to_archive(XI, T, p, Z):
    header = [
        "datetime", "T_K", "p_MPa",
        "CH4", "N2", "CO2", "C2H6", "C3H8", "H20", "H2S", "H2",
        "CO", "O2", "iC4H10", "nC4H10", "iC5H12", "nC5H12",
        "nC6H14", "nC7H16", "nC8H18", "nC9H20", "nC10H22",
        "He", "Ar", "Z"
    ]

    file_path = Path(archive_file)
    write_header = not file_path.exists()

    with file_path.open("a", newline="") as file:
        writer = csv.writer(file)
        if write_header:
            writer.writerow(header)
        writer.writerow([
            datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            T, p, *XI[:21], Z
        ])

class CalculationFrame(tk.Frame):
    def __init__(self, parent, controller):
        super().__init__(parent)

        headline_label = tk.Label(
            self, 
            text='Input following to calculate Z-factor:\n - The absolute temperature T [K] \n - The absolute pressure p [MPa]\n - The mole fraction of each component Xi of the gas mixture', 
            font=('bold', 14), 
            pady=10, 
            justify='left', 
            bg='lightyellow',
            bd=2,
            relief='solid'
            )
        headline_label.grid(row=1, column=0, columnspan=5, sticky='we', padx=10, pady=10)

        header_label = tk.Label(self, text='Mole fractions of mixture', font=('bold', 13))
        header_label.grid(row=2, column=0, columnspan=2, sticky='w', padx=10, pady=5)

        header_label = tk.Label(self, text='Temperature and Pressure', font=('bold', 13))
        header_label.grid(row=2, column=2, columnspan=2, sticky='w', padx=10, pady=5)

        self.grid_columnconfigure(0, weight=0)
        self.grid_columnconfigure(1, weight=0)
        self.grid_columnconfigure(2, weight=0)
        self.grid_columnconfigure(3, weight=0)
        self.grid_columnconfigure(4, weight=1)

        # Equation_of_state_parameters

        def load_coefficients_B1():

            an, bn, cn, kn, un, gn, qn, fn, sn, wn = ([] for _ in range(10))

            with b1_file.open(newline="") as file:
                reader = csv.DictReader(file)
                for row in reader: 
                    an.append(float(row["an"]))
                    bn.append(float(row["bn"]))
                    cn.append(float(row["cn"]))
                    kn.append(float(row["kn"]))
                    un.append(float(row["un"]))
                    gn.append(float(row["gn"]))
                    qn.append(float(row["qn"]))
                    fn.append(float(row["fn"]))
                    sn.append(float(row["sn"]))
                    wn.append(float(row["wn"]))
                
                return (np.array(an), np.array(bn), np.array(cn), np.array(kn),
                np.array(un), np.array(gn), np.array(qn), np.array(fn),
                np.array(sn), np.array(wn))

        self.an, self.bn, self.cn, self.kn, self.un, self.gn, self.qn, self.fn, self.sn, self.wn = load_coefficients_B1()

        # Characterization parameters

        def load_coefficients_B2():

            Xi, Mi, Ei, Ki, Gi, Qi, Fi, Si, Wi = ([] for _ in range(9))

            with b2_file.open(newline="") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    Xi.append(float(row["Xi"]))
                    Mi.append(float(row["Mi"]))
                    Ei.append(float(row["Ei"]))
                    Ki.append(float(row["Ki"]))
                    Gi.append(float(row["Gi"]))
                    Qi.append(float(row["Qi"]))
                    Fi.append(float(row["Fi"]))
                    Si.append(float(row["Si"]))
                    Wi.append(float(row["Wi"]))
                return (np.array(Xi), np.array(Mi), np.array(Ei), np.array(Ki), 
                        np.array(Gi), np.array(Qi), np.array(Fi),
                        np.array(Si), np.array(Wi))

        self.Xi, self.Mi, self.Ei, self.Ki, self.Gi, self.Qi, self.Fi, self.Si, self.Wi = load_coefficients_B2()

        # Binary interaction parameter values

        def load_coefficients_B3():

            Eijast = defaultdict(dict)
            Uij = defaultdict(dict)
            Kij = defaultdict(dict)
            Gijast = defaultdict(dict)

            with b3_file.open(newline="") as file:
                reader = csv.DictReader(file)
                for row in reader:
                    i = int(row["i"])
                    j = int(row["j"])

                    Eijast[i][j] = float(row["Eij"])
                    Uij[i][j] = float(row["Uij"])
                    Kij[i][j] = float(row["Kij"])
                    Gijast[i][j] = float(row["Gij"])
                    

                return (i, j, Eijast, Uij, Kij, Gijast)

        self.i, self.j, self.Eijast, self.Uij, self.Kij, self.Gijast = load_coefficients_B3()

        # Input required values from user

        # Temperature input
        temperature_value = tk.DoubleVar()
        temperature_label = tk.Label(self, text='Temperature [K]:', font=('bold', 12))
        temperature_label.grid(row=3, column=2, sticky='w', padx=10, pady=5)
        temperature_entry = tk.Entry(self, textvariable=temperature_value)
        temperature_entry.grid(row=3, column=3)

        # Pressure input
        pressure_value = tk.DoubleVar()
        pressure_label = tk.Label(self, text='Pressure [MPa]:', font=('bold', 12))
        pressure_label.grid(row=4, column=2, sticky='w', padx=10, pady=5)
        pressure_entry = tk.Entry(self, textvariable=pressure_value)
        pressure_entry.grid(row=4, column=3)

        # Methane mole fraction input
        methane_fraction_value = tk.DoubleVar()
        methane_fraction_label = tk.Label(self, text='Methane [CH4]:', font=('bold', 12))
        methane_fraction_label.grid(row=3, column=0, sticky='w', padx=10, pady=5)
        methane_fraction_entry = tk.Entry(self, textvariable=methane_fraction_value)
        methane_fraction_entry.grid(row=3, column=1, sticky='w')

        # Ethane mole fraction input
        ethane_fraction_value = tk.DoubleVar()
        ethane_fraction_label = tk.Label(self, text='Ethane [C2H6]:', font=('bold', 12))
        ethane_fraction_label.grid(row=4, column=0, sticky='w', padx=10, pady=5)
        ethane_fraction_entry = tk.Entry(self, textvariable=ethane_fraction_value)
        ethane_fraction_entry.grid(row=4, column=1, sticky='w')

        # Propane mole fraction input
        propane_fraction_value = tk.DoubleVar()
        propane_fraction_label = tk.Label(self, text='Propane [C3H8]:', font=('bold', 12))
        propane_fraction_label.grid(row=5, column=0, sticky='w', padx=10, pady=5)
        propane_fraction_entry = tk.Entry(self, textvariable=propane_fraction_value)
        propane_fraction_entry.grid(row=5, column=1, sticky='w')

        # iso-Butane mole fraction input
        isoButane_fraction_value = tk.DoubleVar()
        isoButane_fraction_label = tk.Label(self, text='iso-Butane [C4H10]:', font=('bold', 12))
        isoButane_fraction_label.grid(row=6, column=0, sticky='w', padx=10, pady=5)
        isoButane_fraction_entry = tk.Entry(self, textvariable=isoButane_fraction_value)
        isoButane_fraction_entry.grid(row=6, column=1, sticky='w')

        # n-Butane mole fraction input
        nButane_fraction_value = tk.DoubleVar()
        nButane_fraction_label = tk.Label(self, text='n-Butane [C4H10]:', font=('bold', 12))
        nButane_fraction_label.grid(row=7, column=0, sticky='w', padx=10, pady=5)
        nButane_fraction_entry = tk.Entry(self, textvariable=nButane_fraction_value)
        nButane_fraction_entry.grid(row=7, column=1, sticky='w')

        # iso-Pentane mole fraction input
        isoPentane_fraction_value = tk.DoubleVar()
        isoPentane_fraction_label = tk.Label(self, text='iso-Pentane [C5H12]:', font=('bold', 12))
        isoPentane_fraction_label.grid(row=8, column=0, sticky='w', padx=10, pady=5)
        isoPentane_fraction_entry = tk.Entry(self, textvariable=isoPentane_fraction_value)
        isoPentane_fraction_entry.grid(row=8, column=1, sticky='w')

        # n-Pentane mole fraction input
        nPentane_fraction_value = tk.DoubleVar()
        nPentane_fraction_label = tk.Label(self, text='n-Pentane [C5H12]:', font=('bold', 12))
        nPentane_fraction_label.grid(row=9, column=0, sticky='w', padx=10, pady=5)
        nPentane_fraction_entry = tk.Entry(self, textvariable=nPentane_fraction_value)
        nPentane_fraction_entry.grid(row=9, column=1, sticky='w')

        # n-Hexane mole fraction input
        nHexane_fraction_value = tk.DoubleVar()
        nHexane_fraction_label = tk.Label(self, text='n-Hexane [C6H14]:', font=('bold', 12))
        nHexane_fraction_label.grid(row=10, column=0, sticky='w', padx=10, pady=5)
        nHexane_fraction_entry = tk.Entry(self, textvariable=nHexane_fraction_value)
        nHexane_fraction_entry.grid(row=10, column=1, sticky='w')

        # n-Heptane mole fraction input
        nHeptane_fraction_value = tk.DoubleVar()
        nHeptane_fraction_label = tk.Label(self, text='n-Heptane [C7H16]:', font=('bold', 12))
        nHeptane_fraction_label.grid(row=11, column=0, sticky='w', padx=10, pady=5)
        nHeptane_fraction_entry = tk.Entry(self, textvariable=nHeptane_fraction_value)
        nHeptane_fraction_entry.grid(row=11, column=1, sticky='w')

        # n-Octane mole fraction input
        nOctane_fraction_value = tk.DoubleVar()
        nOctane_fraction_label = tk.Label(self, text='n-Octane [C8H18]:', font=('bold', 12))
        nOctane_fraction_label.grid(row=12, column=0, sticky='w', padx=10, pady=5)
        nOctane_fraction_entry = tk.Entry(self, textvariable=nOctane_fraction_value)
        nOctane_fraction_entry.grid(row=12, column=1, sticky='w')

        # n-Nonane mole fraction input
        nNonane_fraction_value = tk.DoubleVar()
        nNonane_fraction_label = tk.Label(self, text='n-Nonane [C9H20]:', font=('bold', 12))
        nNonane_fraction_label.grid(row=13, column=0, sticky='w', padx=10, pady=5)
        nNonane_fraction_entry = tk.Entry(self, textvariable=nNonane_fraction_value)
        nNonane_fraction_entry.grid(row=13, column=1, sticky='w')

        # n-Decane mole fraction input
        nDecane_fraction_value = tk.DoubleVar()
        nDecane_fraction_label = tk.Label(self, text='n-Decane [C10H22]:', font=('bold', 12))
        nDecane_fraction_label.grid(row=14, column=0, sticky='w', padx=10, pady=5)
        nDecane_fraction_entry = tk.Entry(self, textvariable=nDecane_fraction_value)
        nDecane_fraction_entry.grid(row=14, column=1, sticky='w')

        # Carbon_dioxide mole fraction input
        carbon_dioxide_fraction_value = tk.DoubleVar()
        carbon_dioxide_fraction_label = tk.Label(self, text='Carbon dioxide [CO2]:', font=('bold', 12))
        carbon_dioxide_fraction_label.grid(row=15, column=0, sticky='w', padx=10, pady=5)
        carbon_dioxide_fraction_entry = tk.Entry(self, textvariable=carbon_dioxide_fraction_value)
        carbon_dioxide_fraction_entry.grid(row=15, column=1, sticky='w')

        # Nitrogen mole fraction input
        nitrogen_fraction_value = tk.DoubleVar()
        nitrogen_fraction_label = tk.Label(self, text='Nitrogen [N2]:', font=('bold', 12))
        nitrogen_fraction_label.grid(row=16, column=0, sticky='w', padx=10, pady=5)
        nitrogen_fraction_entry = tk.Entry(self, textvariable=nitrogen_fraction_value)
        nitrogen_fraction_entry.grid(row=16, column=1, sticky='w')

        # Hydrogen sulfide mole fraction input
        hydrogen_sulfide_fraction_value = tk.DoubleVar()
        hydrogen_sulfide_fraction_label = tk.Label(self, text='Hydrogen sulfide [H2S]:', font=('bold', 12))
        hydrogen_sulfide_fraction_label.grid(row=17, column=0, sticky='w', padx=10, pady=5)
        hydrogen_sulfide_fraction_entry = tk.Entry(self, textvariable=hydrogen_sulfide_fraction_value)
        hydrogen_sulfide_fraction_entry.grid(row=17, column=1, sticky='w')

        # Helium mole fraction input
        helium_fraction_value = tk.DoubleVar()
        helium_fraction_label = tk.Label(self, text='Helium [He]:', font=('bold', 12))
        helium_fraction_label.grid(row=18, column=0, sticky='w', padx=10, pady=5)
        helium_fraction_entry = tk.Entry(self, textvariable=helium_fraction_value)
        helium_fraction_entry.grid(row=18, column=1, sticky='w')

        # Water fraction input
        water_fraction_value = tk.DoubleVar()
        water_fraction_label = tk.Label(self, text='Water [H2O]:', font=('bold', 12))
        water_fraction_label.grid(row=19, column=0, sticky='w', padx=10, pady=5)
        water_fraction_entry = tk.Entry(self, textvariable=water_fraction_value)
        water_fraction_entry.grid(row=19, column=1, sticky='w')

        # Oxygen fraction input
        oxygen_fraction_value = tk.DoubleVar()
        oxygen_fraction_label = tk.Label(self, text='Oxygen [O2]:', font=('bold', 12))
        oxygen_fraction_label.grid(row=20, column=0, sticky='w', padx=10, pady=5)
        oxygen_fraction_entry = tk.Entry(self, textvariable=oxygen_fraction_value)
        oxygen_fraction_entry.grid(row=20, column=1, sticky='w')

        # Argon fraction input
        argon_fraction_value = tk.DoubleVar()
        argon_fraction_label = tk.Label(self, text='Argon [Ar]:', font=('bold', 12))
        argon_fraction_label.grid(row=21, column=0, sticky='w', padx=10, pady=5)
        argon_fraction_entry = tk.Entry(self, textvariable=argon_fraction_value)
        argon_fraction_entry.grid(row=21, column=1, sticky='w')

        # Hydrogen fraction input
        hydrogen_fraction_value = tk.DoubleVar()
        hydrogen_fraction_label = tk.Label(self, text='Hydrogen [H2]:', font=('bold', 12))
        hydrogen_fraction_label.grid(row=22, column=0, sticky='w', padx=10, pady=5)
        hydrogen_fraction_entry = tk.Entry(self, textvariable=hydrogen_fraction_value)
        hydrogen_fraction_entry.grid(row=22, column=1, sticky='w')

        # Carbon monoxide fraction input
        carbon_monoxide_fraction_value = tk.DoubleVar()
        carbon_monoxide_fraction_label = tk.Label(self, text='Carbon monoxide [CO]:', font=('bold', 12))
        carbon_monoxide_fraction_label.grid(row=23, column=0, sticky='w', padx=10, pady=5)
        carbon_monoxide_fraction_entry = tk.Entry(self, textvariable=carbon_monoxide_fraction_value)
        carbon_monoxide_fraction_entry.grid(row=23, column=1, sticky='w')

        status_label = tk.Label(self, text="", fg="red",  bg="#FFFDD0", wraplength=400, relief="solid", bd=2, font=('bold', 14))
        status_label.grid(row=24, rowspan=2, column=0, columnspan=3, sticky="w",  padx=5, pady=5, ipadx=5, ipady=5)
        status_label.config(text="Mole fractions must sum to 1")

        # Make working array arranged in correct order from inputs for Xi fractions

        XI = []
        T = float
        p = float

        def take_input_values():
            def get_float(entry):
                try:
                    return float(entry.get())
                except ValueError:
                    return 0.0
            XI = np.array([])
            XI = np.append(XI, [get_float(methane_fraction_entry),
                                get_float(nitrogen_fraction_entry),
                                get_float(carbon_dioxide_fraction_entry),
                                get_float(ethane_fraction_entry),
                                get_float(propane_fraction_entry),
                                get_float(water_fraction_entry),
                                get_float(hydrogen_sulfide_fraction_entry),
                                get_float(hydrogen_fraction_entry),
                                get_float(carbon_monoxide_fraction_entry),
                                get_float(oxygen_fraction_entry),
                                get_float(isoButane_fraction_entry),
                                get_float(nButane_fraction_entry),
                                get_float(isoPentane_fraction_entry),
                                get_float(nPentane_fraction_entry),
                                get_float(nHexane_fraction_entry),
                                get_float(nHeptane_fraction_entry),
                                get_float(nOctane_fraction_entry),
                                get_float(nNonane_fraction_entry),
                                get_float(nDecane_fraction_entry),
                                get_float(helium_fraction_entry),
                                get_float(argon_fraction_entry),
                                ])
            T = get_float(temperature_entry)
            p = get_float(pressure_entry)

            return XI, T, p

        ready_for_Z = False

        def on_submit():
            self.XI, self.T, self.p = take_input_values()

            if self.T == 0 or self.p == 0:
                status_label.config(
                    text="Error: Temperature and Pressure values must be entered correctly",
                    fg="red"
                )
                return

            # Check sum 

            total = np.sum(self.XI[0:21])
            
            if total == 0:
                status_label.config(
                    text="Error: Total mole fraction is 0. Cannot normalize.",
                    fg="red"
                )
                return 
            if total != 1:
                status_label.config(text=f"Warning: Gas composition does not sum to 1. Total mole fraction = {total:.3f}. Proceeding may reduce calculation accuracy.", fg="red")
            else:
                status_label.config(text=f"Total mole fraction = {total:.3f}", fg="green")

            # Normalize components
      
            self.XI = np.array(self.XI)
            self.XI[0:21] = self.XI[0:21]/total
            
            self.ready_for_Z = True
            calculate_Z_button.grid()

        def calculate_Z():
            self.D = 0.0
            self.Rgas = 8.31451e-3
            self.Tol = 0.5e-9
            
            p = self.p
            T = self.T
            Rgas = self.Rgas
            Tol = self.Tol
            D = self.D

            # M calculation
            M = np.sum(self.XI[0:21]*self.Mi[0:21])

            # B - factor calculation
            B = 0.0
            for n in range(18):
                for i in range(21):
                    for j in range(21):
                        Xij = self.XI[i]*self.XI[j] 
                        if Xij != 0:
                            Gprim = self.Gijast.get(i+1, {}).get(j+1, 1.0)
                            Gij = (Gprim * ((self.Gi[i] + self.Gi[j]) / 2 + 1 - self.gn[n]))**self.gn[n]
                            Qpart = (self.Qi[i] * self.Qi[j] + 1 - self.qn[n])**self.qn[n]
                            Fpart = (self.Fi[i]**0.5 * self.Fi[j]**0.5 + 1 - self.fn[n])**self.fn[n]
                            Spart = (self.Si[i] * self.Si[j] + 1 - self.sn[n])**self.sn[n]
                            Wpart = (self.Wi[i] * self.Wi[j] + 1 - self.wn[n])**self.wn[n]
                            Eprim = self.Eijast.get(i+1,{}).get(j+1, 1.0)
                            Eij = (Eprim * (self.Ei[i] * self.Ei[j])**0.5 )
                            B += self.an[n] * self.T**(-self.un[n]) * Xij * Eij**(self.un[n]) * (self.Ki[i] * self.Ki[j])**1.5 * Gij * Qpart * Fpart * Spart * Wpart

            # U - factor calculation
            U = 0.0
            Upart2 = 0.0
            for i in range(20):
                for j in range(i+1, 21):
                    Uij_val = self.Uij.get(i+1, {}).get(j+1, 1.0)
                    Upart2 += self.XI[i] * self.XI[j] * (Uij_val**5 - 1) * (self.Ei[i] * self.Ei[j])**2.5
            Upart1 = np.sum(self.XI[0:21] * self.Ei[0:21]**2.5)
            U = (Upart1**2 + 2 * Upart2 )**0.2
    
            # G - factor calculation
            G = 0.0
            Gpart2 = 0.0
            for i in range(20):
                for j in range(i+1, 21):
                    Gprim = self.Gijast.get(i+1, {}).get(j+1, 1.0)
                    Gpart2 += self.XI[i] * self.XI[j] * (Gprim * (self.Gi[i] + self.Gi[j]) / 2 - 1) * (self.Gi[i] + self.Gi[j])
            Gpart1 = np.sum(self.XI[0:21] * self.Gi[0:21])
            G = Gpart1 + 2 * Gpart2

            # Q - factor calculation
            Q = np.sum(self.XI[0:21] * self.Qi[0:21])

            # F - factor calculation
            F = np.sum(self.XI[0:21]**2 * self.Fi[0:21])

            # K - factor calculation
            K = 0.0
            Kpart2 = 0.0
            for i in range(20):
                for j in range(i+1, 21):
                    Kij_val = self.Kij.get(i+1,{}).get(j+1, 1.0)
                    Kpart2 += self.XI[i] * self.XI[j] * (Kij_val**5 - 1) * (self.Ki[i] * self.Ki[j])**2.5
            Kpart1 = np.sum(self.XI[0:21] * self.Ki[0:21]**2.5)
            K = (Kpart1**2 + 2 * Kpart2 )**0.2
      
            # Begin iterations
            max_iterations = 100
            rho1 = 0.000001
            rho2 = 40.0

            def pressure_function(rho):
                rho = float(rho)
                
                rhor = rho * K**3

                sum_part1 = np.sum(self.an[12:18] * (G + 1 - self.gn[12:18])**self.gn[12:18] * (Q**2 + 1 - self.qn[12:18])**self.qn[12:18] * 
                            (F + 1 - self.fn[12:18])**self.fn[12:18] * U**self.un[12:18] * T**(-self.un[12:18]))
                
                sum_part2 = np.sum(self.an[12:58] * (G + 1 - self.gn[12:58])**self.gn[12:58] * (Q**2 + 1 - self.qn[12:58])**self.qn[12:58] * 
                            (F + 1 - self.fn[12:58])**self.fn[12:58] * U**self.un[12:58] * T**(-self.un[12:58]) * (self.bn[12:58] - self.cn[12:58] * self.kn[12:58] * rhor**self.kn[12:58]) * 
                            rhor**self.bn[12:58] * np.exp(-self.cn[12:58] * rhor**self.kn[12:58]))

                pr = rho * Rgas * T * (1 + B * rho - rhor * sum_part1 + sum_part2)            
                return pr

            for i in range(max_iterations):
                 
                p1 = pressure_function(rho1) - p
                p2 = pressure_function(rho2) - p
                
                if i == 0:
                    if (p1 * p2) >= 0:
                        print("Boundary inputs do not bracket the root of function!")
                        return

                if p2 - p1 == 0:
                    print("Warning: p2 - p1 is zero. Cannot compute rho3 safely.")
                    rho3 = (rho1 + rho2) / 2 
                else:
                    rho3 = rho1 - p1 * (rho2 - rho1) / (p2 - p1)

                p3 = pressure_function(rho3) - p

                denom1 = (p1 - p2) * (p1 - p3)
                denom2 = (p2 - p1) * (p2 - p3)
                denom3 = (p3 - p1) * (p3 - p2)

                D = (rho1 * p2 * p3 / denom1) + (rho2 * p1 * p3 / denom2) + (rho3 * p1 * p2 / denom3)

                if (D - rho1) * (D - rho2) >= 0:
                    D = (rho1 + rho2) / 2

                pfinal = pressure_function(D) - p

                if (abs(p3) < abs(pfinal)) and (pfinal * p3) > 0:
                    if (p3 * p1) > 0:
                        rho1 = rho3
                    else:
                        rho2 = rho3
                else:
                    if (pfinal * p3) < 0:
                        rho1 = D
                        rho2 = rho3
                    elif (p3 * p1):
                        rho1 = D
                    else:
                        rho2 = D

                if abs(pfinal) <= Tol:
                    global Z
                    Z = p / (D * Rgas * T)
                    Z_factor_result.config(text=f"  Z = {Z}")
                    save_to_archive(self.XI, T, p, Z)
                    return 

        Z_factor_result = tk.Label(self, text="  Z = ", bd=4, relief="ridge", width=25, height=2, font=('bold', 13), anchor="w")
        Z_factor_result.grid(row=21, column=2, rowspan=2, columnspan=2, ipadx=5)

        calculate_Z_button = tk.Button(self, text="Calculate Z", command=calculate_Z)
        calculate_Z_button.grid(row=23, column=3, sticky="e")
        calculate_Z_button.grid_remove() 

        submit_button = tk.Button(self, text="Submit entry values", command=on_submit)
        submit_button.grid(row=23, column=2, sticky="W", padx=10)