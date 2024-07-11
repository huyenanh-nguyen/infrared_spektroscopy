import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from pathlib import Path

class Raw_data_Plot:
    def __init__(self, filepath):
        """passing filepath as string

        Args:
            filepath (list): filepath as string in list -> like ["raw_data/H2O.DPT"]
        """
        self.file_path = filepath

    def dptfile_to_dict(self):
        """convert the dptfile that has the Bruker datapoints of the measurements to a dictfile.

        Returns:
            dict: return in the format
                    >>> {0: {"wavenumber": [...], "extinction" : [...]}, 
                         1: {"wavenumber": [...], "extinction" : [...]}}
        """
        filename = self.file_path

        
        data_dict = {}
        for i in range(len(filename)):
            wavenumber = []
            extinction = []
            with open(filename[i]) as file:
                for line in file:
                    (wave, extin) = line.split(",")
    
                    wavenumber.append(float(wave))
                    extinction.append(float(extin))
                
            data_dict[i] = {"wavenumber" : wavenumber, "extinction" : extinction}
        return data_dict
    
    def single_plot_show(self):
        """plotting the data individually
        """
        filename = self.file_path
        for i in range(len(filename)):
            wavenumber = self.dptfile_to_dict()[i]["wavenumber"]
            extinction = self.dptfile_to_dict()[i]["extinction"]

            fig = plt.figure()
            ax = fig.add_subplot()
            plt.plot(wavenumber, extinction)
            plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
            plt.ylabel("SI", fontsize=12)
            plt.grid(axis='x', color='0.95')
            plt.show()

        return None
    

    def multiple_plot_show(self, label:list):
        filename = self.file_path
        fig = plt.figure()
        ax = fig.add_subplot()

        rawdata = self.dptfile_to_dict()
        for i in range(len(filename)):

            plt.plot(rawdata[i]["wavenumber"], rawdata[i]["extinction"], label = label[i] + " Max SI: " + str(round(max(rawdata[i]["extinction"][1633:1936]),3)))

        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='x', color='0.95')
        ax.set(xlim=(1550, 1900), ylim = (-0.02, 0.35))
        plt.legend()
        plt.show()

        return None
    
    def standardcurve(self):
        filename = self.file_path
        fig = plt.figure()
        ax = fig.add_subplot()

        rawdata = self.dptfile_to_dict()
        conc = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        extinction = []
        for i in range(len(filename)):
            extinction.append(max(rawdata[i]["extinction"][1633:1936]))
        
        def linear(x, a):
            return a * x
        
        x_value = np.linspace(0, 2.1)

        popt, pcov = curve_fit(linear, conc, extinction)
        residuals = extinction - linear(np.asarray(conc), *popt)
        ss_res = np.sum(residuals ** 2)
        ss_total = np.sum((extinction - np.mean(extinction)) ** 2)
        r_square = 1 - (ss_res / ss_total)

        std = np.sqrt(np.diag(pcov))

        text = "y = (" + str(round(popt[0],3)) + "$\pm$" + str(round(std[0], 4)) +") x \n"+ "R$^2$ = " + str(round(r_square, 3))

        plt.scatter(conc, extinction, marker = ".")
        plt.plot(x_value, linear(x_value, popt[0]), label = text, color = "orange")
        plt.xlabel("c(Citric Acid) in mM", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.legend()
        plt.show()
        return None

    



water = Raw_data_Plot(["raw_data/H2O.DPT"])
water_D = Raw_data_Plot(["raw_data/D2O.DPT"])
citricacid_standardcurve = Raw_data_Plot([
            "raw_data/02mMCitronensaeure.DPT",
            "raw_data/04mMCitronensaeure.DPT",
            "raw_data/06mMCitronensaeure.DPT",
            "raw_data/08mMCitronensaeure.DPT",
            "raw_data/1mMCitronensaeure.DPT",
            "raw_data/12mMCitronensaeure.DPT",
            "raw_data/14mMCitronensaeure.DPT",
            "raw_data/16mMCitronensaeure.DPT",
            "raw_data/18mMCitronensaeure.DPT",
            "raw_data/2mMCitronensaeure.DPT",
                                          ])


""" print(citricacid_standardcurve.multiple_plot_show([
        "0.2 mM",
        "0.4 mM",
       "0.6 mM",
       "0.8 mM",
        "1.0 mM",
        "1.2 mM",
        "1.4 mM",
        "1.6 mM",
        "1.8 mM",
        "2.0 mM",
    ])) """
print("______")
print(citricacid_standardcurve.standardcurve())
