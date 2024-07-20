import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


class Plotting:
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


    def plotpoints(self, label, xstart, xstop):
        wavenumber = [self.dptfile_to_dict()[i]["wavenumber"][xstart : xstop] for i in range(len(self.dptfile_to_dict()))]
        extinction = [self.dptfile_to_dict()[i]["extinction"][xstart : xstop] for i in range(len(self.dptfile_to_dict()))]

        fig = plt.figure()
        ax =  fig.add_subplot()

        max_SI = []
        max_wave = []

        for i in range(len(wavenumber)):
            labeling = label[i] + ", SI$_{max}$: " + f"{max(extinction[i]):0.3f}"
            plt.plot(wavenumber[i], extinction[i], label = labeling)
            max_SI.append(max(extinction[i]))
            # breakpoint()
            max_wave.append(wavenumber[i][extinction[i].index(max(extinction[i]))])
                            
        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
        ax.set(xlim=(min(wavenumber[0]), max(wavenumber[0])))
        plt.legend()
        plt.show()

        return max_SI, max_wave
    


    def standardcurve(self, label, start, end):
        extinction = self.plotpoints(label, start, end)[0]
        fig = plt.figure()
        conc = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        
        def linear(x, a):
            return a * x
        
        x_value = np.linspace(0, 2.1)

        popt, pcov = curve_fit(linear, conc, extinction)
        residuals = extinction - linear(np.asarray(conc), *popt)
        ss_res = np.sum(residuals ** 2)
        ss_total = np.sum((extinction - np.mean(extinction)) ** 2)
        r_square = 1 - (ss_res / ss_total)

        std = np.sqrt(np.diag(pcov))

        text = "y = (" + f"{popt[0] : 0.3f}" + "$\pm$" + f"{std[0] : 0.3f}" +") x \n"+ "R$^2$ = " + f"{r_square : 0.3f}"

        plt.scatter(conc, extinction, marker = ".")
        plt.plot(x_value, linear(x_value, popt[0]), label = text, color = "orange")
        plt.xlabel("c(Citric Acid) in mM", fontsize=12)
        # plt.errorbar(conc, extinction[0], yerr= extinction[-1], fmt=' ', capsize=3, color = "dimgrey")
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
        plt.legend()
        plt.show()
        return None


beta = Plotting(["raw_data/50grad_Polylysin.DPT"])
print(beta.plotpoints(["bei 50°C"], 1680,  1821))



citricacid_standardcurve = Plotting([
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

print(citricacid_standardcurve.standardcurve(["0.2 mM ",
            "0.4 mM ",
           "0.6 mM ",
          "0.8 mM ",
           "1.0 mM ",
            "1.2 mM ",
            "1.4 mM ",
            "1.6 mM ",
            "1.8 mM ",
            "2.0 mM ",
        ],1633, 1936))