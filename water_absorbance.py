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
        """Displaying multiple Plots in one Figure

        Args:
            label (list): Label of the plots

        Returns:
            Plot: showing multiple Plots in one Figure
        """
        filename = self.file_path
        fig = plt.figure()
        ax = fig.add_subplot()

        rawdata = self.dptfile_to_dict()
        for i in range(len(filename)):

            plt.plot(rawdata[i]["wavenumber"], rawdata[i]["extinction"], label = label[i])# + " Max SI: " + str(round(max(rawdata[i]["extinction"][1633:1936]),3)))

        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
        #ax.set(xlim=(1560, 1710), ylim = (-0.1, 0.1))
        plt.legend()
        plt.show()

        return None


    def fit_gaussian_to_data(self, label, interval):
        """_summary_

        Args:
            label (list): labelname for each graph in a list. name in str
            interval (indexarea): example -> 1633:1936 

        Returns:
            List: wavenumber from each maxima and graph in that choosen area
        """
        filename = self.file_path

        #Define the Gaussian function 
        def Gauss(x,a,x0,sigma):
            return a * np.exp(-(x-x0)**2/(2*sigma**2))

        rawdata = self.dptfile_to_dict()
        wavenumber = [rawdata[i]["wavenumber"][interval] for i in range(len(filename))]
        extinction = [rawdata[i]["extinction"][interval] for i in range(len(filename))]

        fig = plt.figure()
        ax = fig.add_subplot()

        max_extinction = []
        for i in range(len(filename)):
            x = np.asarray(wavenumber[i])
            y = np.asarray(extinction[i])
                        
            mean = x[y.argmax()]                  #note this correction
            sigma = mean - np.where(y > y.max() * np.exp(-.5))[0][0]        #note this correction

            popt, pcov = curve_fit(Gauss, x, y, p0 = [y.max(),mean,sigma])

            # best fit
            fit_extinction = y[extinction[i].index(y.max())-30 : extinction[i].index(y.max())+30]
            fit_wavenumber = x[extinction[i].index(y.max())-30 : extinction[i].index(y.max())+30]

            shorten_popt, shorten_cov = curve_fit(Gauss, fit_wavenumber, fit_extinction, p0 = [*popt])

            

            x_value = np.linspace(x.min(), x.max(),1000)
            raw_text = label[i] + " raw data"
            plt.plot(x, y, label = raw_text)
            text = label[i] + " fit"
            plt.plot(x_value, Gauss(x_value, *shorten_popt), label = text, linestyle= 'dotted')
            max_extinction.append(Gauss(x_value, *shorten_popt).max())

        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        ax.set(xlim=(min(wavenumber[0]), max(wavenumber[0])))
        plt.grid(axis='both', color='0.95')
        plt.legend()
        plt.show()

        return max_extinction


    def standardcurve(self, label):
        extinction = self.fit_gaussian_to_data(label)
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

        text = "y = (" + str(round(popt[0],3)) + "$\pm$" + str(round(std[0], 4)) +") x \n"+ "R$^2$ = " + str(round(r_square, 3))

        plt.scatter(conc, extinction, marker = ".")
        plt.plot(x_value, linear(x_value, popt[0]), label = text, color = "orange")
        plt.xlabel("c(Citric Acid) in mM", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
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

allwater = Raw_data_Plot(["raw_data/H2O.DPT", "raw_data/D2O.DPT", "raw_data/2mMCitronensaeure.DPT"])
print(allwater.fit_gaussian_to_data(["H$_2$O", "D$_2$O", "2mM Citric Acid"]))
print(allwater.multiple_plot_show(["H$_2$O", "D$_2$O", "2mM Citric Acid"]))
sample = Raw_data_Plot(["raw_data/probe1_Citronensaeure.DPT", "raw_data/probe2_Citronensaeure.DPT"]).fit_gaussian_to_data(["Sample 1", "Sample 2"])
print(sample)
print(citricacid_standardcurve.fit_gaussian_to_data(["0.2 mM ",
            "0.4 mM ",
           "0.6 mM ",
          "0.8 mM ",
           "1.0 mM ",
            "1.2 mM ",
            "1.4 mM ",
            "1.6 mM ",
            "1.8 mM ",
            "2.0 mM ",
        ]))



