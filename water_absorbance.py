import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from pathlib import Path


#Define the Gaussian function 
def Gauss(x,a,x0,sigma):
    return a * np.exp(-(x-x0)**2/(2*sigma**2))



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


    def fit_gaussian_to_data(self, label, start, end):
        """_summary_

        Args:
            label (list): labelname for each graph in a list. name in str
            start (int): example -> 1633
            end (int): example -> 1936 

        Returns:
            List: wavenumber from each maxima and graph in that choosen area
        """
        filename = self.file_path

        rawdata = self.dptfile_to_dict()
        wavenumber = [rawdata[i]["wavenumber"][start : end] for i in range(len(filename))]
        extinction = [rawdata[i]["extinction"][start : end] for i in range(len(filename))]

        fig = plt.figure()
        ax = fig.add_subplot()

        max_extinction = []
        a = []
        error = []
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
            a.append(shorten_popt[0])
            error.append(np.linalg.norm(y-Gauss(x, *shorten_popt))) # Frobenius norm
            # breakpoint()

        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        ax.set(xlim=(min(wavenumber[0]), max(wavenumber[0])))
        plt.grid(axis='both', color='0.95')
        plt.legend()
        plt.show()

        return max_extinction, a, error


    def standardcurve(self, label, start, end):
        extinction = self.fit_gaussian_to_data(label, start, end)
        fig = plt.figure()
        conc = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        
        def linear(x, a):
            return a * x
        
        x_value = np.linspace(0, 2.1)

        popt, pcov = curve_fit(linear, conc, extinction[0])
        residuals = extinction[0] - linear(np.asarray(conc), *popt)
        ss_res = np.sum(residuals ** 2)
        ss_total = np.sum((extinction[0] - np.mean(extinction[0])) ** 2)
        r_square = 1 - (ss_res / ss_total)

        std = np.sqrt(np.diag(pcov))

        text = "y = (" + str(round(popt[0],3)) + "$\pm$" + str(round(std[0], 4)) +") x \n"+ "R$^2$ = " + str(round(r_square, 3))

        # plt.scatter(conc, extinction[0], marker = ".")
        plt.plot(x_value, linear(x_value, popt[0]), label = text, color = "orange")
        plt.xlabel("c(Citric Acid) in mM", fontsize=12)
        plt.errorbar(conc, extinction[0], yerr= extinction[-1], fmt=' ', capsize=3, color = "dimgrey")
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
        plt.legend()
        plt.show()
        return None
    

    def gauss_for_secondarystructure(self, label, start, end):
        
        filename = self.file_path

        rawdata = self.dptfile_to_dict()
        wavenumber = [rawdata[i]["wavenumber"][start:end] for i in range(len(filename))]
        extinction = [rawdata[i]["extinction"][start:end] for i in range(len(filename))]

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
            fit_extinction = y[extinction[i].index(y.max())-15 : extinction[i].index(y.max())+15]
            fit_wavenumber = x[extinction[i].index(y.max())-15 : extinction[i].index(y.max())+15]

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
print(allwater.fit_gaussian_to_data(["H$_2$O", "D$_2$O", "2mM Citric Acid"],1633, 1936))
print(allwater.multiple_plot_show(["H$_2$O", "D$_2$O", "2mM Citric Acid"]))
sample = Raw_data_Plot(["raw_data/probe1_Citronensaeure.DPT", "raw_data/probe2_Citronensaeure.DPT"]).fit_gaussian_to_data(["Sample 1", "Sample 2"], 1633, 1936)
print(sample)
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


polylysin = Raw_data_Plot(["raw_data/pHneutral_Polylysin.DPT", "raw_data/pH1162_Polylysin.DPT", "raw_data/50grad_Polylysin.DPT"])
print(polylysin.gauss_for_secondarystructure(["pH neutral, auf Eis", "pH = 11.62", "bei 50°C"], 1712,  1876))

