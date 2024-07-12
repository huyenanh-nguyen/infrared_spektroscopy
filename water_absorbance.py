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


    def fit_gaussian_to_data(self):
        filename = self.file_path

        #Define the Gaussian function 

        def Gauss(x,a,x0,sigma):
            return a * np.exp(-(x-x0)**2/(2*sigma**2))

        
        rawdata = self.dptfile_to_dict()
        wavenumber = []
        extinction = []
        for k in range(len(filename)):
            wavenumber.append(rawdata[k]["wavenumber"][1633:1936])
            extinction.append(rawdata[k]["extinction"][1633:1936])

        fig = plt.figure()
        ax = fig.add_subplot()

        for i in range(len(filename)):
            x = np.asarray(wavenumber[i])
            y = np.asarray(extinction[i])
                        
            mean = x[y.argmax()]                  #note this correction
            sigma = mean - np.where(y > y.max() * np.exp(-.5))[0][0]        #note this correction

            popt, pcov = curve_fit(Gauss, x, y, p0 = [y.max(),mean,sigma])
            # residuals = y - Gauss(x, *popt)
            # ss_res = np.sum(residuals ** 2)
            # ss_total = np.sum((y - np.mean(y)) ** 2)
            # r_square = 1 - (ss_res / ss_total)

            # std = np.sqrt(np.diag(pcov))

            # best fit
            pop,cot = curve_fit(Gauss, x, y, p0 = [*popt])

            max_extinction = y[extinction[i].index(y.max())-30 : extinction[i].index(y.max())+30]
            max_wavenumber = x[extinction[i].index(y.max())-30 : extinction[i].index(y.max())+30]
            # breakpoint()
            sh_mean = max_wavenumber[max_extinction.argmax()]
            sh_sigma = max_wavenumber[max_extinction.argmax()] - np.where(max_extinction > max_extinction.max() * np.exp(-.5))[0][0]

            shorten_popt, shorten_cov = curve_fit(Gauss, max_wavenumber, max_extinction, p0 = [*pop])

            x_value = np.linspace(x.min(), x.max())
            plt.plot(x, y, label = "Raw Data")
            # plt.plot(x_value, Gauss(x_value, *popt), label = "fit", linestyle= "dotted") # fit
            # plt.plot(x_value, Gauss(x_value, *pop), label = "best fit", linestyle= 'dashdot')
            plt.plot(x_value, Gauss(x_value, *shorten_popt), label = "shorten fit", linestyle= 'dashed')

        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        ax.set(xlim=(min(wavenumber[0]), max(wavenumber[0])))
        plt.legend()
        plt.show()
        return None

    
    def fit_two_exponential(self):
        filename = self.file_path

        #Define the Gaussian function 

        def Expo(x, a, b, p, q):
            return 1 + ( 1 / (a * np.exp(p * x) + b * np.exp( -q * x)))

        
        rawdata = self.dptfile_to_dict()
        wavenumber = []
        extinction = []
        for k in range(len(filename)):
            wavenumber.append(np.asarray(rawdata[k]["wavenumber"][1633:1936]))
            extinction.append(np.asarray(rawdata[k]["extinction"][1633:1936]))

        fig = plt.figure()
        ax = fig.add_subplot()

        for i in range(len(filename)):
            x = wavenumber[i]
            y = extinction[i]
                        
            mean = x[y.argmax()]                  #note this correction
            sigma = mean - np.where(y > y.max() * np.exp(-.5))[0][0]        #note this correction

            popt, pcov = curve_fit(Expo, x, y, p0 = [y.min(), y.max()*10 ,0.5, 1.75])
            residuals = y - Expo(x, *popt)
            ss_res = np.sum(residuals ** 2)
            ss_total = np.sum((y - np.mean(y)) ** 2)
            r_square = 1 - (ss_res / ss_total)

            std = np.sqrt(np.diag(pcov))

            x_value = np.linspace(x.min(), x.max())
            # breakpoint()
            plt.plot(x, y, label = "Raw Data")
            plt.plot(x_value, Expo(x_value, *popt), label = "fit", linestyle= "dotted") # fit

        plt.xlabel("wavenumber in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        ax.set(xlim=(min(wavenumber[0]), max(wavenumber[0])))
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
print(allwater.fit_gaussian_to_data())
# print(allwater.fit_two_exponential())
print([        "0.2 mM, ",
            "0.4 mM, ",
           "0.6 mM, ",
          "0.8 mM, ",
           "1.0 mM, ",
            "1.2 mM, ",
            "1.4 mM, ",
            "1.6 mM, ",
            "1.8 mM, ",
            "2.0 mM, ",
        ])
# print("______")
# print(citricacid_standardcurve.standardcurve())


