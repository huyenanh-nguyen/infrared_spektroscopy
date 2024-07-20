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
            plt.plot(wavenumber[i], extinction[i], label = label[i])
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


beta = Plotting(["raw_data/50grad_Polylysin.DPT"])
print(beta.plotpoints(["bei 50°C"], 1680,  1821))
