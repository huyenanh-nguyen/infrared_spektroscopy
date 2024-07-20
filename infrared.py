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
            labeling = label[i] #+ ", SI$_{max}$: " + f"{max(extinction[i]):0.3f}"
            plt.plot(wavenumber[i], extinction[i], label = labeling)
            max_SI.append(max(extinction[i]))
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
        plt.ylabel("SI$_{max}$", fontsize=12)
        plt.grid(axis='both', color='0.95')
        plt.legend()
        plt.show()
        return None




sample = Plotting(["raw_data/probe1_Citronensaeure.DPT", "raw_data/probe2_Citronensaeure.DPT"]).plotpoints(["Sample 1", "Sample 2"], 1633, 1936)
# print(sample)
beta = Plotting(["raw_data/50grad_Polylysin.DPT"])
# print(beta.plotpoints(["bei 50°C"], 1680,  1821))



citricacid_standardcurve = Plotting([
            # "raw_data/02mMCitronensaeure.DPT",
            # "raw_data/04mMCitronensaeure.DPT",
            # "raw_data/06mMCitronensaeure.DPT",
            # "raw_data/08mMCitronensaeure.DPT",
            # "raw_data/1mMCitronensaeure.DPT",
            # "raw_data/12mMCitronensaeure.DPT",
            # "raw_data/14mMCitronensaeure.DPT",
            # "raw_data/16mMCitronensaeure.DPT",
            # "raw_data/18mMCitronensaeure.DPT",
            "raw_data/2mMCitronensaeure.DPT",
                                          ])

print(citricacid_standardcurve.plotpoints([
        #     "0.2 mM ",
        #     "0.4 mM ",
        #    "0.6 mM ",
        #   "0.8 mM ",
        #    "1.0 mM ",
        #     "1.2 mM ",
        #     "1.4 mM ",
        #     "1.6 mM ",
        #     "1.8 mM ",
            "2.0 mM ",
        ],0, 2412))

full_cic_was = Plotting(["raw_data/pHneutral_Polylysin.DPT", "raw_data/pH1162_Polylysin.DPT", "raw_data/50grad_Polylysin.DPT", "raw_data/H2O.DPT", "raw_data/D2O.DPT"])
print(full_cic_was.plotpoints(["Polylysin neutral", "Polylysin basisch", "Polylysin 50°C", "H$_2$O", "D$_2$O"], 1750, 1900))
# allwater = Plotting(["raw_data/H2O.DPT", "raw_data/D2O.DPT", "raw_data/2mMCitronensaeure.DPT"])
# print(allwater.plotpoints(["H$_2$O", "D$_2$O",
#             "2.0 mM ",
#         ],1633, 1936))


# searching Extrema of Polylysin
# pH-neutral
# polylysin = Plotting(["raw_data/pHneutral_Polylysin.DPT", "raw_data/pH1162_Polylysin.DPT", "raw_data/50grad_Polylysin.DPT"])

# extinction = [polylysin.dptfile_to_dict()[i]["extinction"] for i in range(3)]
# wavenumber = [polylysin.dptfile_to_dict()[i]["wavenumber"] for i in range(3)]

# max_ext1 = [max(extinction[i][1818:1867]) for i in range(3)]
# max_ext2 = [max(extinction[i][1790:1815]) for i in range(3)]

# max_wave1 = [wavenumber[i][extinction[i].index(max(extinction[i][1818:1867]))] for i in range(3)]
# max_wave2 = [wavenumber[i][extinction[i].index(max(extinction[i][1790:1815]))] for i in range(3)]

# print("1", max_ext1, max_wave1)
# print('__________')
# print("2", max_ext2, max_wave2)
