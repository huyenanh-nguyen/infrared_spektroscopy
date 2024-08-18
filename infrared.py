import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


class Plotting:
    def __init__(self, filepath):
        """passing filepath as string

        Args:
            filepath (list): filepath as string in list -> like ["/users/.../raw_data/H2O.DPT", "/users/.../raw_data/D2O.DPT"] (if the file is not in the same level like this pythonfile)
                             if the files are in a folder called "ra_data", that is in the same directory, then [raw_data/H2O.DPT, raw_data/D2O.DPT] is enough
        """
        self.file_path = filepath
    
    def dptfile_to_dict(self):
        """
        extracting the wave and extinction data from the .dpt-file.
        This file format is a special form from the IFS 66v/S Spektrometer (Bruker, Berlin - Humboldt Universität zu Berlin am Biophysik Institut) instrument.

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
    

    def full_spectros(self, label):
        """plotting all points from the file

        Args:
            label (list): label for each sample -> as string

        Returns:
            Plt: Plot
        """
        wavenumber = [self.dptfile_to_dict()[i]["wavenumber"] for i in range(len(self.dptfile_to_dict()))]
        extinction = [self.dptfile_to_dict()[i]["extinction"] for i in range(len(self.dptfile_to_dict()))]

        plt.figure()
        plt.xlabel("Wellenzahl in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95') 

        for i in range(len(wavenumber)):
            plt.plot(wavenumber[i], extinction[i], label = label[i])
        
        plt.legend()
        plt.show()

        return None
    
    
    def from_to_spectro(self, label, indexstart, indexstop):
        """
        plotting points from a specific interval

        Args:
            label (list): label for each sample -> as string
            indexstart (int): choose the index where the plot should start
            indexstop (int): choose the index where the plot should end

        Returns:
            Dict: Dictionary contains the max value of the extinction with the depending wavenumber for each sample
        """

        wavenumber = [self.dptfile_to_dict()[i]["wavenumber"][indexstart:indexstop] for i in range(len(self.dptfile_to_dict()))]
        extinction = [self.dptfile_to_dict()[i]["extinction"][indexstart:indexstop] for i in range(len(self.dptfile_to_dict()))]

        plt.figure()
        plt.xlabel("Wellenzahl in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
        

        max_stuff = {}
        for i in range(len(wavenumber)):
            max_stuff[label[i]] = {"max_ext" : max(extinction[i]), "max_wave" : wavenumber[i][extinction[i].index(max(extinction[i]))]}
            plt.plot(wavenumber[i], extinction[i], label = label[i])
        

        plt.legend()
        plt.show()

        return max_stuff
    

    def SI_standardcurve(self, label, indexstart, indexstop):
        """
        fit the linear regression

        Args:
            label (list): label for each sample -> as string
            indexstart (int): choose the index where the plot should start
            indexstop (int): choose the index where the plot should end

        Returns:
            list: with max Extinction within the interval
        """
        wavenumber = [self.dptfile_to_dict()[i]["wavenumber"][indexstart:indexstop] for i in range(len(self.dptfile_to_dict()))]
        extinction = [self.dptfile_to_dict()[i]["extinction"][indexstart:indexstop] for i in range(len(self.dptfile_to_dict()))]


        plt.figure()
        plt.xlabel("Wellenzahl in cm$^{-1}$", fontsize=12)
        plt.ylabel("SI", fontsize=12)
        plt.grid(axis='both', color='0.95')
        

        max_ext = []
        for i in range(len(wavenumber)):
            text = label[i] + ", SI$_{max}$ : " + f"{max(extinction[i]) : .3f}"
            max_ext.append(max(extinction[i]))
            plt.plot(wavenumber[i], extinction[i], label = text)
        

        plt.legend()
        plt.show()

        return max_ext
    

    def standardcurve(self, label, indexstart, indexstop):
        """_summary_

        Args:
            label (_type_): _description_
            indexstart (_type_): _description_
            indexstop (_type_): _description_

        Returns:
            _type_: _description_
        """

        extinction = np.array(self.SI_standardcurve(label, indexstart, indexstop))
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
        y_label = "SI$_{max}$ bei  1720 cm$^{-1}$"
        plt.ylabel(y_label, fontsize=12)
        plt.grid(axis='both', color='0.95')
        plt.legend()
        plt.show()
        return None


#___________________________________________________________________________________________________________________________________________________________________
# standard curve 
standard_path = [
            "raw_data/02mMCitronensaeure.DPT",
            "raw_data/04mMCitronensaeure.DPT",
            "raw_data/06mMCitronensaeure.DPT",
            "raw_data/08mMCitronensaeure.DPT",
            "raw_data/1mMCitronensaeure.DPT",
            "raw_data/12mMCitronensaeure.DPT",
            "raw_data/14mMCitronensaeure.DPT",
            "raw_data/16mMCitronensaeure.DPT",
            "raw_data/18mMCitronensaeure.DPT",
            "raw_data/2mMCitronensaeure.DPT"    
]

standard_start = 1633
standard_stop = 1936

standard_label = [
            "0.2 mM ",
            "0.4 mM ",
            "0.6 mM ",
            "0.8 mM ",
            "1.0 mM ",
            "1.2 mM ",
            "1.4 mM ",
            "1.6 mM ",
            "1.8 mM ",
            "2.0 mM "
]

citronacidsample_path = ["raw_data/probe1_Citronensaeure.DPT", "raw_data/probe2_Citronensaeure.DPT"]

standard = Plotting(standard_path)
print(standard.standardcurve(standard_label, standard_start, standard_stop))

# unkownsample

unknown_path = [
    "raw_data/probe1_Citronensaeure.DPT",
    "raw_data/probe2_Citronensaeure.DPT"
]

unknown_label = [
    "Sample 1",
    "Sample 2"
]
sample_start = 1633
sample_stop = 1936

unkown = Plotting(unknown_path)
print(unkown.SI_standardcurve(unknown_label, sample_start, sample_stop))

# citric adcid and h2o spectrum
h2o_acid_path = ["raw_data/H2O.DPT", "raw_data/2mMCitronensaeure.DPT" ]
h2o_acid_label = ["H$_2$O", "2.0 mM Citronensäure"]
h2o_acid = Plotting(h2o_acid_path)

print(h2o_acid.full_spectros(h2o_acid_label))


# h2o
h2o_path = ["raw_data/H2O.DPT"]
h2o_label = ["H$_2$O"]

h2o = Plotting(h2o_path)
print(h2o.full_spectros(h2o_label))



# polylysin
polylysin_path = [
    "raw_data/pHneutral_Polylysin.DPT",
    "raw_data/pH1162_Polylysin.DPT",
    "raw_data/50grad_Polylysin.DPT",
]

polylysin_label = [
    "Polylysin neutral", "Polylysin basisch", "Polylysin 50°C"
]

polylysin_start = 1750
polylysin_stop = 1900

polylysin = Plotting(polylysin_path)
print(polylysin.from_to_spectro(polylysin_label, polylysin_start, polylysin_stop))


# polylysin and water

poly_h2o_path = [
    "raw_data/pHneutral_Polylysin.DPT",
    "raw_data/pH1162_Polylysin.DPT",
    "raw_data/50grad_Polylysin.DPT",
    "raw_data/H2O.DPT"
]

poly_h2o_label = [
    "Polylysin neutral", "Polylysin basisch", "Polylysin 50°C", "H$_2$O"
]


poly_h2o = Plotting(poly_h2o_path)
print(poly_h2o.from_to_spectro(poly_h2o_label, polylysin_start, polylysin_stop))