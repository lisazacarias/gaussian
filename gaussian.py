from gaussian_ui import Ui_MainWindow
from PyQt4 import QtGui, QtCore

from pylab import array, plt, floor, show
from numpy import argsort, power, exp, zeros
from scipy.io import loadmat
from scipy.optimize import curve_fit
from operator import itemgetter
from sys import argv, exit
from matplotlib.colors import ColorConverter

from read_xcor_data import *

class GaussianUi(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)
        
        self.ui.uploadButton.pressed.connect(self.getData)
        self.ui.numFits.activated.connect(self.updateGuess)
        self.ui.guessToEdit.activated.connect(self.updateGuessToEdit)
        self.ui.filterPedestal.stateChanged.connect(self.setUseZeros)
        self.ui.getFitButton.pressed.connect(self.getFit)
        
        self.axes = self.ui.gaussPlot.canvas.fig.add_subplot(111, axisbg='white')
        self.useZeros = not self.ui.filterPedestal.isChecked()
        
        self.changeStates(False)
        self.ui.filterPedestal.setVisible(False)
        
    def updateGuessToEdit(self):
        gauss = int(self.ui.guessToEdit.currentText())
        self.ui.heightSlider.setSliderPosition(1000*self.guess[3 + 3*gauss])
        self.ui.centerSlider.setSliderPosition(self.guess[2 + 3*gauss])
        self.ui.widthSlider.setSliderPosition(self.guess[4 + 3*gauss])
    
    def changeStates(self, isEnabled):
        self.ui.numFits.setEnabled(isEnabled)
        self.ui.filterPedestal.setEnabled(isEnabled)
        self.ui.getFitButton.setEnabled(isEnabled)
        self.ui.centerSlider.setEnabled(isEnabled)
        self.ui.heightSlider.setEnabled(isEnabled)
        self.ui.widthSlider.setEnabled(isEnabled)
        self.ui.guessToEdit.setEnabled(isEnabled)
    
    def setUseZeros(self):
        if self.ui.filterPedestal.isChecked():
            self.useZeros = False
        else:
            self.useZeros = True
        
        try:
            self.cleanAndPlotData()
            self.getGuess()
            self.updateOptions(self.ui.numFits, 1, self.potentialPeaks + 1, True)
        except AttributeError:
            pass
        
    def updateGuess(self):
        self.cleanAndPlotData()
        self.getGuess()
        self.updateOptions(self.ui.guessToEdit, 0, 
                           int(self.ui.numFits.currentText()))
        self.updateGuessToEdit()
    
    def cleanAndPlotData(self):
        self.ui.results.clear()
        self.axes.clear()
        self.axes.plot(self.ampList, '.', marker='o')
        self.ui.gaussPlot.canvas.draw()
    
    def updateOptions(self, comboBox, start, end):
        comboBox.clear()
        for i in xrange(start, end):
            comboBox.addItem(str(i))
    
    def updateNumFitsOptions(self):
        self.ui.numFits.setCurrentIndex(0)
        self.ui.guessToEdit.setCurrentIndex(0)
            
        self.cleanAndPlotData()
        
        self.ui.numFits.clear()
        self.ui.guessToEdit.clear()
        
        for i in xrange(1, self.potentialPeaks + 1):
            self.ui.numFits.addItem(str(i))
            
        for i in xrange(0, self.potentialPeaks):
            self.ui.guessToEdit.addItem(str(i))
        
    def plotFit(self, popt, isGuess):
        self.colors = []
    
        self.ui.results.append("adjustment: " + str(self.totalAdjustment))
                               
        # Print and plot the optmized line fit.
        # Note that popt has the same format as the guess, meaning that the first
        # two parameters are the m and b of the line, respectively
        self.ui.results.append("line: " + "m = " + str(popt[0]) + ", b = " + str(popt[1]))
        color = self.axes.plot([popt[0]*j + popt[1] + self.totalAdjustment
                               for j in self.x], '--')[0].get_color()
        self.colors.append(ColorConverter().to_rgb(color))
        
        # Print and plot the optimized gaussian fit(s)
        # Again, the first two elements were the line, and each gaussian is a
        # subsequent group of 3 elements (hence starting at index 2 and
        # incrementing by 3)
        for i in xrange(2, len(popt), 3):
            self.ui.results.append("gaussian " + str(i//3) + ": center = " + str(popt[i])
                   + ", amplitude = " + str(popt[i+1]) + ", width = "
                   + str(popt[i+2]))
            color = self.axes.plot([gaussian(j, popt[i], popt[i + 2], popt[i + 1])
                    + self.totalAdjustment for j in self.x], '--')[0].get_color()
            self.colors.append(ColorConverter().to_rgb(color))
        
        if not isGuess:
            fit = genGaussSum(self.x, *popt)
            self.axes.plot(fit + self.totalAdjustment, linewidth=2)
            
        self.ui.gaussPlot.canvas.draw()
    
    def printError(self):
        QtGui.QMessageBox.warning(self, "Error Loading File",
                                  "Please choose a valid Matlab file")
        self.ui.results.append("Invalid input")
    
    def populateSliders(self):
        self.ui.centerSlider.setRange(0, len(self.data))
        self.ui.heightSlider.setRange(0, 1000*max(self.data))
        self.ui.widthSlider.setRange(1, len(self.data)//2)
    
    def getData(self):
        filepath = str(QtGui.QFileDialog.getOpenFileName(self, 'Open File',
                                                         '/u1/lcls/matlab/data'))
        try:
            axdata = loadmat(filepath)
            self.ampList = extract(axdata, 'ampList')
            self.data, self.totalAdjustment, self.step = processData(self.ampList)
            self.x = range(0, len(self.ampList))
            
            self.cleanAndPlotData()
            self.getGuess()
            self.updateOptions(self.ui.numFits, 1, self.potentialPeaks + 1)
            self.updateOptions(self.ui.guessToEdit, 0, 
                               int(self.ui.numFits.currentText()))
            self.changeStates(True)
            self.populateSliders()
            self.updateGuessToEdit()
                
        # Handle case where user cancels file selection
        except IOError:
            if filepath == "":
                pass
            else:
                self.printError()
        
        except ValueError:
            self.printError()
    
    def getGuess(self):
        self.numPeaks = int(self.ui.numFits.currentText())
        self.guess, self.potentialPeaks = getGuess(self.data, self.step, 
                                                   self.useZeros, self.numPeaks)
        
        self.ui.results.append("----------GUESS----------")
        self.plotFit(self.guess, True)
    
    def getFit(self):
        try:
            self.fit = getFit(self.data, self.x, self.guess)
            checkBounds(self.fit, self.data, self.x)
            self.axes.clear()
            self.axes.plot(self.ampList, '.', marker='o')
            self.ui.results.append("----------RESULT----------")
            self.plotFit(self.fit, False)
        except (RuntimeError, AssertionError):
            self.ui.results.append("----------OPTIMIZATION NOT FOUND----------")
            QtGui.QMessageBox.information(self, "Unable to Find Optimized Fit",
                                          "Auto detecting the pedestal or"
                                          + " selecting a different number of "
                                          + "Gaussians to fit might help")
        except AttributeError:
            QtGui.QMessageBox.information(self, "Error Getting Fit",
                                          "Please upload a Matlab XCor data file"
                                          + " to fit")
    
if __name__ == "__main__":
    app = QtGui.QApplication(argv)
    window = GaussianUi()
    window.show()
    exit(app.exec_())
