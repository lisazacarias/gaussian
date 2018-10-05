#!/usr/local/lcls/package/python/current/bin/python
from __future__ import division

from gaussian_ui import Ui_XCORGaussianFit
from PyQt4 import QtGui

from read_xcor_data import *

import pyqtgraph as pg

from PyQt4.QtGui import QGridLayout

class GaussianUi(QtGui.QMainWindow):
    def __init__(self):
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        QtGui.QMainWindow.__init__(self)
        self.ui = Ui_XCORGaussianFit()
        self.ui.setupUi(self)
        
        self.ui.uploadButton.pressed.connect(self.getData)
        self.ui.numFits.activated.connect(self.updateGuess)
        self.ui.guessToEdit.activated.connect(self.updateGuessToEdit)
        self.ui.filterPedestal.stateChanged.connect(self.setUseZeros)
        self.ui.getFitButton.pressed.connect(self.getFit)
        self.ui.heightSlider.valueChanged.connect(self.changeGaussAmp)
        self.ui.centerSlider.valueChanged.connect(self.changeGaussCenter)
        self.ui.widthSlider.valueChanged.connect(self.changeGaussWidth)

        self.plot = pg.PlotWidget()
        layout = QGridLayout()
        self.ui.gaussPlot.setLayout(layout)
        layout.addWidget(self.plot)
        self.plot.showGrid(True, True)
        
        self.useZeros = not self.ui.filterPedestal.isChecked()
        
        self.changeStates(False)
        self.ui.filterPedestal.setVisible(False)
        
    def changeGaussParam(self, offset, val):
        gauss = int(self.ui.guessToEdit.currentText())
        self.guess[offset + 3*gauss] = val
        self.cleanAndPlotData()
        self.ui.results.append("----------GUESS----------")
        self.plotFit(self.guess, True)
    
    def changeGaussAmp(self):
        self.changeGaussParam(3, int(self.ui.heightSlider.value())/1000)
        
    def changeGaussCenter(self):
        self.changeGaussParam(2, int(self.ui.centerSlider.value()))
        
    def changeGaussWidth(self):
        self.changeGaussParam(4, int(self.ui.widthSlider.value()))
    
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
        self.plot.clear()
        data = pg.ScatterPlotItem(range(0, len(self.ampList)),
                                        self.ampList, symbol='o', size=5)
        self.plot.addItem(data)
    
    def updateOptions(self, comboBox, start, end):
        comboBox.clear()
        for i in xrange(start, end):
            comboBox.addItem(str(i))
        
    def plotFit(self, popt, isGuess):
        self.colors = []
        self.gaussPlots = []
    
        self.ui.results.append("vertical offset: " + str(self.totalAdjustment))
                               
        # Print and plot the optmized line fit.
        # Note that popt has the same format as the guess, meaning that the
        # first two parameters are the m and b of the line, respectively
        self.ui.results.append("line: " + "m = " + str(popt[0]) + ", b = "
                               + str(popt[1]))
        line = pg.PlotCurveItem([popt[0]*j + popt[1] + self.totalAdjustment
                                      for j in self.x],
                                pen=pg.mkPen(style=pg.QtCore.Qt.DotLine, width=3))
        self.plot.addItem(line)
        self.gaussPlots.append(line)

        # Print and plot the optimized gaussian fit(s)
        # Again, the first two elements were the line, and each gaussian is a
        # subsequent group of 3 elements (hence starting at index 2 and
        # incrementing by 3)
        for i in xrange(2, len(popt), 3):
            self.ui.results.append("gaussian " + str(i//3) + ": center = " 
                   + str(popt[i]) + ", amplitude = " + str(popt[i+1]) 
                   + ", width = " + str(popt[i+2]))
            gaussPlot = pg.PlotCurveItem([gaussian(j, popt[i], popt[i + 2],
                                                   popt[i + 1])
                                          + self.totalAdjustment for j in self.x],
                                         pen=pg.mkPen(style=pg.QtCore.Qt.DotLine,
                                                      width=3))
            self.gaussPlots.append(gaussPlot)
            self.plot.addItem(gaussPlot)

        if not isGuess:
            fit = genGaussSum(self.x, *popt)
            self.plot.addItem(pg.PlotCurveItem(fit + self.totalAdjustment,
                                               pen=pg.mkPen(width=3)))

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
            self.ui.numFits.setCurrentIndex(0)
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
            self.plot.clear()
            data = pg.ScatterPlotItem(range(0, len(self.ampList)),
                                      self.ampList, symbol='o', size=5)
            self.plot.addItem(data)

            self.ui.results.append("----------RESULT----------")
            self.plotFit(self.fit, False)
        except (RuntimeError, AssertionError):
            self.ui.results.append("----------OPTIMIZATION NOT FOUND----------")
            QtGui.QMessageBox.information(self, "Unable to Find Optimized Fit",
                                          "Try messing with the guess parameters")
        except AttributeError:
            QtGui.QMessageBox.information(self, "Error Getting Fit",
                                          "Please upload a Matlab XCor data file"
                                          + " to fit")
    
if __name__ == "__main__":
    app = QtGui.QApplication(argv)
    window = GaussianUi()
    window.show()
    exit(app.exec_())
