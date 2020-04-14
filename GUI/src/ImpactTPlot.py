#!/usr/bin/env python
#This code is to plot the result from ImpactZ
#Input : fort.xx
#Output: figures about beam size and emittance
# plots are saved at '/post'

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import tkinter as tk
from tkinter import ttk,filedialog
import time,os,sys
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.stats import gaussian_kde
import numpy as np


import ParticlePlot, SlicePlot
_height=300
_width =200

IMPACT_T_ADVANCED_PLOT_TYPE= {'Centriod location (mm)'    :2,
                     'Rms size (mm)'             :3,
                     'Centriod momentum (MC)'    :4,
                     'Rms momentum (MC)'         :5,
                     'Twiss'                     :6,
                     'Emittance (mm-mrad)'       :7}

IMPACT_T_SciFormatter = FormatStrFormatter('%2.1E')
IMPACT_T_sciMaxLimit  = 99999 *2
IMPACT_T_sciMinLimit  = 0.0001*2

class AdvancedPlotControlFrame(tk.Toplevel):
    """Output"""

    def __init__(self, master=None, cnf={}, **kw):
        tk.Toplevel.__init__(self, master, cnf, **kw)
        self.title('ImpactT Plot')
        self.focus_set()
        """Plot Control"""
        self.frame_plotButton = tk.Frame(self)
        self.frame_plotButton.grid(column=0, row = 0, pady=5 ,padx=10, sticky="we")

        self.frame_radio = tk.Frame(self.frame_plotButton)
        self.frame_radio.pack(side='top')

        self.plotDirct = tk.IntVar()
        self.plotDirct.set(0)
        self.frame_radio.x = tk.Radiobutton(self.frame_radio, variable=self.plotDirct,
                                           text="X", value=0)
        self.frame_radio.x.pack(side='left')
        self.frame_radio.y = tk.Radiobutton(self.frame_radio, variable=self.plotDirct,
                                           text="Y", value=1)
        self.frame_radio.y.pack(side='left')
        self.frame_radio.z = tk.Radiobutton(self.frame_radio, variable=self.plotDirct,
                                           text="Z", value=2)
        self.frame_radio.z.pack(side='left')

        self.plotTypeComx = tk.StringVar(self.frame_plotButton,'Rms size (mm)')
        self.plotType = ttk.Combobox(self.frame_plotButton,text=self.plotTypeComx,
                                     width = 20,
                                     values=list(IMPACT_T_ADVANCED_PLOT_TYPE.keys()))
        self.plotType.pack(side = 'top')
        self.plot = tk.Button(self.frame_plotButton,text='plot',command=self.makePlot)
        self.plot.pack(fill = 'both',expand =1,side = 'top',padx=10)

        self.t = ttk.Separator(self, orient=tk.HORIZONTAL).grid(column=0, row = 1, sticky="we")


        self.frame2 = tk.Frame(self, height =_height/5, width = _width)
        self.frame2.grid(column=0, row = 2, pady=5 ,padx=10, sticky="nswe")

        rowN=0

        self.button_overall = tk.Button(self.frame2,text='Overall',
                               command = self.overallPlot)
        self.button_overall.grid(row = rowN, column=0,  pady=5 ,padx=5, columnspan = 2, sticky="nswe")
        rowN+=1

        self.button_emitGrowth      = tk.Button(self.frame2,text='EmitGrowth',
                                                command = self.emitGrowthPlot)
        self.button_emitGrowth      .grid(row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_Ek              = tk.Button(self.frame2,text='Kinetic Energy',
                                                command = lambda: self.energyPlot(3,'Kinetic Energy (MeV)'))
        self.button_Ek              .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        '''
        self.button_beta            = tk.Button(self.frame2,text='Beta',
                                                command = lambda: self.energyPlot(4,'Beta'))
        self.button_beta            .grid(row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_gamma           = tk.Button(self.frame2,text='Gamma',
                                                command = lambda: self.energyPlot(2,'Gamma'))
        self.button_gamma           .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1
        '''
        self.button_rmax            = tk.Button(self.frame2,text='Rmax',
                                                command = lambda: self.energyPlot(5,'Rmax (mm)'))
        self.button_rmax            .grid(row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_dw              = tk.Button(self.frame2,text='Rms delta E',
                                                command = lambda: self.energyPlot(6,'Rms delta E (MC^2)'))
        self.button_dw              .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1

        self.button_Temperature         = tk.Button(self.frame2,text='Temperature Plot',
                                                    command = self.makeTemperaturePlot)
        self.button_Temperature         .grid(row = rowN, column=0,  pady=5 ,padx=5, sticky="nswe")
        self.button_Loss                = tk.Button(self.frame2,text='live Particle #',
                                                    command = self.liveParticlePlot)
        self.button_Loss                .grid(row = rowN, column=1,  pady=5 ,padx=5, sticky="nswe")
        rowN+=1

        self.t = ttk.Separator(self.frame2, orient=tk.HORIZONTAL).grid(column=0, row = rowN, columnspan=2,sticky="we")
        rowN+=1

        self.max                        = tk.Button(self.frame2,text='Max amplitude',
                                                    command = self.maxPlot)
        self.max                        .grid(row = rowN, column=0,  pady=5 ,padx=5, columnspan=2,sticky="nswe")
        rowN+=1

        self.button_3order              = tk.Button(self.frame2,text='3 order parameter',
                                                    command = self.make3orderPlot)
        self.button_3order              .grid(row = rowN, column=0,  pady=5 ,padx=5, sticky="nswe")
        self.button_4order              = tk.Button(self.frame2,text='4 order parameter',
                                                    command = self.make4orderPlot)
        self.button_4order              .grid(row = rowN, column=1,  pady=5 ,padx=5, sticky="nswe")
        rowN+=1

        self.t = ttk.Separator(self.frame2, orient=tk.HORIZONTAL).grid(column=0, row = rowN, columnspan=2,sticky="we")
        rowN+=1


        self.button_Particle            = tk.Button(self.frame2,text='Phase Space Plot',
                                                    command = self.ParticlePlot)
        self.button_Particle            .grid(row = rowN, column=0,  pady=5 ,padx=5, sticky="nswe")
        self.button_ParticleDesity1D    = tk.Button(self.frame2,text='Density1D',
                                                    command = self.ParticleDensityPlot1D)
        self.button_ParticleDesity1D    .grid(row = rowN, column=1,  pady=5 ,padx=5, sticky="nswe")
        rowN+=1

        self.button_ParticleDensity     = tk.Button(self.frame2,text='Density2D (by Grid)',
                                                    command = self.ParticleDensityPlot)
        self.button_ParticleDensity     .grid( row = rowN, column=0, pady=5 ,padx=5, sticky="nswe")
        self.button_ParticleDensity2    = tk.Button(self.frame2,text='Density2D (by Ptc)',
                                                    command = self.ParticleDensityPlot2)
        self.button_ParticleDensity2    .grid(row = rowN, column=1, pady=5 ,padx=5, sticky="nswe")
        rowN+=1

        self.t = ttk.Separator(self.frame2, orient=tk.HORIZONTAL).grid(column=0, row = rowN, columnspan=2,sticky="we")
        rowN+=1

        self.button_SlicePlot     = tk.Button(self.frame2,text='Slice plot',
                                                    command = self.SlicePlot)
        self.button_SlicePlot     .grid( row = rowN, column=0, columnspan=2, pady=5 ,padx=5, sticky="nswe")

        rowN+=1



    def overallPlot(self):
        print(self.__class__.__name__)

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Overall plot')

        l=OverallFrame(plotWindow)
        l.pack()

    def energyPlot(self,y,ylabel):
        print(sys._getframe().f_back.f_code.co_name)

        plotWindow = tk.Toplevel(self)
        plotWindow.title(sys._getframe().f_back.f_code.co_name)

        # Column 6, rms energy spread, is not available in per-bunch output
        if y == 6:
             per_bunch=False
        else:
             per_bunch=True

        l=PlotFrame(plotWindow, 'fort.18', y, ylabel, per_bunch)
        l.pack()

    def emitGrowthPlot(self):
        print(sys._getframe().f_back.f_code.co_name)

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')

        l=EmitGrowthFrame(plotWindow)
        l.pack()

    def makeTemperaturePlot(self):
        print((self.plotType))

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')

        l=TemperatureFrame(plotWindow)
        l.pack()

    def liveParticlePlot(self):
        print(sys._getframe().f_back.f_code.co_name)

        plotWindow = tk.Toplevel(self)
        plotWindow.title(sys._getframe().f_back.f_code.co_name)

        l=PlotFrame(plotWindow, 'fort.28', 4, 'Live particle number')
        l.pack()

    def ParticlePlot(self):
        print(self.__class__.__name__)
        filename = filedialog.askopenfilename(parent=self)
        try:
            t=open(filename)
            t.close()
        except:
            return

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Phase Space Plot')

        l=ParticlePlot.ParticleFrame(plotWindow,filename,1.0,'ImpactT')
        l.pack()

    def ParticleDensityPlot(self):
        print(self.__class__.__name__)
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')

        l=ParticlePlot.ParticleDensityFrame_weight2D(plotWindow,fileName,1.0,'ImpactT')
        l.pack()

    def ParticleDensityPlot1D(self):
        print(self.__class__.__name__)
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')

        l=ParticlePlot.ParticleDensityFrame_weight1D(plotWindow,fileName,1.0,'ImpactT')
        l.pack()

    def ParticleDensityPlot2(self):
        print(self.__class__.__name__)
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return
        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')

        l=ParticlePlot.ParticleDensityFrame2D_slow(plotWindow,fileName,1.0,'ImpactT')
        l.pack()

    def SlicePlot(self):
        fileName=filedialog.askopenfilename(parent=self)
        try:
            t=open(fileName)
            t.close()
        except:
            return

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Slice Plot')

        l=SlicePlot.SliceBaseFrame(plotWindow,fileName)
        l.pack()

    def makePlot(self):
        print(self.__class__.__name__)

        PlotFileName='fort.'+str(self.plotDirct.get()+24)
        yx=IMPACT_T_ADVANCED_PLOT_TYPE[self.plotType.get()]
        yl=yx if self.plotDirct.get()!=2 else yx-1

        plotWindow = tk.Toplevel(self)
        plotWindow.title('Plot')

        l=PlotFrame(plotWindow, PlotFileName, yl, self.plotType.get())
        l.pack()


    def maxPlot(self):
        print(self.__class__.__name__)
        filename = 'fort.27'
        try:
            t=open(filename)
            t.close()
        except:
            return

        plotWindow = tk.Toplevel(self)
        plotWindow.title('maxPlot')

        l=PlotMaxFrame(plotWindow,filename)
        l.pack()
    def make3orderPlot(self):
        print(self.__class__.__name__)
        filename = 'fort.29'
        try:
            t=open(filename)
            t.close()
        except:
            return

        plotWindow = tk.Toplevel(self)
        plotWindow.title('make3orderPlot')

        l=Plot3orderFrame(plotWindow,filename)
        l.pack()

    def make4orderPlot(self):
        print(self.__class__.__name__)
        filename = 'fort.30'
        try:
            t=open(filename)
            t.close()
        except:
            return

        plotWindow = tk.Toplevel(self)
        plotWindow.title('make4orderPlot')

        l=Plot4orderFrame(plotWindow,filename)
        l.pack()

class PlotBaseFrame(tk.Frame):
    """Basic plot object that other objects inherit from"""
    def __init__(self, parent, per_bunch=True):
        tk.Frame.__init__(self, parent)
        self.Nbunch = int(parent.master.master.Nbunch.get())
        self.create_option_frame(self.Nbunch, per_bunch)
        self.create_figure()
        self.create_canvas()
        self.create_toolbar()
    def create_option_frame(self, Nbunch, per_bunch=True):
        """Creates a frame to select options of the plot"""
        self.option_frame = tk.Frame(self)
        self.option_frame.pack()
        self.create_bunch_selector(Nbunch, per_bunch)
        self.create_xaxis_selector()
        self.create_plot_button()
    def create_bunch_selector(self, Nbunch, per_bunch=True):
        """Creates bunch selector as part of the option frame"""
        if per_bunch and Nbunch > 1:
            self.bunch_list = ['All']
            self.bunch_list.extend(range(1, Nbunch + 1))
            self.bunch_default = tk.StringVar(self.option_frame, 'All')
            self.bunch_label = tk.Label(self.option_frame,
                                        text="Select bunch: ")
            self.bunch_label.pack(side='left')
            self.bunch_select = ttk.Combobox(self.option_frame,
                                             text=self.bunch_default,
                                             width=6,
                                             values=self.bunch_list)
            self.bunch_select.pack(fill='both', expand=1, side='left')
    def create_xaxis_selector(self):
        """Creates x-axis selector as part of the option frame"""
        self.xaxis_list = ['z', 't']
        self.xaxis_default = tk.StringVar(self.option_frame, 'z')
        self.xaxis_label = tk.Label(self.option_frame, text="Select x-axis: ")
        self.xaxis_label.pack(side='left')
        self.xaxis_select = ttk.Combobox(self.option_frame,
                                         text=self.xaxis_default,
                                         width=6,
                                         values=self.xaxis_list)
        self.xaxis_select.pack(fill='both', expand=1, side='left')
    def create_plot_button(self):
        """Creates plot button as part of the option frame"""
        self.plot_button = tk.Button(self.option_frame,
                                     text="Plot",
                                     foreground="blue",
                                     bg="red",
                                     font=("Verdana", 12),
                                     command=self.plot)
        self.plot_button.pack(fill='both', expand=1, side='left')
    def create_figure(self):
        """Creates the main figure section which will later hold the plot"""
        self.fig = Figure(figsize=(7,5), dpi=100)
        self.subfig = self.fig.add_subplot(111)
    def create_canvas(self):
        """Create and draw the canvas object"""
        self.canvas = FigureCanvasTkAgg(self.fig, self)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side=tk.BOTTOM,
                                         fill=tk.BOTH,
                                         expand=True)
    def create_toolbar(self):
        """Create the bottom toolbar"""
        self.toolbar = NavigationToolbar2Tk(self.canvas, self)
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
    def get_selected_bunch(self):
        """Find which bunch is selected in the option frame"""
        if hasattr(self, "bunch_select"):
            return self.bunch_select.get()
        else:
            return 'All'
    def get_selected_xaxis(self):
        """Find which x axis is selected in the option frame"""
        if hasattr(self, "xaxis_select"):
            return self.xaxis_select.get()
        else:
            return 'z'
    def get_filelist(self):
        selected_bunch = self.get_selected_bunch()
        default_filelist = self.get_default_filelist()
        if selected_bunch == 'All':
            return default_filelist
        else:
            return self.get_bunch_filelist(default_filelist,
                                           int(selected_bunch))
    def get_bunch_filelist(self, filelist, bunch):
        for idx, filename in enumerate(filelist):
            basename = '.'.join(filename.split('.')[0:-1])
            extension = filename.split('.')[-1]
            new_extension = str(int(extension) + bunch*1000)
            new_filename = '.'.join([basename, new_extension])
            filelist[idx] = new_filename
        return filelist
    def get_default_filelist(self):
        raise NotImplementedError()
    def get_xaxis_parameters(self, tcol, zcol):
        selected_xaxis = self.get_selected_xaxis()
        if selected_xaxis == 't':
            return [tcol, 'time (s)']
        else:
            return [zcol, 'z direction (m)']

class PlotFrame(PlotBaseFrame):
    def __init__(self, parent, PlotFileName, yl, labelY, per_bunch=True):
        PlotBaseFrame.__init__(self, parent, per_bunch)
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.45,
                                  box.y0*1.1,
                                  box.width,
                                  box.height])
        self.plot_file = int(PlotFileName.split('.')[-1])
        self.ycol, self.ylabel = yl, labelY
        self.plot()
    def plot(self):
        # Get selected options
        PlotFileName = self.get_filelist()[0]
        self.xcol, self.xlabel = self.get_xaxis_parameters()
        # Open file
        try:
            fin = open(PlotFileName,'r')
        except:
            print(( "  ERROR! Can't open file '" + PlotFileName + "'"))
            return
        # Read in data
        linesList  = fin.readlines()
        fin.close()
        linesList  = [line.split() for line in linesList ]
        x   = np.array([float(xrt[self.xcol]) for xrt in linesList])
        y   = np.array([float(xrt[self.ycol]) for xrt in linesList])
        # Convert units
        if self.ylabel in ['Centriod location (mm)','Rms size (mm)','Rmax (mm)']:
            y = y*1.0e3       # unit convert from m to mm
        elif self.ylabel in ['Emittance (mm-mrad)']:
            y = y*1.0e6       # unit convert from (m-rad) to (mm-mrad)
        # Plot data
        self.subfig.clear()
        self.subfig.plot(x,y)
        self.subfig.set_xlabel(self.xlabel)
        self.subfig.set_ylabel(self.ylabel)
        # Set axes scales
        xMax = np.max(x)
        xMin = np.min(x)
        yMax = np.max(y)
        yMin = np.min(y)
        if (xMax-xMin)>IMPACT_T_sciMaxLimit or (xMax-xMin)<IMPACT_T_sciMinLimit:
            self.subfig.xaxis.set_major_formatter(IMPACT_T_SciFormatter)
        if (yMax-yMin)>IMPACT_T_sciMaxLimit or (yMax-yMin)<IMPACT_T_sciMinLimit:
            self.subfig.yaxis.set_major_formatter(IMPACT_T_SciFormatter)
        # Redraw
        self.canvas.draw()
    def get_default_filelist(self):
        return [f'fort.{self.plot_file}']
    def get_xaxis_parameters(self):
        return PlotBaseFrame.get_xaxis_parameters(self, tcol=0, zcol=1)
    def quit(self):
        self.destroy()

class OverallFrame(PlotBaseFrame):
    def __init__(self, parent):
        PlotBaseFrame.__init__(self, parent)
        self.plot()
    def create_figure(self):
        """Overrides the base method with four subplots"""
        self.fig = Figure(figsize=(12,5), dpi=100)
        self.subfig = []
        self.subfig.append(self.fig.add_subplot(221))
        self.subfig.append(self.fig.add_subplot(222))
        self.subfig.append(self.fig.add_subplot(223))
        self.subfig.append(self.fig.add_subplot(224))
        for i in range(0, 4):
            box = self.subfig[i].get_position()
            self.subfig[i].set_position([box.x0*1.1,
                                         box.y0*1.1,
                                         box.width,
                                         box.height*0.88])
    def plot(self):
        filelist = self.get_filelist()
        xcol, xlabel = self.get_xaxis_parameters()

        picNum = 4
        fileList    = [[]*2]*picNum
        saveName    = []
        labelList   = [[]*2]*picNum
        xdataList   = [[]*2]*picNum
        ydataList   = [[]*2]*picNum
        xyLabelList = [[]*2]*picNum

        saveName.append('sizeX')
        fileList[0]     = filelist[0]
        labelList[0]    = ['rms.X','max.X']
        xdataList[0]    = [xcol, xcol]
        ydataList[0]    = [4,3]
        xyLabelList[0]  = [xlabel, 'beam size in X (mm)']

        saveName.append('sizeY')
        fileList[1]     = filelist[1]
        labelList[1]    = ['rms.Y','max.Y']
        xdataList[1]    = [xcol, xcol]
        ydataList[1]    = [4,5]
        xyLabelList[1]  = [xlabel, 'beam size in Y (mm)']

        saveName.append('sizeZ')
        fileList[2]     = filelist[2]
        labelList[2]    = ['rms.Z','max.Z']
        xdataList[2]    = [xcol, xcol]
        ydataList[2]    = [3,7]
        xyLabelList[2]  = [xlabel, 'beam size in Z (mm)']

        saveName.append('emitXY')
        fileList[3]     = filelist[3]
        labelList[3]    = ['emit.nor.X','emit.nor.Y']
        xdataList[3]    = [xcol, xcol]
        ydataList[3]    = [8,8]
        xyLabelList[3]  = [xlabel, 'emittance at X and Y (mm*mrad)']

        lineType = ['r-','b--']

        for i in range(0,picNum):
            self.subfig[i].clear()
            for j in range(0,2):
                try:
                    fin = open(fileList[i][j],'r')
                except:
                    print("ERRPR Can't open file ' " + fileList[i][j] + "'")
                    return
                linesList  = fin.readlines()
                fin .close()
                linesList  = [line.split() for line in linesList ]
                xId = xdataList[i][j]-1
                yId = ydataList[i][j]-1
                x   = np.array([float(xrt[xId]) for xrt in linesList])
                y   = np.array([float(xrt[yId]) for xrt in linesList])
                if i in range(0,picNum-1):
                    y=y*1.0e3
                elif i == picNum-1:
                    y=y*1.0e6
                self.subfig[i].plot(x, y, lineType[j], linewidth=2, label=labelList[i][j])


            self.subfig[i].set_xlabel(xyLabelList[i][0])
            self.subfig[i].set_ylabel(xyLabelList[i][1])

            xMax = np.max(x)
            xMin = np.min(x)
            yMax = np.max(y)
            yMin = np.min(y)
            if (xMax-xMin)>IMPACT_T_sciMaxLimit or (xMax-xMin)<IMPACT_T_sciMinLimit:
                self.subfig[i].xaxis.set_major_formatter(IMPACT_T_SciFormatter)
            if (yMax-yMin)>IMPACT_T_sciMaxLimit or (yMax-yMin)<IMPACT_T_sciMinLimit:
                self.subfig[i].yaxis.set_major_formatter(IMPACT_T_SciFormatter)

            self.subfig[i].legend(loc='upper center', bbox_to_anchor=(0.5, 1.21),fancybox=True, shadow=True, ncol=5)
        self.canvas.draw()
    def get_default_filelist(self):
        filelist = [[]*2]*4
        filelist[0] = ['fort.24', 'fort.27']
        filelist[1] = ['fort.25', 'fort.27']
        filelist[2] = ['fort.26', 'fort.27']
        filelist[3] = ['fort.24', 'fort.25']
        return filelist
    def get_bunch_filelist(self, filelist, bunch):
        """Extends the base method for the four sublists"""
        for idx, sublist in enumerate(filelist):
            new_sublist = PlotBaseFrame.get_bunch_filelist(self, sublist, bunch)
            filelist[idx] = new_sublist
        return filelist
    def get_xaxis_parameters(self):
        return PlotBaseFrame.get_xaxis_parameters(self, tcol=1, zcol=2)

class EmitGrowthFrame(PlotBaseFrame):
    def __init__(self, parent):
        PlotBaseFrame.__init__(self, parent)
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.4, box.y0, box.width, box.height])
        self.plot()
    def plot(self):
        fileList = self.get_filelist()
        xcol, xlabel = self.get_xaxis_parameters()
        xdataList = [xcol, xcol]
        ydataList = [8, 8]
        xyLabelList = [xlabel, 'Avg emit growth in X and Y']
        lineType = ['r-','b--']
        try:
            fin1 = open(fileList[0],'r')
        except:
            print("  ERROR! Can't open file '" + fileList[0] + "'")
            return
        try:
            fin2 = open(fileList[1],'r')
        except:
            print("  ERROR! Can't open file '" + fileList[1] + "'")
            return
        linesList1  = fin1.readlines()
        linesList2  = fin2.readlines()
        fin1 .close()
        fin2 .close()
        linesList1  = [line.split() for line in linesList1 ]
        linesList2  = [line.split() for line in linesList2 ]
        xId = xdataList[0]-1
        yId = ydataList[0]-1
        try:
            x   = [float(xrt[xId]) for xrt in linesList1]
            start = (float(linesList1[0][yId]) + float(linesList2[0][yId]))/2
            if start < 1.0e-16:
                start=1.0e-16
            y   = [(float(linesList1[k][yId]) + float(linesList2[k][yId]))/2 / start -1 for k in range(len(linesList1))]
        except:
            print("  ERROR! Can't read data '" + fileList[1] + "'")
        self.subfig.cla()
        self.subfig.plot(x, y, lineType[0], linewidth=2, label='emit.growth')
        self.subfig.set_xlabel(xyLabelList[0])
        self.subfig.set_ylabel(xyLabelList[1])
        self.subfig.legend()
        self.canvas.draw()
    def get_default_filelist(self):
        return ['fort.24', 'fort.25']
    def get_xaxis_parameters(self):
        return PlotBaseFrame.get_xaxis_parameters(self, tcol=1, zcol=2)

class TemperatureFrame(PlotBaseFrame):
    def __init__(self, parent):
        PlotBaseFrame.__init__(self, parent)
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.2, box.y0, box.width, box.height])
        self.plot()
    def plot(self):
        filelist = self.get_filelist()
        xcol, xlabel = self.get_xaxis_parameters()
        labelList= ['X','Y','Z']
        lineType = ['-','--',':']
        col      = ['b','g','r']
        linew    = [2,2,3]
        picNum = len(filelist)
        plotPath = './post'
        if os.path.exists(plotPath) == False:
            os.makedirs(plotPath)
        self.subfig.cla()
        for i in range(0, picNum):
            try:
                fin = open(filelist[i],'r')
            except:
                print( "  ERROR! Can't open file '" + filelist[i] + "'")
                return
            linesList  = fin.readlines()
            fin.close()
            linesList  = [line.split() for line in linesList ]
            x = [float(xrt[xcol]) for xrt in linesList]
            if i==2:
                ycol = 4
            else:
                ycol = 5
            y = [float(xrt[ycol])*float(xrt[ycol]) for xrt in linesList]
            self.subfig.plot(x, y,
                             color=col[i],
                             linestyle=lineType[i],
                             linewidth=linew[i],
                             label=labelList[i])
        self.subfig.set_xlabel(xlabel)
        self.subfig.set_ylabel('Temperature')
        self.subfig.legend()
        self.canvas.draw()
    def get_default_filelist(self):
        return ['fort.24', 'fort.25', 'fort.26']
    def get_xaxis_parameters(self):
        return PlotBaseFrame.get_xaxis_parameters(self, tcol=0, zcol=1)

class PlotHighOrderBaseFrame(PlotBaseFrame):
    particle_directions = {'X (mm)'    :2,
                           'Px (MC)'   :3,
                           'Y (mm)'    :4,
                           'Py (MC)'   :5,
                           'Z (mm)'    :6,
                           'Pz (MC)'   :7}
    data = np.array([])
    def __init__(self, parent, PlotFileName):
        PlotBaseFrame.__init__(self, parent)
        box = self.subfig.get_position()
        self.subfig.set_position([box.x0*1.4, box.y0, box.width, box.height])
        self.plot_file = int(PlotFileName.split('.')[-1])
        self.current_file = ''
        self.plot()
    def create_option_frame(self, Nbunch, per_bunch=True):
        """Overrides base class, adding a y-axis selector"""
        self.option_frame = tk.Frame(self)
        self.option_frame.pack()
        self.create_bunch_selector(Nbunch, per_bunch)
        self.create_xaxis_selector()
        self.create_yaxis_selector()
        self.create_plot_button()
    def create_yaxis_selector(self):
        self.yaxis_list = [direction for direction in self.particle_directions]
        self.yaxis_default = tk.StringVar(self.option_frame, 'X (mm)')
        self.yaxis_label = tk.Label(self.option_frame,
                                    text="Select y-axis dimension: ")
        self.yaxis_label.pack(side='left')
        self.yaxis_select = ttk.Combobox(self.option_frame,
                                         text=self.yaxis_default,
                                         width=6,
                                         values=self.yaxis_list)
        self.yaxis_select.pack(fill='both', expand=1, side='left')
    def load(self, PlotFileName):
        """Load data from file into plot object data array"""
        try:
            self.data = np.loadtxt(PlotFileName)
        except:
            print(( "  ERROR! Can't open file '" + PlotFileName + "'"))
            return
        self.data = np.transpose(self.data)
        for i in range(0,6,2):
            self.data[i] = self.data[i] * 1e3  # from m to mm
    def plot(self):
        # Get selected options
        PlotFileName = self.get_filelist()[0]
        self.xcol, self.xlabel = self.get_xaxis_parameters()
        self.ycol, self.ylabel = self.get_yaxis_parameters()
        # Load data
        if PlotFileName != self.current_file:
            print(f'Loading data from {PlotFileName}')
            self.load(PlotFileName)
            self.current_file = PlotFileName
        # Plot data
        self.subfig.cla()
        self.subfig.plot(self.data[self.xcol], self.data[self.ycol])
        self.format_axes(self.xcol, self.ycol)
        self.subfig.set_xlabel(self.xlabel)
        self.subfig.set_ylabel(self.ylabel)
        self.canvas.draw()
    def get_default_filelist(self):
        return [f'fort.{self.plot_file}']
    def get_xaxis_parameters(self):
        return PlotBaseFrame.get_xaxis_parameters(self, tcol=0, zcol=1)
    def get_yaxis_parameters(self):
        selected_yaxis = self.yaxis_select.get()
        return [self.particle_directions[selected_yaxis], selected_yaxis]
    def format_axes(self, xcol, ycol):
        xMax = np.max(self.data[xcol])
        xMin = np.min(self.data[xcol])
        yMax = np.max(self.data[ycol])
        yMin = np.min(self.data[ycol])
        if (   (xMax - xMin) > IMPACT_T_sciMaxLimit
            or (xMax - xMin) < IMPACT_T_sciMinLimit):
            self.subfig.xaxis.set_major_formatter(IMPACT_T_SciFormatter)
        if (   (yMax - yMin) > IMPACT_T_sciMaxLimit
            or (yMax - yMin) < IMPACT_T_sciMinLimit):
            self.subfig.yaxis.set_major_formatter(IMPACT_T_SciFormatter)

class PlotMaxFrame(PlotHighOrderBaseFrame):
    def __init__(self, parent, PlotFileName):
        PlotHighOrderBaseFrame.__init__(self, parent, PlotFileName)
    def get_yaxis_parameters(self):
        ycol, ylabel = PlotHighOrderBaseFrame.get_yaxis_parameters(self)
        return [ycol, 'Max ' + ylabel]

class Plot3orderFrame(PlotHighOrderBaseFrame):
    def __init__(self, parent, PlotFileName):
        PlotHighOrderBaseFrame.__init__(self, parent, PlotFileName)
    def get_yaxis_parameters(self):
        ycol, ylabel = PlotHighOrderBaseFrame.get_yaxis_parameters(self)
        return [ycol, 'Cubic root of 3rd moment of ' + ylabel]

class Plot4orderFrame(PlotHighOrderBaseFrame):
    def __init__(self, parent, PlotFileName):
        PlotHighOrderBaseFrame.__init__(self, parent, PlotFileName)
    def get_yaxis_parameters(self):
        ycol, ylabel = PlotHighOrderBaseFrame.get_yaxis_parameters(self)
        return [ycol, 'Square square root of 4th moment of ' + ylabel]
