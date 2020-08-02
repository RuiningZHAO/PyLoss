# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'wavecal.ui'
#
# Created by: PyQt5 UI code generator 5.9
#
# WARNING! All changes made in this file will be lost!

import sys

from PyQt5 import QtCore, QtGui, QtWidgets

import numpy as np

from scipy import signal

from astropy.stats import sigma_clip

from matplotlib.backends.qt_compat import QtCore, QtWidgets, is_pyqt5
if is_pyqt5():
    from matplotlib.backends.backend_qt5agg import (
        FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
else:
    from matplotlib.backends.backend_qt4agg import (
        FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure

import matplotlib.pyplot as plt
# Set plot parameters
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

class Ui_MainWindow(object):
        
    def __init__( self, pix, cnt, inverse ):
        '''
        '''
        self.pix = pix * 1
        self.cnt = cnt * 1
        self.inverse = inverse
        
        if self.inverse:
            self.pix = self.pix[::-1]
            self.cnt = self.cnt[::-1]

    def setupUi(self, MainWindow):

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1000, 800)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
        MainWindow.setSizePolicy(sizePolicy)
        MainWindow.setMinimumSize(QtCore.QSize(1000, 800))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName("gridLayout")
        self.GroupBox_peak = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.GroupBox_peak.sizePolicy().hasHeightForWidth())
        self.GroupBox_peak.setSizePolicy(sizePolicy)
        self.GroupBox_peak.setMinimumSize(QtCore.QSize(100, 276))
        self.GroupBox_peak.setObjectName("GroupBox_peak")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.GroupBox_peak)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")

        self.VerticalLayout_arc = QtWidgets.QVBoxLayout()
        self.VerticalLayout_arc.setObjectName("VerticalLayout_arc")

        # 1. Create New Canvas
        self.Figure_arc = plt.figure()
        self.FigureCanvas_arc = FigureCanvas( self.Figure_arc )
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.FigureCanvas_arc.sizePolicy().hasHeightForWidth())
        self.FigureCanvas_arc.setSizePolicy(sizePolicy)
        self.FigureCanvas_arc.setMinimumSize(QtCore.QSize(100, 276))
        self.FigureCanvas_arc.setObjectName("FigureCanvas_arc")
        # 2. Create corresponding Tool Bar
        self.ToolBar_arc = NavigationToolbar( self.FigureCanvas_arc, MainWindow )
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ToolBar_arc.sizePolicy().hasHeightForWidth())
        self.ToolBar_arc.setSizePolicy(sizePolicy)
        self.ToolBar_arc.setMinimumSize(QtCore.QSize(100, 30))
        self.ToolBar_arc.setMaximumSize(QtCore.QSize(16777215, 30))
        self.ToolBar_arc.setObjectName("ToolBar_arc")
        # 3. Add to Layout
        self.VerticalLayout_arc.addWidget( self.ToolBar_arc )
        self.VerticalLayout_arc.addWidget( self.FigureCanvas_arc )

        self.horizontalLayout_4.addLayout(self.VerticalLayout_arc)

        self.VerticalLayout_peak = QtWidgets.QVBoxLayout()
        self.VerticalLayout_peak.setContentsMargins(5, -1, 5, -1)
        self.VerticalLayout_peak.setSpacing(5)
        self.VerticalLayout_peak.setObjectName("VerticalLayout_peak")
        # `Peak seeking parameters` Label
        self.Label_peak = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_peak.sizePolicy().hasHeightForWidth())
        self.Label_peak.setSizePolicy(sizePolicy)
        self.Label_peak.setMinimumSize(QtCore.QSize(180, 25))
        self.Label_peak.setMaximumSize(QtCore.QSize(180, 25))
        self.Label_peak.setObjectName("Label_peak")
        self.VerticalLayout_peak.addWidget(self.Label_peak)
        # `Peak seeking parameters` Line
        self.Line_peak = QtWidgets.QFrame(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Line_peak.sizePolicy().hasHeightForWidth())
        self.Line_peak.setSizePolicy(sizePolicy)
        self.Line_peak.setMinimumSize(QtCore.QSize(180, 0))
        self.Line_peak.setMaximumSize(QtCore.QSize(180, 16777215))
        self.Line_peak.setFrameShape(QtWidgets.QFrame.HLine)
        self.Line_peak.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.Line_peak.setObjectName("Line_peak")
        self.VerticalLayout_peak.addWidget(self.Line_peak)
        # PEAK SEEKING PARAMETERS
        # =======================
        self.FormLayout_peak = QtWidgets.QFormLayout()
        self.FormLayout_peak.setSizeConstraint(QtWidgets.QLayout.SetFixedSize)
        self.FormLayout_peak.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        self.FormLayout_peak.setLabelAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.FormLayout_peak.setFormAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignTop)
        self.FormLayout_peak.setContentsMargins(0, -1, 0, -1)
        self.FormLayout_peak.setHorizontalSpacing(10)
        self.FormLayout_peak.setObjectName("FormLayout_peak")
        # 1.1 Height Label
        self.Label_height = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_height.sizePolicy().hasHeightForWidth())
        self.Label_height.setSizePolicy(sizePolicy)
        self.Label_height.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_height.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_height.setObjectName("Label_height")
        self.FormLayout_peak.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.Label_height)
        # 1.1 Height LineEdit
        self.LineEdit_height = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_height.sizePolicy().hasHeightForWidth())
        self.LineEdit_height.setSizePolicy(sizePolicy)
        self.LineEdit_height.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_height.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_height.setObjectName("LineEdit_height")
        self.FormLayout_peak.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.LineEdit_height)
        
        self.Label_threshold = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_threshold.sizePolicy().hasHeightForWidth())
        self.Label_threshold.setSizePolicy(sizePolicy)
        self.Label_threshold.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_threshold.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_threshold.setObjectName("Label_threshold")
        self.FormLayout_peak.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.Label_threshold)
        self.LineEdit_threshold = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_threshold.sizePolicy().hasHeightForWidth())
        self.LineEdit_threshold.setSizePolicy(sizePolicy)
        self.LineEdit_threshold.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_threshold.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_threshold.setObjectName("LineEdit_threshold")
        self.FormLayout_peak.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.LineEdit_threshold)
        self.Label_distance = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_distance.sizePolicy().hasHeightForWidth())
        self.Label_distance.setSizePolicy(sizePolicy)
        self.Label_distance.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_distance.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_distance.setObjectName("Label_distance")
        self.FormLayout_peak.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.Label_distance)
        self.LineEdit_distance = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_distance.sizePolicy().hasHeightForWidth())
        self.LineEdit_distance.setSizePolicy(sizePolicy)
        self.LineEdit_distance.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_distance.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_distance.setBaseSize(QtCore.QSize(0, 0))
        self.LineEdit_distance.setObjectName("LineEdit_distance")
        self.FormLayout_peak.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.LineEdit_distance)
        self.Label_prominence = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_prominence.sizePolicy().hasHeightForWidth())
        self.Label_prominence.setSizePolicy(sizePolicy)
        self.Label_prominence.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_prominence.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_prominence.setObjectName("Label_prominence")
        self.FormLayout_peak.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.Label_prominence)
        self.LineEdit_prominence = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_prominence.sizePolicy().hasHeightForWidth())
        self.LineEdit_prominence.setSizePolicy(sizePolicy)
        self.LineEdit_prominence.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_prominence.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_prominence.setObjectName("LineEdit_prominence")
        self.FormLayout_peak.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.LineEdit_prominence)
        self.Label_width = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_width.sizePolicy().hasHeightForWidth())
        self.Label_width.setSizePolicy(sizePolicy)
        self.Label_width.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_width.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_width.setObjectName("Label_width")
        self.FormLayout_peak.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.Label_width)
        self.LineEdit_width = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_width.sizePolicy().hasHeightForWidth())
        self.LineEdit_width.setSizePolicy(sizePolicy)
        self.LineEdit_width.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_width.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_width.setObjectName("LineEdit_width")
        self.FormLayout_peak.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.LineEdit_width)
        self.Label_wlen = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_wlen.sizePolicy().hasHeightForWidth())
        self.Label_wlen.setSizePolicy(sizePolicy)
        self.Label_wlen.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_wlen.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_wlen.setObjectName("Label_wlen")
        self.FormLayout_peak.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.Label_wlen)
        self.LineEdit_wlen = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_wlen.sizePolicy().hasHeightForWidth())
        self.LineEdit_wlen.setSizePolicy(sizePolicy)
        self.LineEdit_wlen.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_wlen.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_wlen.setObjectName("LineEdit_wlen")
        self.FormLayout_peak.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.LineEdit_wlen)
        self.Label_rel = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_rel.sizePolicy().hasHeightForWidth())
        self.Label_rel.setSizePolicy(sizePolicy)
        self.Label_rel.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_rel.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_rel.setObjectName("Label_rel")
        self.FormLayout_peak.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.Label_rel)
        self.LineEdit_rel = QtWidgets.QLineEdit(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_rel.sizePolicy().hasHeightForWidth())
        self.LineEdit_rel.setSizePolicy(sizePolicy)
        self.LineEdit_rel.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_rel.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_rel.setObjectName("LineEdit_rel")
        self.FormLayout_peak.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.LineEdit_rel)
        self.Label_plateau = QtWidgets.QLabel(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_plateau.sizePolicy().hasHeightForWidth())
        self.Label_plateau.setSizePolicy(sizePolicy)
        self.Label_plateau.setMinimumSize(QtCore.QSize(97, 20))
        self.Label_plateau.setMaximumSize(QtCore.QSize(97, 20))
        self.Label_plateau.setObjectName("Label_plateau")
        self.FormLayout_peak.setWidget(7, QtWidgets.QFormLayout.LabelRole, self.Label_plateau)
        self.LineEdit_plateau = QtWidgets.QLineEdit(self.GroupBox_peak)
        self.LineEdit_plateau.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_plateau.sizePolicy().hasHeightForWidth())
        self.LineEdit_plateau.setSizePolicy(sizePolicy)
        self.LineEdit_plateau.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_plateau.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_plateau.setObjectName("LineEdit_plateau")
        self.FormLayout_peak.setWidget(7, QtWidgets.QFormLayout.FieldRole, self.LineEdit_plateau)
        self.VerticalLayout_peak.addLayout(self.FormLayout_peak)
        
        # Button 'Seek'
        self.Button_seek = QtWidgets.QPushButton(self.GroupBox_peak)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_seek.sizePolicy().hasHeightForWidth())
        self.Button_seek.setSizePolicy(sizePolicy)
        self.Button_seek.setMinimumSize(QtCore.QSize(180, 20))
        self.Button_seek.setMaximumSize(QtCore.QSize(180, 20))
        self.Button_seek.setObjectName("Button_seek")
        self.VerticalLayout_peak.addWidget(self.Button_seek)

        self.horizontalLayout_4.addLayout(self.VerticalLayout_peak)
        self.horizontalLayout_4.setStretch(0, 1)
        self.gridLayout.addWidget(self.GroupBox_peak, 0, 0, 1, 1)

        self.GroupBox_line = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.GroupBox_line.sizePolicy().hasHeightForWidth())
        self.GroupBox_line.setSizePolicy(sizePolicy)
        self.GroupBox_line.setMinimumSize(QtCore.QSize(220, 100))
        self.GroupBox_line.setMaximumSize(QtCore.QSize(200, 16777215))
        self.GroupBox_line.setObjectName("GroupBox_line")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.GroupBox_line)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")

        self.TableWidget_line = QtWidgets.QTableWidget(self.GroupBox_line)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.TableWidget_line.sizePolicy().hasHeightForWidth())
        self.TableWidget_line.setSizePolicy(sizePolicy)
        self.TableWidget_line.setMinimumSize(QtCore.QSize(200, 100))
        self.TableWidget_line.setMaximumSize(QtCore.QSize(200, 16777215))
        self.TableWidget_line.setObjectName("TableWidget_line")
        self.TableWidget_line.setColumnCount(2)
        self.TableWidget_line.setRowCount(0)

        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        item.setFont(font)
        self.TableWidget_line.setHorizontalHeaderItem(0, item)
        self.TableWidget_line.setColumnWidth( 0, 120 )

        item = QtWidgets.QTableWidgetItem()
        item.setTextAlignment(QtCore.Qt.AlignCenter)
        font = QtGui.QFont()
        font.setPointSize(8)
        font.setBold(True)
        font.setWeight(75)
        item.setFont(font)
        self.TableWidget_line.setHorizontalHeaderItem(1, item)
        self.TableWidget_line.setColumnWidth( 1, 43 )

        self.horizontalLayout_3.addWidget(self.TableWidget_line)

        self.gridLayout.addWidget(self.GroupBox_line, 0, 1, 2, 1)
        self.GroupBox_wave = QtWidgets.QGroupBox(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.GroupBox_wave.sizePolicy().hasHeightForWidth())
        self.GroupBox_wave.setSizePolicy(sizePolicy)
        self.GroupBox_wave.setMinimumSize(QtCore.QSize(100, 276))
        self.GroupBox_wave.setObjectName("GroupBox_wave")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.GroupBox_wave)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.VerticalLayout_fit = QtWidgets.QVBoxLayout()
        self.VerticalLayout_fit.setObjectName("VerticalLayout_fit")

        # PLOT
        # ====
        # 1. Create New Canvas
        self.Figure_fit = plt.figure()
        self.FigureCanvas_fit = FigureCanvas( self.Figure_fit )
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.MinimumExpanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.FigureCanvas_fit.sizePolicy().hasHeightForWidth())
        self.FigureCanvas_fit.setSizePolicy(sizePolicy)
        self.FigureCanvas_fit.setMinimumSize(QtCore.QSize(100, 276))
        self.FigureCanvas_fit.setObjectName("FigureCanvas_fit")
        # 2. Create corresponding Tool Bar
        self.ToolBar_fit = NavigationToolbar( self.FigureCanvas_fit, MainWindow )
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.ToolBar_fit.sizePolicy().hasHeightForWidth())
        self.ToolBar_fit.setSizePolicy(sizePolicy)
        self.ToolBar_fit.setMinimumSize(QtCore.QSize(100, 30))
        self.ToolBar_fit.setMaximumSize(QtCore.QSize(16777215, 30))
        self.ToolBar_fit.setObjectName("ToolBar_fit")
        # 3. Add to Layout
        self.VerticalLayout_fit.addWidget( self.ToolBar_fit )
        self.VerticalLayout_fit.addWidget( self.FigureCanvas_fit )

        self.horizontalLayout.addLayout(self.VerticalLayout_fit)

        self.VerticalLayout_wave = QtWidgets.QVBoxLayout()
        self.VerticalLayout_wave.setContentsMargins(5, -1, 5, -1)
        self.VerticalLayout_wave.setObjectName("VerticalLayout_wave")
        self.Label_wave = QtWidgets.QLabel(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_wave.sizePolicy().hasHeightForWidth())
        self.Label_wave.setSizePolicy(sizePolicy)
        self.Label_wave.setMinimumSize(QtCore.QSize(180, 25))
        self.Label_wave.setMaximumSize(QtCore.QSize(180, 25))
        self.Label_wave.setObjectName("Label_wave")
        self.VerticalLayout_wave.addWidget(self.Label_wave)
        self.Line_wave = QtWidgets.QFrame(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Line_wave.sizePolicy().hasHeightForWidth())
        self.Line_wave.setSizePolicy(sizePolicy)
        self.Line_wave.setMinimumSize(QtCore.QSize(180, 0))
        self.Line_wave.setMaximumSize(QtCore.QSize(180, 16777215))
        self.Line_wave.setFrameShape(QtWidgets.QFrame.HLine)
        self.Line_wave.setFrameShadow(QtWidgets.QFrame.Sunken)
        self.Line_wave.setObjectName("Line_wave")
        self.VerticalLayout_wave.addWidget(self.Line_wave)
        self.FormLayout_wave = QtWidgets.QFormLayout()
        self.FormLayout_wave.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.FormLayout_wave.setContentsMargins(0, -1, 0, -1)
        self.FormLayout_wave.setHorizontalSpacing(17)
        self.FormLayout_wave.setVerticalSpacing(5)
        self.FormLayout_wave.setObjectName("FormLayout_wave")
        self.Label_order = QtWidgets.QLabel(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_order.sizePolicy().hasHeightForWidth())
        self.Label_order.setSizePolicy(sizePolicy)
        self.Label_order.setMinimumSize(QtCore.QSize(90, 20))
        self.Label_order.setMaximumSize(QtCore.QSize(90, 20))
        self.Label_order.setObjectName("Label_order")
        self.FormLayout_wave.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.Label_order)
        self.LineEdit_order = QtWidgets.QLineEdit(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_order.sizePolicy().hasHeightForWidth())
        self.LineEdit_order.setSizePolicy(sizePolicy)
        self.LineEdit_order.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_order.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_order.setObjectName("LineEdit_order")
        self.FormLayout_wave.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.LineEdit_order)
        self.Label_sig_lower = QtWidgets.QLabel(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_sig_lower.sizePolicy().hasHeightForWidth())
        self.Label_sig_lower.setSizePolicy(sizePolicy)
        self.Label_sig_lower.setMinimumSize(QtCore.QSize(90, 20))
        self.Label_sig_lower.setMaximumSize(QtCore.QSize(90, 20))
        self.Label_sig_lower.setObjectName("Label_sig_lower")
        self.FormLayout_wave.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.Label_sig_lower)
        self.LineEdit_sig_lower = QtWidgets.QLineEdit(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_sig_lower.sizePolicy().hasHeightForWidth())
        self.LineEdit_sig_lower.setSizePolicy(sizePolicy)
        self.LineEdit_sig_lower.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_sig_lower.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_sig_lower.setObjectName("LineEdit_sig_lower")
        self.FormLayout_wave.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.LineEdit_sig_lower)
        self.Label_sig_upper = QtWidgets.QLabel(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_sig_upper.sizePolicy().hasHeightForWidth())
        self.Label_sig_upper.setSizePolicy(sizePolicy)
        self.Label_sig_upper.setMinimumSize(QtCore.QSize(90, 20))
        self.Label_sig_upper.setMaximumSize(QtCore.QSize(90, 20))
        self.Label_sig_upper.setObjectName("Label_sig_upper")
        self.FormLayout_wave.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.Label_sig_upper)
        self.LineEdit_sig_upper = QtWidgets.QLineEdit(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_sig_upper.sizePolicy().hasHeightForWidth())
        self.LineEdit_sig_upper.setSizePolicy(sizePolicy)
        self.LineEdit_sig_upper.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_sig_upper.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_sig_upper.setObjectName("LineEdit_sig_upper")
        self.FormLayout_wave.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.LineEdit_sig_upper)
        self.Label_maxiter = QtWidgets.QLabel(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Label_maxiter.sizePolicy().hasHeightForWidth())
        self.Label_maxiter.setSizePolicy(sizePolicy)
        self.Label_maxiter.setMinimumSize(QtCore.QSize(90, 20))
        self.Label_maxiter.setMaximumSize(QtCore.QSize(90, 20))
        self.Label_maxiter.setObjectName("Label_maxiter")
        self.FormLayout_wave.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.Label_maxiter)
        self.LineEdit_maxiter = QtWidgets.QLineEdit(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.LineEdit_maxiter.sizePolicy().hasHeightForWidth())
        self.LineEdit_maxiter.setSizePolicy(sizePolicy)
        self.LineEdit_maxiter.setMinimumSize(QtCore.QSize(73, 20))
        self.LineEdit_maxiter.setMaximumSize(QtCore.QSize(73, 20))
        self.LineEdit_maxiter.setObjectName("LineEdit_maxiter")
        self.FormLayout_wave.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.LineEdit_maxiter)
        self.VerticalLayout_wave.addLayout(self.FormLayout_wave)
        self.Button_fit = QtWidgets.QPushButton(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_fit.sizePolicy().hasHeightForWidth())
        self.Button_fit.setSizePolicy(sizePolicy)
        self.Button_fit.setMinimumSize(QtCore.QSize(180, 20))
        self.Button_fit.setMaximumSize(QtCore.QSize(180, 20))
        self.Button_fit.setObjectName("Button_fit")
        self.VerticalLayout_wave.addWidget(self.Button_fit)
        self.Button_save = QtWidgets.QPushButton(self.GroupBox_wave)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.Button_save.sizePolicy().hasHeightForWidth())
        self.Button_save.setSizePolicy(sizePolicy)
        self.Button_save.setMinimumSize(QtCore.QSize(180, 20))
        self.Button_save.setMaximumSize(QtCore.QSize(180, 20))
        self.Button_save.setObjectName("Button_save")
        self.VerticalLayout_wave.addWidget(self.Button_save)
        self.horizontalLayout.addLayout(self.VerticalLayout_wave)
        self.horizontalLayout.setStretch(0, 1)
        self.gridLayout.addWidget(self.GroupBox_wave, 1, 0, 1, 1)
        self.gridLayout.setRowMinimumHeight(0, 1)
        self.gridLayout.setRowMinimumHeight(1, 1)
        self.gridLayout.setColumnStretch(0, 1)
        MainWindow.setCentralWidget(self.centralwidget)
        self.statusBar = QtWidgets.QStatusBar(MainWindow)
        self.statusBar.setObjectName("statusBar")
        MainWindow.setStatusBar(self.statusBar)

        self.retranslateUi(MainWindow)
        self.initializeUi(MainWindow) # Add by RNZ 08/01/2020
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        MainWindow.setTabOrder(self.LineEdit_height, self.LineEdit_threshold)
        MainWindow.setTabOrder(self.LineEdit_threshold, self.LineEdit_distance)
        MainWindow.setTabOrder(self.LineEdit_distance, self.LineEdit_prominence)
        MainWindow.setTabOrder(self.LineEdit_prominence, self.LineEdit_width)
        MainWindow.setTabOrder(self.LineEdit_width, self.LineEdit_wlen)
        MainWindow.setTabOrder(self.LineEdit_wlen, self.LineEdit_rel)
        MainWindow.setTabOrder(self.LineEdit_rel, self.LineEdit_plateau)
        MainWindow.setTabOrder(self.LineEdit_plateau, self.Button_seek)
        MainWindow.setTabOrder(self.Button_seek, self.LineEdit_order)
        MainWindow.setTabOrder(self.LineEdit_order, self.LineEdit_sig_lower)
        MainWindow.setTabOrder(self.LineEdit_sig_lower, self.LineEdit_sig_upper)
        MainWindow.setTabOrder(self.LineEdit_sig_upper, self.LineEdit_maxiter)
        MainWindow.setTabOrder(self.LineEdit_maxiter, self.Button_fit)
        MainWindow.setTabOrder(self.Button_fit, self.Button_save)
        MainWindow.setTabOrder(self.Button_save, self.TableWidget_line)

    def retranslateUi(self, MainWindow):

        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Wavelength Calibration UI"))
        self.GroupBox_peak.setTitle(_translate("MainWindow", "Peak Seeking"))
        self.Label_peak.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:12pt;\">Peak parameters</span></p></body></html>"))
        self.Label_height.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">height</span></p></body></html>"))
        self.Label_threshold.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">threshold</span></p></body></html>"))
        self.Label_distance.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">distance</span></p></body></html>"))
        self.Label_prominence.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">prominence</span></p></body></html>"))
        self.Label_width.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">width</span></p></body></html>"))
        self.Label_wlen.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">wlen</span></p></body></html>"))
        self.Label_rel.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">rel_height</span></p></body></html>"))
        self.Label_plateau.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">plateau_size</span></p></body></html>"))
        self.Button_seek.setText(_translate("MainWindow", "Seek"))
        self.GroupBox_line.setTitle(_translate("MainWindow", "Line Identification"))
        item = self.TableWidget_line.horizontalHeaderItem(0)
        item.setText(_translate("MainWindow", "WAVELENGTH"))
        item = self.TableWidget_line.horizontalHeaderItem(1)
        item.setText(_translate("MainWindow", "MASK"))
        self.GroupBox_wave.setTitle(_translate("MainWindow", "Wavelength Solution"))
        self.Label_wave.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:12pt;\">Fit parameters</span></p></body></html>"))
        self.Label_order.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">order</span></p></body></html>"))
        self.Label_sig_lower.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">sigma_lower</span></p></body></html>"))
        self.Label_sig_upper.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">sigma_upper</span></p></body></html>"))
        self.Label_maxiter.setText(_translate("MainWindow", "<html><head/><body><p><span style=\" font-size:10pt;\">maxiter</span></p></body></html>"))
        self.Button_fit.setText(_translate("MainWindow", "Fit"))
        self.Button_save.setText(_translate("MainWindow", "Save"))
    
    def initializeUi(self, MainWindow):

        # Initial Plot
        # ------------
        ax = self.Figure_arc.add_subplot( 111 )
        ax.step( self.pix, self.cnt, 'k-', lw = 0.8, where = 'mid' )
        ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 5, width = 1.5, labelsize = 12 )
        ax.tick_params( which = 'minor', direction = 'in', top = True, right = True, length = 3, width = 1.5, labelsize = 12 )
        ax.set_xlim( self.pix.min(), self.pix.max() )
        ax.set_ylim( -0.05, 1.15 )
        ax.set_xlabel( 'Pixels', fontsize = 16 )
        ax.set_ylabel( 'Normalized Counts', fontsize = 16 )
        self.Figure_arc.tight_layout()
        self.FigureCanvas_arc.draw()

        # Seeking Parameters
        # ------------------
        # 1. Default
        self.LineEdit_height.setText( '0' )
        self.LineEdit_threshold.setText( '0' )
        self.LineEdit_distance.setText( '10' )
        self.LineEdit_prominence.setText( '0.05' )
        self.LineEdit_width.setText( '8' )
        self.LineEdit_wlen.setText( '20' )
        self.LineEdit_rel.setText( '0.5' )
        self.LineEdit_plateau.setText( '0' )
        # 2. Validate
        rx_0_1   = QtCore.QRegExp( "^(1|(0(\\.\\d{1,2})?))$" )
        rx_0_100 = QtCore.QRegExp( "^(100|((\\d{1}|[1-9]{1}[0-9]{1})(\\.\\d{1,2})?))$" )
        rx_1_100 = QtCore.QRegExp( "^(100|(([1-9]{1}[0-9]?)(\\.\\d{1,2})?))$" )
        self.LineEdit_height.setValidator(     QtGui.QRegExpValidator( rx_0_1   ) )
        self.LineEdit_threshold.setValidator(  QtGui.QRegExpValidator( rx_0_1   ) )
        self.LineEdit_distance.setValidator(   QtGui.QRegExpValidator( rx_1_100 ) )
        self.LineEdit_prominence.setValidator( QtGui.QRegExpValidator( rx_0_1   ) )
        self.LineEdit_width.setValidator(      QtGui.QRegExpValidator( rx_0_100 ) )
        self.LineEdit_wlen.setValidator(       QtGui.QRegExpValidator( rx_0_100 ) )
        self.LineEdit_rel.setValidator(        QtGui.QRegExpValidator( rx_0_1   ) )
        self.LineEdit_plateau.setValidator(    QtGui.QRegExpValidator( rx_0_100 ) )
        # 3. None empty
        self.LineEdit_height.textChanged.connect(     lambda: self._disable_button_seek( self.LineEdit_height     ) )
        self.LineEdit_threshold.textChanged.connect(  lambda: self._disable_button_seek( self.LineEdit_threshold  ) )
        self.LineEdit_distance.textChanged.connect(   lambda: self._disable_button_seek( self.LineEdit_distance   ) )
        self.LineEdit_prominence.textChanged.connect( lambda: self._disable_button_seek( self.LineEdit_prominence ) )
        self.LineEdit_width.textChanged.connect(      lambda: self._disable_button_seek( self.LineEdit_width      ) )
        self.LineEdit_wlen.textChanged.connect(       lambda: self._disable_button_seek( self.LineEdit_wlen       ) )
        self.LineEdit_rel.textChanged.connect(        lambda: self._disable_button_seek( self.LineEdit_rel        ) )
        self.LineEdit_plateau.textChanged.connect(    lambda: self._disable_button_seek( self.LineEdit_plateau    ) )

        # Fitting Parameters
        # ------------------
        # 1. Default
        self.LineEdit_order.setText( '5' )
        self.LineEdit_sig_lower.setText( '3' )
        self.LineEdit_sig_upper.setText( '3' )
        self.LineEdit_maxiter.setText( '3' )
        # 2. Validate
        rx_1_20 = QtCore.QRegExp( "^(20|([2-9]{1})|(1\\d?))$" )
        rx_0_5  = QtCore.QRegExp( "^(5|([0-4](\\.\\d{1,2})?))$" )
        rx_0_10 = QtCore.QRegExp( "^(10|\\d{1})$" )
        self.LineEdit_order.setValidator(     QtGui.QRegExpValidator(   rx_1_20 ) )
        self.LineEdit_sig_lower.setValidator(  QtGui.QRegExpValidator(  rx_0_5  ) )
        self.LineEdit_sig_upper.setValidator(   QtGui.QRegExpValidator( rx_0_5  ) )
        self.LineEdit_maxiter.setValidator( QtGui.QRegExpValidator(     rx_0_10 ) )
        # 3. None empty
        self.LineEdit_order.textChanged.connect(     lambda: self._disable_button_fit( self.LineEdit_order     ) )
        self.LineEdit_sig_lower.textChanged.connect( lambda: self._disable_button_fit( self.LineEdit_sig_lower ) )
        self.LineEdit_sig_upper.textChanged.connect( lambda: self._disable_button_fit( self.LineEdit_sig_upper ) )
        self.LineEdit_maxiter.textChanged.connect(   lambda: self._disable_button_fit( self.LineEdit_maxiter   ) )
        # 4. Disable
        self.LineEdit_order.setEnabled( False )
        self.LineEdit_sig_lower.setEnabled( False )
        self.LineEdit_sig_upper.setEnabled( False )
        self.LineEdit_maxiter.setEnabled( False )

        # Button
        # ------
        # 1. Seek
        self.Button_seek.clicked.connect( self._find_peaks )
        # 2. Fit
        self.Button_fit.setEnabled( False )
        self.Button_fit.clicked.connect( self._fitting )
        # 3. Save
        self.Button_save.setEnabled( False )
        self.Button_save.clicked.connect( self._save )

    def _find_peaks( self ):
        '''
        '''
        
        # Get peak parameters
        # -------------------
        self.height     = float( self.LineEdit_height.text()     )
        self.threshold  = float( self.LineEdit_threshold.text()  )
        self.distance   = float( self.LineEdit_distance.text()   )
        self.prominence = float( self.LineEdit_prominence.text() )
        self.width      = float( self.LineEdit_width.text()      )
        self.wlen       = float( self.LineEdit_wlen.text()       )
        self.rel        = float( self.LineEdit_rel.text()        )
        self.plateau    = float( self.LineEdit_plateau.text()    )

        self.nrow_old = self.TableWidget_line.rowCount()
        if self.nrow_old == 0:
            # Seek peaks
            # ----------
            self.pkids, _ = signal.find_peaks( x            = self.cnt, 
                                               height       = self.height, 
                                               threshold    = self.threshold, 
                                               distance     = self.distance, 
                                               prominence   = self.prominence, 
                                               width        = self.width, 
                                               wlen         = self.wlen, 
                                               rel_height   = self.rel, 
                                               plateau_size = self.plateau )
            # Write to table
            # --------------
            self.TableWidget_line.setRowCount( len(self.pkids) )
            for i, pkid in enumerate( self.pkids ):
                # Column `WAVELENGTH`
                item = QtWidgets.QTableWidgetItem( '' )
                item.setTextAlignment( QtCore.Qt.AlignCenter )
                font = QtGui.QFont()
                font.setPointSize(8)
                self.TableWidget_line.setItem( i, 0, item )
                # Column `MASK`
                item = QtWidgets.QTableWidgetItem( '0' )
                item.setTextAlignment( QtCore.Qt.AlignCenter )
                font = QtGui.QFont()
                font.setPointSize(8)
                self.TableWidget_line.setItem( i, 1, item )
                # Row height
                self.TableWidget_line.setRowHeight( i, 25 )

        elif self.nrow_old > 0:
            # Get old contents
            # ----------------
            self.pkid = [ pkid for pkid in self.pkids ]
            self.wave = [ self.TableWidget_line.item( i, 0 ).text() for i in range( self.nrow_old ) ]
            self.mask = [ self.TableWidget_line.item( i, 1 ).text() for i in range( self.nrow_old ) ]
            # Seek peaks
            # ----------
            self.pkids, _ = signal.find_peaks( x            = self.cnt, 
                                               height       = self.height, 
                                               threshold    = self.threshold, 
                                               distance     = self.distance, 
                                               prominence   = self.prominence, 
                                               width        = self.width, 
                                               wlen         = self.wlen, 
                                               rel_height   = self.rel, 
                                               plateau_size = self.plateau )
            # Write to table
            # --------------
            self.TableWidget_line.setRowCount( len(self.pkids) )
            for i, pkid in enumerate( self.pkids ):
                if pkid in self.pkid:
                    idx = self.pkid.index( pkid )
                    # Column `WAVELENGTH`
                    item = QtWidgets.QTableWidgetItem( self.wave[idx] )
                    item.setTextAlignment( QtCore.Qt.AlignCenter )
                    font = QtGui.QFont()
                    font.setPointSize(8)
                    self.TableWidget_line.setItem( i, 0, item )
                    # Column `MASK`
                    item = QtWidgets.QTableWidgetItem( self.mask[idx] )
                    item.setTextAlignment( QtCore.Qt.AlignCenter )
                    font = QtGui.QFont()
                    font.setPointSize(8)
                    self.TableWidget_line.setItem( i, 1, item )
                else:
                    # Column `WAVELENGTH`
                    item = QtWidgets.QTableWidgetItem( '' )
                    item.setTextAlignment( QtCore.Qt.AlignCenter )
                    font = QtGui.QFont()
                    font.setPointSize(8)
                    self.TableWidget_line.setItem( i, 0, item )
                    # Column `MASK`
                    item = QtWidgets.QTableWidgetItem( '0' )
                    item.setTextAlignment( QtCore.Qt.AlignCenter )
                    font = QtGui.QFont()
                    font.setPointSize(8)
                    self.TableWidget_line.setItem( i, 1, item )
                # Row height
                self.TableWidget_line.setRowHeight( i, 25 )

        # Plot
        # ----
        self.Figure_arc.clear()
        ax = self.Figure_arc.add_subplot( 111 )
        ax.step( self.pix, self.cnt, 'k-', lw = 0.8, where = 'mid' )
        for i, pkid in enumerate( self.pkids ):
            ax.plot( [self.pix[pkid], self.pix[pkid]], [self.cnt[pkid], 1.05], 'r--', lw = 0.8 )
            ax.plot( [self.pix[pkid], self.pix[pkid]], [          1.12, 1.15], 'r--', lw = 0.8 )
            ax.annotate( F'{i+1:d}', ( self.pix[pkid], 1.08 ), xycoords = 'data', ha = 'center', va = 'center', fontsize = 12, color = 'r' )
        ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 5, width = 1.5, labelsize = 12 )
        ax.tick_params( which = 'minor', direction = 'in', top = True, right = True, length = 3, width = 1.5, labelsize = 12 )
        ax.set_xlim( self.pix.min(), self.pix.max() )
        ax.set_ylim( -0.05, 1.15 )
        ax.set_xlabel( 'Pixels', fontsize = 16 )
        ax.set_ylabel( 'Normalized Counts', fontsize = 16 )
        self.Figure_arc.tight_layout()
        self.FigureCanvas_arc.draw()

        if self.TableWidget_line.rowCount() > 0:
            # Enable Fitting Parameters
            # -------------------------
            self.LineEdit_order.setEnabled( True )
            self.LineEdit_sig_lower.setEnabled( True )
            self.LineEdit_sig_upper.setEnabled( True )
            self.LineEdit_maxiter.setEnabled( True )
            # Enable Button `Fit`
            # -------------------
            self.Button_fit.setEnabled( True )
        else:
            # Enable Fitting Parameters
            # -------------------------
            self.LineEdit_order.setEnabled( False )
            self.LineEdit_sig_lower.setEnabled( False )
            self.LineEdit_sig_upper.setEnabled( False )
            self.LineEdit_maxiter.setEnabled( False )
            # Enable Button `Fit`
            # -------------------
            self.Button_fit.setEnabled( True )       

    def _fitting( self ):
        '''
        '''
        
        self.fit_failed = False

        # Get fit parameters
        # ------------------
        self.order     = int(   self.LineEdit_order.text()     )
        self.sig_lower = float( self.LineEdit_sig_lower.text() )
        self.sig_upper = float( self.LineEdit_sig_upper.text() )
        self.maxiter   = int(   self.LineEdit_maxiter.text()   )
        
        # Get labelled lines
        # ------------------
        wave = np.zeros( self.TableWidget_line.rowCount() )
        mask = np.zeros( self.TableWidget_line.rowCount(), dtype = bool )
        for i in range( self.TableWidget_line.rowCount() ):
            try:
                w = float( self.TableWidget_line.item( i, 0 ).text() )
                m = bool( int( self.TableWidget_line.item( i, 1 ).text() ) )
            except:
                w = np.nan
                m = 1
            wave[i] = w; mask[i] = m

        # Fit
        # ---
        if ( ~mask ).sum() <= self.order:

            self.Button_save.setEnabled( False )

        else:

            for i in range( self.maxiter ):
                z = np.polyfit( self.pix[self.pkids][~mask], wave[~mask], self.order )
                p = np.poly1d( z )
                mask[~mask] = sigma_clip( wave[~mask] - p( self.pix[self.pkids][~mask] ), 
                                          sigma_lower = self.sig_lower,
                                          sigma_upper = self.sig_upper,
                                          maxiters    = 1,
                                          masked      = True ).mask

            z = np.polyfit( self.pix[self.pkids][~mask], wave[~mask], self.order )
            p = np.poly1d( z )
            self.wave_fit = p( self.pix )

            # Plot
            # ----
            self.Figure_fit.clear()
            self.Figure_fit.subplots_adjust( hspace = 0 )
            ax = self.Figure_fit.add_subplot( 211 )
            ax.plot( self.pix[self.pkids][~mask], wave[~mask], 'kx', ms = 6 )
            ax.plot( self.pix[self.pkids][ mask], wave[ mask], 'rx', ms = 6 )
            ax.plot( self.pix, self.wave_fit, 'k-', lw = 0.8 )
            ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 5, width = 1.5, labelsize = 12 )
            ax.tick_params( which = 'minor', direction = 'in', top = True, right = True, length = 3, width = 1.5, labelsize = 12 )
            ax.set_xticklabels([])
            ax.set_xlim( self.pix.min(), self.pix.max() )
            ax.set_ylabel( 'Wavelengths', fontsize = 16 )

            ax = self.Figure_fit.add_subplot( 212 )
            ax.plot( self.pix[self.pkids][~mask], wave[~mask] - p( self.pix[self.pkids][~mask] ), 'kx', ms = 8 )
            ax.plot( self.pix[self.pkids][ mask], wave[ mask] - p( self.pix[self.pkids][ mask] ), 'rx', ms = 8 )
            ax.axhline( y = 0, ls = '--', color = 'k', lw = 0.8 )
            for i, pkid in enumerate( self.pkids ):
                ax.annotate( F'{i+1:d}', ( self.pix[pkid], 0 ), xycoords = 'data', ha = 'center', va = 'center', fontsize = 12, color = 'r' )
            ax.tick_params( which = 'major', direction = 'in', top = True, right = True, length = 5, width = 1.5, labelsize = 12 )
            ax.tick_params( which = 'minor', direction = 'in', top = True, right = True, length = 3, width = 1.5, labelsize = 12 )
            ax.set_xlim( self.pix.min(), self.pix.max() )
            ax.set_xlabel( 'Pixels', fontsize = 16 )
            ax.set_ylabel( 'Residuals', fontsize = 16 )
            self.Figure_fit.align_ylabels()
            self.Figure_fit.tight_layout()
            self.FigureCanvas_fit.draw()

            # Enable Button `Save`
            # --------------------
            self.Button_save.setEnabled( True )

    def _save( self ):
        '''
        '''
        np.savetxt( './wavecal.dat', np.vstack([self.wave_fit, self.cnt]).T, fmt = '%15.8e' )
        self.statusBar.showMessage( 'Saved to `wavecal.dat` successfully.', 2 )

    def _disable_button_seek( self, LineEdit ):
        '''
        '''
        if len( LineEdit.text() ) == 0:
            self.Button_seek.setEnabled( False )
        else:
            self.Button_seek.setEnabled( True )
    
    def _disable_button_fit( self, LineEdit ):
        '''
        '''
        if len( LineEdit.text() ) == 0:
            self.Button_fit.setEnabled( False )
        else:
            self.Button_fit.setEnabled( True )

def main( pix, cnt, inverse ):

    app = QtWidgets.QApplication(sys.argv)

    # Setup UI
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow( pix, cnt, inverse )
    ui.setupUi( MainWindow )
    MainWindow.show()

    sys.exit( app.exec_() )

if __name__ == '__main__':
    
    pix, cnt = np.loadtxt( './wavecal.dat' ).T
    
    inverse = False

    main( pix, cnt, inverse )