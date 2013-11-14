# -*- coding: utf-8 -*-
"""
/***************************************************************************
 Siam
                                 A QGIS plugin
 For now only hypsometry around inlets
                              -------------------
        begin                : 2013-11-05
        copyright            : (C) 2013 by Mikhail Titov
        email                : mlt@gmx.us
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
# Import the PyQt and QGIS libraries
from PyQt4.QtCore import *
from PyQt4.QtGui import *
from qgis.core import *
# Initialize Qt resources from file resources.py
import resources_rc
# Import the code for the dialog
from hypsometrydialog import HypsometryDialog
import multiprocessing, logging.config
import os.path, sys


class Siam:

    def __init__(self, iface):
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)
        # initialize locale
        locale = QSettings().value("locale/userLocale")[0:2]
        localePath = os.path.join(self.plugin_dir, 'i18n', 'siam_{}.qm'.format(locale))

        if os.path.exists(localePath):
            self.translator = QTranslator()
            self.translator.load(localePath)

            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)

        if sys.platform == 'win32':
            # OSGeo4W does not bundle python in exec_prefix for python
            path = os.path.abspath(os.path.join(sys.exec_prefix, '../../bin/pythonw.exe'))
            multiprocessing.set_executable(path)
            sys.argv = [ None ]

        # # Create the dialog (after translation) and keep reference
        # self.dlgHypsometry = HypsometryDialog()
        logging.config.fileConfig(
            os.path.join(os.path.dirname(__file__), 'logging.conf'),
            disable_existing_loggers=False)

    def initGui(self):
        # Create action that will start plugin configuration
        self.action = QAction(
            QIcon(":/plugins/siam/icon.png"),
            u"Hypsometry", self.iface.mainWindow())
        # connect the action to the run method
        self.action.triggered.connect(self.runHypsometry)

        # Add toolbar button and menu item
        self.iface.addToolBarIcon(self.action)
        self.iface.addPluginToMenu(u"&Side Inlets Analysis Methods", self.action)

    def unload(self):
        # Remove the plugin menu item and icon
        self.iface.removePluginMenu(u"&Side Inlets Analysis Methods", self.action)
        self.iface.removeToolBarIcon(self.action)

    def runHypsometry(self):
        # show the dialog
        dlgHypsometry = HypsometryDialog()
        dlgHypsometry.show()
        # Run the dialog event loop
        result = dlgHypsometry.exec_()
        # See if OK was pressed
        if result == 1:
            # do something useful (delete the line containing pass and
            # substitute with your code)
            pass
