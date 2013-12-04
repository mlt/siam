# -*- coding: utf-8 -*-

# /***************************************************************************
#  HypsometryDialog
#                                  A QGIS plugin
#  For now only hypsometry around inlets
#                              -------------------
#         begin                : 2013-11-05
#         copyright            : (C) 2013 by Mikhail Titov
#         email                : mlt@gmx.us
#  ***************************************************************************/
#
# /***************************************************************************
#  *                                                                         *
#  *   This program is free software; you can redistribute it and/or modify  *
#  *   it under the terms of the GNU General Public License as published by  *
#  *   the Free Software Foundation; either version 2 of the License, or     *
#  *   (at your option) any later version.                                   *
#  *                                                                         *
#  ***************************************************************************/


from PyQt4 import QtCore, QtGui, uic
from PyQt4.QtCore import QCoreApplication
from qgis.core import *
import logging, logging.config
from os.path import dirname, join
from persistent import PersistentDialog
from hypsometry import Starter
import multiprocessing as mp
import os, sys
import logging
import time

try:
    from osgeo import gdal, osr, ogr
    from osgeo.gdalconst import *
except ImportError:
    import gdal, osr, ogr
    from gdalconst import *

class HypsometryDialog(PersistentDialog):
    def __init__(self):
        PersistentDialog.__init__(self)
        uic.loadUi(join(dirname(__file__), 'ui_hypsometry.ui'), self)
        self.buttonBox.accepted.disconnect()
        self.buttonBox.rejected.disconnect()
        self.buttonRun = self.buttonBox.button(QtGui.QDialogButtonBox.Ok)
        self.buttonAbort = self.buttonBox.button(QtGui.QDialogButtonBox.Abort)
        self.buttonClose = self.buttonBox.button(QtGui.QDialogButtonBox.Close)

        self._log = logging.getLogger('siam.hypsometry')
        # self.handler = LogHandler()
        # self.handler.setLevel(logging.INFO)
        # formatter = logging.Formatter('%(processName)s: %(message)s')
        # self.handler.setFormatter(formatter)

        self._load_settings()

        self.timer = QtCore.QTimer(self)
        self.timer.setInterval(2000)
        # self.timer.setTimerType(QtCore.VeryCoarseTimer) # QT5
        self.timer.timeout.connect(self.on_timer_timeout)
        self.statusBar = QtGui.QStatusBar(self)
        # self.statusBar.setWidth(self.width())
        self.labelTime = QtGui.QLabel(self)
        # self.labelTime.setFrameShape(QtGui.QFrame.NoFrame)
        self.labelTime.setText('Ready')
        self.statusBar.addWidget(self.labelTime)
        self.verticalLayout.addWidget(self.statusBar)
        self.tic = 0
        # self.statusBar.showMessage('ready')

    # def event(self, e):
    #     if e.type() == QtCore.QEvent.StatusTip:
    #         self.statusBar.showMessage(e.tip())
    #         return True
    #     return True

    # @QtCore.pyqtSlot()
    def on_timer_timeout(self):
        self.labelTime.setText('Elapsed: {:d}s'.format(int(time.time()-self.tic)))

    @QtCore.pyqtSlot(str)
    def on_cbConnection_activated(self, name):
        settings = QtCore.QSettings()
        settings.beginGroup('PostgreSQL')
        settings.beginGroup('connections')
        settings.beginGroup(name)
        self.host = settings.value('host')
        self.port = settings.value('port')
        self.dbname = settings.value('database')
        self.user = settings.value('username')
        fname = "PG:dbname='{database:s}' host={host:s} port={port:s} user='{user:s}'".format(
            database = self.dbname,
            host = self.host,
            port = self.port,
            user = self.user
        )
        settings.endGroup()
        settings.endGroup()
        settings.endGroup()

        self.cbInlets.clear()
        ds = ogr.Open(fname)
        if ds is None:
            QtGui.QMessageBox.critical(self, 'Connection failed', 'Failed to connect to a data source "{:s}"'.format(fname))
            return

        for i in range(ds.GetLayerCount()):
            layer = ds.GetLayer(i)
            if ogr.wkbPoint == layer.GetGeomType():
                self.cbInlets.addItem(layer.GetName())

        self.cbDEM.clear()
        layer = ds.ExecuteSQL('select r_table_name from raster_columns order by r_table_name')
        for f in layer:
            self.cbDEM.addItem(f.GetField(0))

        self.cbMap.clear()
        self.cbOutput.clear()
        self.cbPartitions.clear()
        layer = ds.ExecuteSQL("""
SELECT table_schema || '.' || table_name as table
FROM information_schema.tables
--WHERE table_type = 'BASE TABLE'
WHERE table_schema NOT IN ('information_schema', 'pg_catalog')
ORDER BY table_schema,table_name;""")
        for f in layer:
            name = f.GetField(0)
            self.cbMap.addItem(name)
            self.cbOutput.addItem(name)
            self.cbPartitions.addItem(name)

        self.cbMap.setCurrentIndex(-1)
        self.cbMap.setEditText('side_inlets_parts')
        self.cbOutput.setCurrentIndex(-1)
        self.cbOutput.setEditText('hypsometry')
        self.cbPartitions.setCurrentIndex(-1)
        self.cbPartitions.setEditText('dem_parts')

    def _load_settings(self):
        settings = QtCore.QSettings()
        selected = settings.value('PostgreSQL/connections/selected') or settings.value('siam/hypsometry')
        settings.beginGroup('PostgreSQL')
        settings.beginGroup('connections')
        self.cbConnection.addItems( settings.childGroups() )
        settings.endGroup()
        settings.endGroup()
        if selected:
            self.cbConnection.setCurrentIndex( self.cbConnection.findText(selected) )
            self.on_cbConnection_activated( selected )

        settings.beginGroup('siam')
        settings.beginGroup('hypsometry')
        self.loadall(settings)
        settings.endGroup()
        settings.endGroup()

    def _save_settings(self):
        settings = QtCore.QSettings()
        settings.beginGroup('siam')
        settings.beginGroup('hypsometry')
        self.saveall(settings)
        settings.endGroup()
        settings.endGroup()

    # @QtCore.pyqtSlot()
    def on_buttonAbort_clicked(self):
        if QtGui.QMessageBox.Yes == QtGui.QMessageBox.question(
                self,
                QCoreApplication.translate('hypsometry', 'Are you sure?'),
                QCoreApplication.translate('hypsometry', 'abort?'),
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No):
            self.on_finished()
            self.extractor.kill()
            # self.buttonAbort.setEnabled(False)
            # self.tic = 0
            # self.buttonRun.setEnabled(True)

    # @QtCore.pyqtSlot()
    def on_buttonClose_clicked(self):
        if self.tic == 0 or QtGui.QMessageBox.Yes == QtGui.QMessageBox.question(
                self,
                QCoreApplication.translate('hypsometry', 'Are you sure?'),
                QCoreApplication.translate('hypsometry', 'close?'),
                QtGui.QMessageBox.Yes | QtGui.QMessageBox.No):
            self.reject()

    # @QtCore.pyqtSlot(QtGui.QAbstractButton)
    def on_buttonBox_clicked(self, button):
        sb = self.buttonBox.standardButton(button)
        dispatcher = {
            QtGui.QDialogButtonBox.Ok: self.on_buttonRun_clicked,
            QtGui.QDialogButtonBox.Close: self.on_buttonClose_clicked,
            QtGui.QDialogButtonBox.Abort: self.on_buttonAbort_clicked
        }
        if dispatcher.has_key(sb):
            dispatcher[sb]()

    # def reject(self):
    #     QtGui.QMessageBox.question(self, 'hello', 'reject()')

    # def accept(self):
    #     # super(HypsometryDialog, self).accept()

    @QtCore.pyqtSlot()
    def on_buttonRun_clicked(self):
        self._save_settings()
        self.thread = QtCore.QThread(self) # self
        args = dict(
            host=self.host,
            port=self.port,
            dbname=self.dbname,
            user=self.user,
            dem_table=str(self.cbDEM.currentText()),
            dem_parts=str(self.cbPartitions.currentText()),
            parts_map=str(self.cbMap.currentText()),
            layer=str(self.cbInlets.currentText()),
            out="PG:host={host:s} port={port:s} dbname='{dbname:s}' user='{user:s}'".format(
                host=self.host, port=self.port, dbname=self.dbname, user=self.user),
            table=str(self.cbOutput.currentText()),
            max_height=self.sbMaxStage.value(),
            max_area=self.sbMaxArea.value(),
            step=self.sbStep.value(),
            radius=self.sbRadius.value(),
            fixup=False,
            where=None,
            threads=mp.cpu_count(),
            mp=True,
            _loglevel=self._log.getEffectiveLevel()
        )
        self.extractor = Worker(args)
        self.thread.started.connect(self.extractor.start)
        self.thread.finished.connect(self.on_finished)
        # self.extractor.finished.connect(self.thread.quit)
        QtCore.QObject.connect(self.extractor, QtCore.SIGNAL("finished()"),
                               self.thread.quit)
        self.extractor.moveToThread(self.thread)
        self.buttonRun.setEnabled(False)
        self.buttonAbort.setEnabled(True)
        self.thread.start()
        self.tic = time.time()
        self.timer.start()

    def on_finished(self):
        self.timer.stop()
        # self._log.removeHandler(self.handler)
        self.tic = 0
        self.buttonRun.setEnabled(True)
        self.buttonAbort.setEnabled(False)

class Worker(QtCore.QObject, Starter):
    """ Worker object to fetch hypsometry polygons in separate thread """
    finished = QtCore.pyqtSignal()

    def __init__(self, args):
        super(Worker, self).__init__()
        Starter.__init__(self, args)

    @QtCore.pyqtSlot()
    def start(self):
        self.run()
        self.finished.emit()
