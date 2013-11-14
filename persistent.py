# -*- coding: utf-8 -*-
# ***************************************************************************
# PersistentDialog
# ---------------------
# Date : November 2013
# Copyright : (C) 2013 by Mikhail Titov
# Email : mlt at gmx dot us
# ***************************************************************************
# *                                                                         *
# * This program is free software; you can redistribute it and/or modify    *
# * it under the terms of the GNU General Public License as published by    *
# * the Free Software Foundation; either version 2 of the License, or       *
# * (at your option) any later version.                                     *
# *                                                                         *
# ***************************************************************************

"""
Persistent entries dialog synopsis::

    from persistent import PersistentDialog
    from PyQt4.uic import loadUi
    from PyQt4.QtCore import QSettings
    from os.path import dirname, join

    class MyDialog(PersistentDialog):
        def __init__(self):
            PersistentDialog.__init__(self)
            loadUi(join(dirname(__file__), 'ui_mydialog.ui'), self)
            self._load_settings()

        def _load_settings(self):
            ... # extra initialization of, e.g., combo box items goes here
            settings = QSettings()
            settings.beginGroup('myplugin')
            settings.beginGroup('mydialog')
            self.loadall(settings) # magic happens here
            settings.endGroup()
            settings.endGroup()

        def _save_settings(self):
            settings = QSettings()
            settings.beginGroup('myplugin')
            settings.beginGroup('mydialog')
            self.saveall(settings) # here is where magic happens
            settings.endGroup()
            settings.endGroup()

        def accept(self):
            super(MyDialog, self).accept()
            self._save_settings()

.. todo:: Emit signals upon changing entries
"""

__author__ = 'Mikhail Titov'
__date__ = 'November 2013'
__copyright__ = '(C) 2013, Mikhail Titov'

# This will get replaced with a git SHA1 when you do a git archive

__revision__ = '$Format:%H$'


from PyQt4.QtGui import QDialog

class PersistentDialog(QDialog):
    """:qt:`QDialog <qdialog>` with basic user entries persistence using
:qt:`QSettings <qsettings>`

    """

    _setters = { 'QComboBox' : 'setEditText',
                 'QDoubleSpinBox' : 'setValue',
                 'QLineEdit' : 'setText',
                 'QCheckBox' : 'setCheckState'}

    _getters = { 'QComboBox' : 'currentText',
                 'QDoubleSpinBox' : 'value',
                 'QLineEdit' : 'text',
                 'QCheckBox' : 'checkState'}

    _types = { 'QDoubleSpinBox' : float, 'QCheckBox' : int }

    def __cbselect(o, val):
        """Re-select item in combobox if there is such val"""

        index = o.findText(val)
        if index != -1:
            o.setCurrentIndex(index)


    _selectors = { 'QComboBox' : __cbselect }

    def _loadme(self, settings, what):
        o = getattr(self, what)
        c = o.__class__.__name__
        val = settings.value(what, type=self._types.get(c, unicode))
        f = getattr(o, self._setters[c])
        f( val )
        if self._selectors.has_key(c):
            self._selectors[c](o, val)

    def loadall(self, settings):
        """Load dialog entries once descended into necessary settings group"""

        for key in settings.childKeys():
            if self.__dict__.has_key(key):
                o = getattr(self, key)
                self._loadme(settings, key)

    def _saveme(self, settings, what):
        o = getattr(self, what)
        f = getattr(o, self._getters[o.__class__.__name__])
        settings.setValue(what, f())

    def saveall(self, settings):
        """Save dialog entries once descended into necessary settings group"""

        for name in self.__dict__:
            o = getattr(self, name)
            if self._getters.has_key(o.__class__.__name__):
                self._saveme(settings, name)
