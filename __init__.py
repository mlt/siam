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
 This script initializes the plugin, making it known to QGIS.
"""

from qgis.core import QgsMessageLog
import logging

def classFactory(iface):
    # load Siam class from file Siam
    from siam import Siam
    return Siam(iface)

class LogHandler(logging.Handler):
    """Funnel logging into QGIS log console"""
    _levels = { 'DEBUG' : 0,
                'INFO' : 0,
                'WARNING' : 1,
                'ERROR' : 2,
                'CRITICAL' : 2}
    def emit(self, record):
        s = self.format(record)
        QgsMessageLog.logMessage(s, record.name, self._levels[record.levelname])
