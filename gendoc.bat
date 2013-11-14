@ECHO OFF

set OSGEO4W_ROOT=C:\OSGeo4W
set PYTHONPATH=%OSGEO4W_ROOT%\apps\qgis\python;%OSGEO4W_ROOT%\apps\qgis\python\plugins;%HOME%/.qgis/python/plugins;%OSGEO4W_ROOT%\apps\python27\lib\site-packages
set PATH=%OSGEO4W_ROOT%\bin;%OSGEO4W_ROOT%\apps\qgis\bin;%OSGEO4W_ROOT%\apps\qgis\plugins;%OSGEO4W_ROOT%\apps\Python27\Scripts

sphinx-apidoc --force -o help/source .
cd help
call make.bat html
