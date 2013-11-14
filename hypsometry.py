#!/usr/bin/env python

"""Quick & Dirty way to retrieve a hypsometry given a set of starting
points and DEM.

"""

from __future__ import print_function
import logging, logging.config
from os.path import dirname, join
from time import sleep
import numpy as np
import argparse
import multiprocessing as mp
from copy import deepcopy
from MultiProcessingLog import QueueHandler, QueueListener

try:
    from osgeo import gdal, osr, ogr
    from osgeo.gdalconst import *
except ImportError:
    import gdal, osr, ogr
    from gdalconst import *

class Hypsometry:
    "Hypsometry creator"

    def __init__(self, args):
        "Set up DB connection"

        self.__dict__.update(args)
        self._log = logging.getLogger(__name__).getChild(self.__class__.__name__)
        self._log.info("Using GDAL '%s'", gdal.VersionInfo("RELEASE_NAME"))

    def _read(self):
        """Read in input data into memory"""

        connstr = "PG:host={:s} port={:s} dbname={:s} user={:s}".format(self.host, self.port, self.dbname, self.user)
        conn_ogr = ogr.Open(connstr)
        if not self.layer is None:
            self.pts = conn_ogr.GetLayerByName(self.layer)
        else:
            self.pts = conn_ogr.GetLayer()
        self.srs = self.pts.GetSpatialRef()
        if self.srs is None:
            self._log.warn("Failed to get spatial reference for inlets. I hope you do use the same projection!")
        else:
            self.SRID = int(self.srs.GetAttrValue("AUTHORITY", 1))
            self._log.info("Points will be read from '%s' with SRID=%d", connstr, self.SRID)
        if not self.part is None:
            fid = self.pts.GetFIDColumn()
            if fid is None or fid == '':
                self._log.critical("FID column is unknown!! I assume it is 'gid'")
                fid = 'gid'
            where = "{:s} in (select gid from {:s} where pid={:d})".format(fid, self.parts_map, self.part)
            # it is partition table that shall be created with
            # necessary points only if those need to be restricted
            if self.pts.SetAttributeFilter(where):
                self._log.critical("Failed to restrict features for partition #%d", self.part)

        connstr_dem = connstr + " table={:s} mode=2".format(self.dem_table)
        if not self.part is None:
            # we do have to explicitly use table name

            # though shorter it does not leverage indices
            # connstr_dem += " where='exists ( select 1 from dem_parts where ST_Contains(geom, rast::geometry) and gid = {:d})'".format(self.part)
            # version below depends on dem_parts
            connstr_dem += """ where='rid in (
select rid from {:s}, {:s}
where ST_Contains(geom, rast::geometry)
  and pid = {:d})'""".format(self.dem_table, self.dem_parts, self.part)
            # connstr_dem += " where='rid in ( select rid from {:s}, dem_parts where ST_Contains(geom, rast::geometry) and gid = {:d})'".format(self.dem_table, self.part)
#             connstr_dem += """ where='
# rid in (
#   select rid
#   from {:s}, dem_parts
#   where ST_Contains(geom, rast::geometry)
#     and gid = {:d})'""".format('dem50', self.part)
            # But the one below requires range (500).

            # Since temporary tables are to be created from within the
            # code in the future and there will be a distance
            # parameter, it is probably okay

            # From another POV, we can prefer the shorter version
            # above for the same reason!
#             connstr_dem += """ where='
# rid in (
#   select rid
#   from dem50
#   join {:s} si on ST_DWithin(rast::geometry, geom, 500)
#   join {:s} on si.gid = inlet_gid
#   where part_gid = {:d})'""".format(self.layer, self.parts_map, self.part)
        dataset = gdal.Open(connstr_dem, GA_ReadOnly )
        if dataset is None:
            self._log.critical("Can't open %s", connstr_dem)
            raise Exception()
        self.geotransform = dataset.GetGeoTransform()
        self.invgeotransform = gdal.InvGeoTransform(self.geotransform)[1]
        self.proj = dataset.GetProjection()
        srs = osr.SpatialReference(self.proj)
        SRID = int(srs.GetAttrValue("AUTHORITY", 1))
        band = dataset.GetRasterBand(1)
        self.NODATA = band.GetNoDataValue()
        self._log.info("Raster from '%s' is used with SRID=%d and NoDATA=%f", connstr_dem, SRID, self.NODATA)
        if self.SRID != SRID:
            raise Exception("Points & DEM should use same coordinate system")
        self.raster = band.ReadAsArray()

        self.pts.ResetReading()
        self.pts_dict = dict()
        self.z = []
        for pt in self.pts:
            geom = pt.GetGeometryRef()
            self._log.debug("%d => %.3f, %.3f", pt.GetFID(), geom.GetX(), geom.GetY())

            iPixel = int(self.invgeotransform[0] +
                         self.invgeotransform[1] * geom.GetX() +
                         self.invgeotransform[2] * geom.GetY())
            iLine = int(self.invgeotransform[3] +
                        self.invgeotransform[4] * geom.GetX() +
                        self.invgeotransform[5] * geom.GetY() );
            try:
                val = self.raster[iLine, iPixel]
                if val != self.NODATA: # GDAL 1.9 is broken also should not(?) compare floats
                    self._log.debug("%d => %d, %d : %f", pt.GetFID(), iPixel, iLine, val)
                    self.pts_dict[pt.GetFID()] = [pt, val]
                    self.z.append(val)
            except IndexError:
                pass
        self.z.sort()

    def _mkmem(self):
        """Set up in-memory stuff"""

        self._log.debug('Creating in-memory stuff')
        mrd = gdal.GetDriverByName('MEM') #GTiff')#
        self.inmem = mrd.Create('', self.raster.shape[1], self.raster.shape[0], 1, gdal.GDT_Byte)
        if gdal.CE_Failure == self.inmem.SetProjection(self.proj): #dataset.GetProjectionRef()
            self._log.critical('Failed to set projection for in-memory raster')
        self.inmem.SetGeoTransform(self.geotransform)
        self.imb = self.inmem.GetRasterBand(1)
        drv = ogr.GetDriverByName('Memory')
        self.dst_ds = drv.CreateDataSource('out') # is it a must?

    def _getout(self):
        """Connect to shared output layer"""
        self.out_ogr = ogr.Open(self.out, True)
        self.polys = self.out_ogr.GetLayerByName(self.table)
        if self.polys is None:
            self._log.critical('Failed to get an output layer %s', self.table)

    def _process(self, z, z_max):
        """Elevation steps for a 'given' inlet and possibly others"""
        found_any = True
        while z <= z_max and found_any:
            self._log.debug('Filtering for z=%.2f m ...', z)
            tmp = self.raster <= z
            self.imb.WriteArray(tmp.astype(np.byte))
            self._log.debug('Polygonizing...')
            lyr = self.dst_ds.CreateLayer('polys', self.srs, ogr.wkbPolygon)
            fd = ogr.FieldDefn('below', ogr.OFTInteger)
            lyr.CreateField( fd )
            gdal.Polygonize(self.imb, None, lyr, 0)
            if lyr.SetAttributeFilter('below>0'):
                self._log.critical("Failed to restrict polygons")
            found_any = self._findpoly(z, lyr)
            self.dst_ds.DeleteLayer(1)
            z += self.step
        return z

    def _process_queue(self):
        """Main loop to deal with elevations queue"""
        z = -9999
        for z_inlet in self.z:
            if z < z_inlet:
                z = z_inlet
                self._log.debug("Moving to next inlet elevation")
            z_max = z_inlet + self.max_height
            z = self._process(z, z_max)

    def _findpoly(self, z, lyr):
        """Find polygons for inlets"""
        self._log.debug('Searching for one in %d polygons...', lyr.GetFeatureCount())
        found_any = False
        for k, v in self.pts_dict.items():
            if v[1] <= z:
                lyr.ResetReading()
                for p in lyr:
                    poly_geom = p.GetGeometryRef()
                    if v[0].GetGeometryRef().Within(poly_geom):
                        self._log.info('Found polygon for point %d at %.2f', k, z)
                        feat = ogr.Feature(self.polys.GetLayerDefn())
                        feat.SetField('polygon', p.GetFID())
                        feat.SetField('point', k)
                        # Somehow without float() it says NotImplementedError: Wrong number of arguments for overloaded function 'Feature_SetField'.
                        feat.SetField('z', float(z))
                        feat.SetField('stage', float(z - v[1]))
                        feat.SetGeometry(poly_geom)
                        area = poly_geom.GetArea()
                        feat.SetField('area', area)
                        self.polys.CreateFeature(feat)
                        if area > self.max_area or v[1] < z - self.max_height:
                            self._log.debug('Either area (%.1f) is getting too big for %d or it is way below z. Removing.', area, k)
                            del self.pts_dict[k]
                        else:
                            # No need to mark if it was the last for given elevation
                            found_any = True
                        break
        return found_any

    def run(self):
        self._read()
        self._mkmem()
        self._getout()
        self._process_queue()

    # def __del__(self):
    #     "Disconnect from DB"
    #     self._log.debug("Closing ODBC/PG connection")
    #     self.conn.close()

def run_part(args):
    """Entry point for sub-processes.

    Must be plain function on Windows. It can't be un-/bound function:("""
    _log = logging.getLogger(__name__)
    _log.setLevel(args.get('_loglevel', logging.NOTSET))
    hndlr = QueueHandler(args['_queue'])
    _log.addHandler(hndlr)
    _log.info("Working on part %d", args['part'])
    h = Hypsometry(args)
    h.run()
    _log.info("Done working on part %d", args['part'])
    _log.removeHandler(hndlr)

class Starter:
    """Performs initial work before spawning workers."""

    def __init__(self, args):
        self.__dict__.update(args)
        self.args = args        # to pass around lean copy
        self._log = logging.getLogger(__name__).getChild(self.__class__.__name__)
        manager = mp.Manager()
        self.queue = manager.Queue(-1)
        listener = QueueListener(self.queue)

    def mkout(self):
        """Set up output stuff"""

        self._log.info("Making output table '%s'", self.table)
        self.conn_ogr = ogr.Open("PG:host={:s} port={:s} dbname={:s} user={:s}".format(self.host, self.port, self.dbname, self.user))
        if not self.layer is None:
            pts = self.conn_ogr.GetLayerByName(self.layer)
        else:
            pts = self.conn_ogr.GetLayer()
        self.srs = pts.GetSpatialRef()

        out_ogr = ogr.Open(self.out, True)
        polys = out_ogr.CreateLayer(self.table, self.srs, ogr.wkbPolygon, ['OVERWRITE=YES','GEOMETRY_NAME=geom','FID=gid','PG_USE_COPY=YES'])
        if polys is None:
            self._log.critical('Failed to create an output layer %s', self.table)
        fd = ogr.FieldDefn('z', ogr.OFTReal)
        polys.CreateField(fd)
        # polygon id (for a given z) + z make globally unique polygon id
        fd = ogr.FieldDefn('polygon', ogr.OFTInteger)
        polys.CreateField(fd)
        fd = ogr.FieldDefn('point', ogr.OFTInteger)
        polys.CreateField(fd)
        fd = ogr.FieldDefn('stage', ogr.OFTReal)
        polys.CreateField(fd)
        fd = ogr.FieldDefn('area', ogr.OFTReal)
        polys.CreateField(fd)

    def _prepare(self):
        """Set up partitions for parallel processing"""

        self._log.info("Creating partitions in '%s'", self.dem_parts)
        parts = self.conn_ogr.CreateLayer(self.dem_parts, self.srs, ogr.wkbPolygon,
                                          ['OVERWRITE=YES','GEOMETRY_NAME=geom','FID=pid','PG_USE_COPY=YES'])
        if parts is None:
            self._log.critical("Failed to create a layer '%s' for paritions", self.dem_parts)
        # Indices will be implicitly created for this layer
        where = '' if self.where is None else 'and ' + self.where
        lyr = self.conn_ogr.ExecuteSQL("""
    insert into {dem_parts:s}(pid, geom)
    select row_number() over () as pid, (g).geom as geom
    from (
      select ST_Dump(ST_Union(rast::geometry)) as g
      from {side_inlets:s}, {dem_table:s}
        where ST_DWithin(rast::geometry, geom, {radius:f}) {where:s}
        ) foo returning pid;""".format(dem_parts=self.dem_parts, side_inlets=self.layer,
                                       dem_table=self.dem_table, radius=self.radius, where=where))
        cnt = lyr.GetFeatureCount()
        self._log.debug("Having %d partitions", cnt)

        self._log.info("Establishing mapping between inlets and partitions")
        where = '' if self.where is None else 'where ' + self.where
        self.conn_ogr.ExecuteSQL("""
drop table if exists {side_inlets_parts:s};
create table {side_inlets_parts:s} (
  pid int not null,
  gid int primary key
);
insert into {side_inlets_parts:s}(gid, pid)
select gid, pid
from (
  select gid, pid, row_number() over inlet
  from (
    select variance(ST_Value(rast, i.geom)) over part, gid, pid
    from {dem_parts:s} p
        join {side_inlets:s} i on ST_Dwithin(i.geom, p.geom, {radius:f})
    join {dem_table:s} on ST_Intersects(i.geom, rast) {where:s}
    window part as (partition by pid)
  ) foo
  window inlet as (partition by gid order by variance)
) bar
where row_number = 1;
create index on {side_inlets_parts:s}(pid);
        """.format(side_inlets_parts=self.parts_map, side_inlets=self.layer,
                   dem_parts=self.dem_parts, dem_table=self.dem_table,
                   radius=self.radius, where=where))
        return cnt

    def run(self):
        if self.fixup:
            # We have to keep output table in same database to be able
            # to selectively remove records. Alternatively we can
            # fetch list of gid using WHERE and delete those in output
            # table using point field.
            pass
        else:
            self.mkout()
        if self.mp:
            parts = self._prepare()
            self.pool = mp.Pool(processes = min(parts, self.threads))
            # args = vars(args)
            def fix_part(part):
                out = deepcopy(self.args)
                out['part'] = part
                out['_queue'] = self.queue
                return out
            self.pool.map(run_part, [fix_part(x+1) for x in range(parts)])
        else:
            h = Hypsometry(self.args)
            h.run()

    def kill(self):
        """To be used from GUI thread to abort operations"""
        self.pool.terminate()
        self.pool.join()

if __name__ == '__main__':
    logging.config.fileConfig(join(dirname(__file__), 'logging.conf'))
    _log = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--host',
                        default='localhost',
                        help='PGHOST')
    parser.add_argument('--port',
                        default='5432',
                        help='PGPORT')
    parser.add_argument('--dbname',
                        default='gis',
                        help='Database name')
    parser.add_argument('--user',
                        default='user',
                        help='User name to use while connecting to database. Use dot pgpass to supply a password.')
    parser.add_argument('--dem-table',
                        default='dem50',
                        help='Raster table name with DEM')
    parser.add_argument('--dem-parts',
                        default='dem_parts',
                        help='Table to be created with partition polygons')
    parser.add_argument('--parts-map',
                        default='side_inlets_parts',
                        help='Table with inlet-partition mapping')
    parser.add_argument('--layer',
                        default='side_inlets',
                        help='A layer in --points')
    parser.add_argument('--out',
                        default='PG:host=localhost port=5432 dbname=gis user=user',
                        help='An output recognizeable by OGR')
    parser.add_argument('--table', default='hypsometry',
                        help='Table name in DB defined by out. Existing table if any will be dropped!')
    parser.add_argument('--max-height', type=float,
                        default=2,
                        help='Maximum stage')
    parser.add_argument('--max-area', type=float,
                        default=1e5,
                        help='Maximum area')
    parser.add_argument('--step', type=float,
                        default=.01,
                        help='Stage step')
    parser.add_argument('--radius', type=float,
                        default=500,
                        help='Clip/search radius for partitions')
    parser.add_argument('--part', type=int,
                        help='Specific partition to process. Used internally')
    parser.add_argument('--where',
                        help="SQL restriction on side inlets to consider, e.g., 'gid in (1,2)'")
    parser.add_argument('--fixup', action='store_true',
                        help="Don't overwrite output table but delete records for points matching WHERE")
    parser.add_argument('--threads', type=int,
                        default=mp.cpu_count(),
                        help='Maximum number of parallel processes')
    parser.add_argument('--mp', action='store_true',
                        help='Use multiprocess & partitioning')
    args = vars( parser.parse_args() )
    args['_loglevel'] = _log.getEffectiveLevel()
    s = Starter(args)
    s.run()
    _log.info("Exiting")
