'''
Misc functions
'''

import Image, ImageDraw
import numpy as np
try:
    from osgeo import gdal, ogr
except ImportError:
    import gdal, ogr

class Transformer():
    "Transform raster & world coordinates"

    def __init__(self, geotransform):
        self.geotransform = geotransform
        self.invgeotransform = gdal.InvGeoTransform(self.geotransform)[1]

    @staticmethod
    def _transform(transform, x, y):
        xx = transform[0] + transform[1] * x + transform[2] * y
        yy = transform[3] + transform[4] * x + transform[5] * y
        return xx, yy

    def world2pixel(self, x, y):
        return self._transform(self.invgeotransform, x, y)

    def pixel2world(self, x, y):
        return self._transform(self.geotransform, x + .5, y + .5)


def lowest_point_within_polygon(raster, polygon, transformer, NODATA=float('-Inf'), fun=np.argmin):
    """Find the lowest location within a polygon.

    TODO: Need to do something about polygons with holes & multipart polygons.
    Perhaps keep drawing other filled polygons using 1 for holes and using 0 for
    other parts.

    :param numpy.array raster: Raster to analyze
    :param Geometry polygon: The polygon to use as a mask when looking for minimum
    :param Transformer transformer: Coordinate transformation object
    :param float nodata: NODATA value to quickly check on validity
    :return: OGR point
    :rtype: ogr.Geometry
    :raises AssertionError: if stumbled upon NODATA. Supplied polygon is expected to be derived previously from valid data.
    """

    img = Image.new('L', raster.shape[::-1], 1)
    draw = ImageDraw.Draw(img)
    cnt = polygon.GetGeometryCount()
    for i in range(cnt):
        points = polygon.GetGeometryRef(i)
#             if not points.GetGeometryType() in (ogr.wkbLineString,):#ogr.wkbLinearRing,):
        if not (points.GetGeometryType() in (ogr.wkbLineString, ogr.wkbLineString25D)):
            points = points.GetGeometryRef(0)
        rng = xrange(points.GetPointCount())
#         pixels = np.vstack(self._world2pixel(
#                                    np.vectorize(points.GetX, otypes=[np.float])(rng),
#                                    np.vectorize(points.GetY, otypes=[np.float])(rng))).transpose().flatten().tolist() # 67.12 sec
#         pixels = np.vstack(self._world2pixel(
#                                    np.vectorize(points.GetX)(rng),
#                                    np.vectorize(points.GetY)(rng))).transpose().flatten().tolist() # 66.8 sec
#         xx, yy = self._world2pixel(
#                                    np.vectorize(points.GetX)(rng),
#                                    np.vectorize(points.GetY)(rng))
#         pixels = zip(xx, yy) # 67.26 sec
#         pixels = [self._world2pixel(points.GetX(p), points.GetY(p)) for p in range(points.GetPointCount())] # 65.87 sec
        pixels = [transformer.world2pixel(points.GetX(p), points.GetY(p)) for p in rng] # 66.44 sec
        draw.polygon(pixels, fill=0, outline=0)
#         img.save('mask-{:d}.png'.format(self.part))
    mask = np.fromstring(img.tostring(), 'b')
    mask.shape = raster.shape
    masked = np.ma.masked_array(raster, mask=mask)
    y, x = np.unravel_index(fun(masked), mask.shape)
    z = raster[y, x]
    assert not np.isclose(z, NODATA)
    pt = ogr.Geometry(ogr.wkbPoint25D)
    easting, northing = transformer.pixel2world(x, y)
    pt.SetPoint(0, easting, northing, float(z))
    return pt
