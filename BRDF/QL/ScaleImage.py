from IPython import embed
import numpy
from scipy.misc import bytescale

def ScaleImage( Image ):
    '''
    Create a byte scaled image with specified histogram
    '''
    #Range = numpy.array([0.0, 0.8])
    Range = numpy.array([0.0, 1.0])
    ScaledImage = bytescale( Image, cmin=Range[0], cmax=Range[1] )

    # MODIS Rapid Response Enhancement for True Colour
    # x array of input values
    x = numpy.array([0,  30,  60, 120, 190, 255])

    # y array of output values
    y = numpy.array([0, 110, 160, 210, 240, 255])

    # Create output array
    rows = Image.shape[0]
    cols = Image.shape[1]

    Scaled = numpy.zeros( (rows,cols), numpy.uint8 )

    for i in range(x.shape[0] - 1):
        x1 = x[i]
        x2 = x[i + 1]
        y1 = y[i]
        y2 = y[i + 1]

        m =  (y2 - y1) / float((x2 - x1))
        b = y2 - (m * x2)

        mask = numpy.where( (ScaledImage >= x1) & (ScaledImage < x2) )
        Scaled[mask] = (m * ScaledImage + b)[mask]

    mask = numpy.where( ScaledImage >= x2 )
    Scaled[mask] = 255

    return Scaled

