import pyfits
import sys


def write_ascii(filename, ext):
    ofile = open(filename.split('.')[0]+'.asc', 'w')
    tbdata = pyfits.getdata(filename, 1)
    wl = tbdata['wavelength'].ravel()
    net = tbdata[ext].ravel()
    ofile.write('Wavelength             %s\n' %(ext))
    for w, n in zip(wl, net):
        ofile.write('%4.8f\t%4.8e\n' %(w, n))
    ofile.close()
                 

if __name__ == "__main__":
    filename = sys.argv[1]
    try:
        ext = sys.argv[2]
    except:
        ext = 'net'
    write_ascii(filename, ext)
