import glob
import sys
import os
import convertFOFPropertiesToVTK

if __name__ == '__main__':
    for f in glob.glob(os.path.join(sys.argv[1], '*.fofproperties')):
        prefix, timestamp, _ = os.path.split(f)[-1].split('.')
        timestamp = int(timestamp)
        
        newfile = os.path.join(sys.argv[2], '%s-fofproperties.%04d.vtp' % (prefix, timestamp))
        convertFOFPropertiesToVTK.main(f, newfile)
