from optparse import OptionParser
from PlotAnalysisLib import *

# Notes: There could be issues if the times to be analyzed overlap themselves

# command line options
parser = OptionParser()
parser.add_option("-f", "--file", dest = "pickledfile",
                  help = "Path to pickled general information file",
                  metavar = "FILE")
parser.add_option("-b", "--bands", dest = "pickledbands",
                  help = "path to pickled list of anlaysis and secondary \
channels", metavar = "FILE")
parser.add_option("-t", "--times", dest = "pickledtimes",
                  help = "path to pickled list of times used for analysis",
                  metavar = "FILE")
parser.add_option("-d", "--dir", dest = "plotDir",
                  help = "Path to plot directory", metavar = "DIRECTORY")
parser.add_option("-l", "--lowData", dest = "lowData",
                  help = "Analysis uses VALUE% of the lowest points as ordered \
by secondary channel amplitude", metavar = "VALUE")
parser.add_option("-u", "--upperData", dest = "highData",
                  help = "Analysis uses VALUE% of the highest points as ordered\
by secondary channel amplitude. If both '-u' and '-l' options are used, will \
use overlapping points.", metavar = "VALUE")
parser.add_option("-x", "--xMax", dest = "xMaxGen",
                  help = "Allows manual selection of x axis plot maximum",
                  metavar = "VALUE")
parser.add_option("-y", "--yMax", dest = "yMaxGen",
                  help = "Allows manual selection of y axis plot maximum. Use \
'L' as VALUE to focus on fit line.", metavar = "VALUE")
parser.add_option("-i", "--intervals", dest = "timeIntervals",
                  help = "Allows data to be broken up into time intervals of \
length VALUE seconds. Can enter 'W' to break time up into week long chunks \
to be analyzed. Can enter path to file to break time up by time sections in \
file. (file version not currently implemented)", metavar = "VALUE")
parser.add_option("-q", "--iterator1", dest = "i",
                  help = "First iterator", metavar = "VALUE")
parser.add_option("-w", "--iterator2", dest = "n",
                  help = "Second iterator", metavar = "VALUE")

(options, args) = parser.parse_args()

# Number of Bins
NUMBER_BINS = 20

# load conf file
file = open(options.pickledfile, "rb")
conf_names = pickle.load(file)
file.close()

# load band tracker file
bandFile = open(options.pickledbands, "rb")
band_name_dict = pickle.load(bandFile)
bandFile.close()

# load times file
with open(options.pickledtimes, "rb") as infile:
    time_sections = pickle.load(infile)

# load frame files
channelData = load_frames(conf_names)

# some kind of check to make sure that all the times overlap

# run analysis
#analyzedData = blrms_upconv_analysis(channelData, band_name_dict, options)
#sections_make_plots(channelData, band_name_dict, options)
i = int(options.i)
n = int(options.n)
sections_make_plots_v2(channelData, band_name_dict, time_sections, options, i, n)
