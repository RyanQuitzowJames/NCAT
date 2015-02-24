from optparse import OptionParser
#from PrototypeNCATPlotAnalysisLib_V18 import *
import DisplayPageLib as ndpage
import os


print('done importing')
# Notes: There could be issues if the times to be analyzed overlap themselves

# command line options
parser = OptionParser()
parser.add_option("-d", "--dir", dest = "plotDir",
                  help = "Path to plot directory", metavar = "DIRECTORY")
parser.add_option("-i", "--iterator1", dest = "i",
                  help = "First Iterator", metavar = "ITERATOR1")
parser.add_option("-n", "--iterator2", dest = "n",
                  help = "Second Iterator", metavar = "ITERATOR2")
(options, args) = parser.parse_args()
print('dont with options crap')

# make webpage
directory = options.plotDir
if directory[-1] == "/":
    directory = directory[:-1]
directory = directory[::-1]
directory = directory[directory.find("/")+1:]
directory = directory[::-1]

print('messed with directory')

iMax = int(options.i)
nMax = int(options.n)

print('iMax and nMax are defined')
rank_string = ""
sub_string = ""
for i in range(iMax):
    for n in range(nMax):      
        name = directory + "/plots/minimalLinearBoundPlots/rankings_" + str(i) + "_" + str(n) + ".txt"
        if os.path.exists(name):
            file = open(name, 'r')
            file_string = file.read()
            rank_string = rank_string + file_string + '\n'
            file.close()
        name2 = directory + "/plots/timeSqrdPlots/subplots_" + str(i) + "_" + str(n) + ".txt"
        if os.path.exists(name2):
            file2 = open(name2, 'r')
            file2_string = file2.read()
            sub_string = sub_string + file2_string + '\n'
            file2.close()
print('for loop complete')
file =  open(directory + "/plots/minimalLinearBoundPlots/rankings.txt", 'w')
file.write(rank_string)
file.close()
file2 =  open(directory + "/plots/timeSqrdPlots/subplots.txt", 'w')
file2.write(sub_string)
file.close()
print('wrote new files')
ndpage.make_display_page(directory)
