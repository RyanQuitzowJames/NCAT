from optparse import OptionParser
import os
from NCATLib import *

# command line options
parser = OptionParser()
parser.add_option("-c", "--conf", dest = "configfile",
                  help = "Path to config file detailing analysis",
                  metavar = "FILE")
parser.add_option("-t", "--timefile", dest = "timefile",
                  help = "Path to file deatailing times to be analyzed",
                  metavar = "FILE")
parser.add_option("-d", "--dir", dest = "plotDir",
                  help = "Path to plot directory", metavar = "DIRECTORY")
parser.add_option("-b", "--endBuffer", dest = "endbuffer",
                  help = "Amount of time to drop from the end of time segments",
                  metavar = "TIME")
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
parser.add_option("-m", "--maxjobs", dest = "maxjobs",
                  help = "Maximum number of jobs ever submitted at once \
through condor", metavar = "NUMBER")
parser.add_option("-i", "--intervals", dest = "timeIntervals",
                  help = "Allows data to be broken up into time intervals of \
length VALUE seconds. Can enter 'W' to break time up into week long chunks \
to be analyzed. Can enter path to file to break time up by time sections in \
file. (file version not currently implemented)", metavar = "VALUE")

(options, args) = parser.parse_args()

# constants

    # GPS Range
GPS_lower = 800000000
GPS_upper = 1200000000

    # blrms2frame path
blrmsExc = "/home/quitzow/newBLRMS_9_19_2014/blrms2frame"

    # python script path
ncatPath = os.path.dirname(os.path.realpath(__file__))
pythonScript = ncatPath + "/PlotAnalysis.py"
#pythonScript = "/home/vincent.roma/public_html/upconversionCode/prototypeNCAT6/testScript"
# ULNATPlotAnalysis path
    # shell path
shellPath = "#!/bin/bash"
#shellPath = "#!/usr/bin/env"
# check for minimum commands line arguments to function
if not options.configfile or not options.timefile or not options.plotDir:
    print("\nMissing arguments: please specify at least a configuration file, a time \
list \nfile and an output plot directory to run this program.\n\n")
    quit_program = True
else:
    quit_program = False

# Get inputs

    # Get config file input
catDictionary, quit_program = get_input_from_file_V_C_line(quit_program,
                                                           options.configfile)

    # Get start and end times
if options.endbuffer:
    times, quit_program = load_time_file_V_C_line(quit_program, GPS_lower,
                                              GPS_upper, options.timefile,
                                                  options.endbuffer)
else:
    times, quit_program = load_time_file_V_C_line(quit_program, GPS_lower,
                                              GPS_upper, options.timefile)

# probably want to move copies of any weighting functions to the support file
# directory

# Build file system

if not quit_program:
    # Build base analysis directory
    if options.plotDir[-1] == "/":
        baseDir = dated_dir(options.plotDir + "blrms_analysis")
    else:
        baseDir = dated_dir(options.plotDir + "/blrms_analysis")
    print(baseDir)#debug
#    baseDir = os.getcwd() + "/" + dated_dir("blrms_analysis")#, iterate_name = False)

    # Change to base analysis directory
#    os.chdir(baseDir)
#basePath = os.getcwd()

    # Build support file sub directory
    supportDir = create_dir(baseDir + "/support_files")
    print(supportDir)#debug

    # Build directory to hold copies of input parameter files
    copyDir = create_dir(supportDir + "/input_files")

    # Build output frame file sub directory
    frameDir = create_dir(baseDir + "/frames")
    print(frameDir)#debug

    # Build plot sub directory
    plotDir = create_dir(baseDir + "/plots")
    print(plotDir)#debug

    # Build support file sub directory for frame lists
#    os.chdir(supportDir) # change it so this isn't neccessary by making where the config files are written take an absolute path.
    frameListDir = create_dir(supportDir + "/frame_lists")
    print(frameListDir)#debug

    # Build support file sub directory for dag logs
    #logDir = create_dir(supportDir + "/logs")
    dagLogDir = create_dir(supportDir + "/dagLogs")
    logDir = create_dir(supportDir + "/jobLogs")
    print(logDir)#debug
else:
    supportDir = None
    frameListDir = None
    frameDir = None
    plotDir = None

# copy input files
if not quit_program:
    copy_input_file(options.configfile, copyDir)
    copy_input_file(options.timefile, copyDir)

# Create config file for blrms2frame executable in support directory
#os.chdir(supportDir)
conf_names, quit_program, band_name_dict = create_config_file(catDictionary,
                                                              supportDir,
                                                              quit_program)
#os.chdir(basePath)

# Create frame file lists
# can change this so that the function writes directly to the directory
# regardless of the directory you're in
#if not quit_program:
#    os.chdir(frameListDir)
conf_names, quit_program = create_frame_file_list_files(conf_names, quit_program, times, frameListDir)

# Build frame sub directories
#conf_names, quit_program
create_frame_output_directories(conf_names, frameDir, quit_program)

# Pickle conf_name object
pickled_conf_names = pickle_object(conf_names, supportDir, "conf_names",
                                   quit_program)

# Pickle band names object
pickled_band_names = pickle_object(band_name_dict, supportDir, "band_names",
                                   quit_program)

# Pickle time segments for multiple time graphs
pickled_times = pickle_object(times, supportDir, "analysis_times",
                              quit_program)

plot_dir = plotDir

minBoundPlotDir = create_dir(plot_dir + "/minimalLinearBoundPlots")
timeSqrdPlotDir = create_dir(plot_dir + "/timeSqrdPlots")
bandFile = open(pickled_band_names, "rb")
band_name_dict = pickle.load(bandFile)
bandFile.close()
anChanList = band_name_dict["analysischannels"]
secChanList = band_name_dict["secondarychannels"]

file = open(pickled_times, "rb")
times_list = pickle.load(file)
if len(times_list) > 1:
    #print(len(times_list))#debug
    multipleTimePlotDirs = dict((create_dir(timeSqrdPlotDir + "/" + x[0] + "-" + x[1]), x) for x in times_list)
    multipleVsPlotDirs = dict((create_dir(minBoundPlotDir + "/" + x[0] + "-" + x[1]), x) for x in times_list)
file.close()

with open(plot_dir + "/totals.txt", 'w') as outfile:
    outfile.write(str(len(anChanList)) + " " + str(len(secChanList)))

# Create analysis shell script
shellScripts = create_shell_script(pickled_conf_names, pickled_band_names,
                                  pickled_times, supportDir, plotDir,
                                  shellPath, pythonScript, options,
                                  quit_program)
shellScript = shellScripts[0]
shellScript2 = shellScripts[1]

# Set shell to executable
os.chmod(shellScript, 0764)
os.chmod(shellScript2, 0764)# 0755)

# Create DAG
dagFile = create_dag(conf_names, pickled_conf_names, supportDir, plotDir, blrmsExc,
           shellScript, pickled_band_names, quit_program)

# Run DAG
run_dag(quit_program, dagFile, options.maxjobs)
