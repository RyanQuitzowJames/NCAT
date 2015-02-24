from __future__ import division
from numpy import arange, array, ones, argsort, mean, sqrt
from scipy import stats
import pickle, os
from pylal.dq.dqFrameUtils import fromframefile
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from math import ceil, log10, floor
import DisplayPageLib as ndpage

# Notes: There could be issues if the times to be analyzed overlap themselves

# functions

def avg(data):
    if len(data) == 0:
        return 0
    else:
        return float(sum(data)) / len(data)

def testRank(slope, intercept, x_data, y_data, x_min, x_max, y_min, y_max):
    if slope < 0:
        return 0
    high_thresh = 1.2 * 10**-43
    high_float = (high_thresh - intercept) / slope
    total = 0
    for i in range(len(x_data)):
        if x_data[i] >= high_float and y_data[i] >= high_thresh:
            total += y_data[i]
            #total += 1
    rank = (total / len(x_data)) * -1
    return rank

def student_T_var(bin):
    mean = avg(bin)
    n = len(bin)
    if n <= 3:
        return avg(bin) / len(bin)**.5
    s2 = [((bin[i] - mean)**2) / (n - 1) for i in range(n)]
    s2 = sum(s2)
    s = s2**(.5)
    sigma = s / (n**.5)
    r = n - 1
    variance = (sigma**2) * (r / (r - 2))
    return variance

def chiSqrdFit(x, y, sigmas):
    N = len(x)
    A = sum([x[i] / sigmas[i]**2 for i in range(N)])
    B = sum([1 / sigmas[i]**2 for i in range(N)])
    C = sum([y[i] / sigmas[i]**2 for i in range(N)])
    D = sum([x[i]**2 / sigmas[i]**2 for i in range(N)])
    E = sum([x[i] * y[i] / sigmas[i]**2 for i in range(N)])
    F = sum([y[i]**2 / sigmas[i]**2 for i in range(N)])

    print("D = " + str(D))
    print("B = " + str(B))
    print("A = " + str(A))

    slope = (E * B - C * A) / (D * B - A**2)
    intercept = (D * C - E * A) / (D * B - A**2)

    sigmaSlope = (B / (B * D - A**2))**.5
    sigmaIntercept = (D / (B * D - A**2))**.5

    #print(sigmas)

    X2_list = [((y[i] - slope * x[i] - intercept)**2) / sigmas[i]**2 for i in range(N)]
    X2 = sum(X2_list)
    v = N - 2

    reduced_X2 = X2 / v
    p = reduced_X2

    print('X2 = ' + str(X2))
    print('v = ' + str(v))
    print("Reduced X2 = " + str(reduced_X2))
    print('Chi p = ' + str(p))

    return slope, intercept, reduced_X2, p, sigmaSlope, sigmaIntercept

def stdDev(data):
    n = float(len(data))
    if n <= 1:
        return 0
    mean = sum(data) / n
    diffs = [dataPoint - mean for dataPoint in data]
    diffsSquared = [dataPoint**2 for dataPoint in diffs]
    sumOfSquares = sum(diffsSquared)
    sigma2 = sumOfSquares / (n - 1)
    sigma = sigma2**(.5)
    return sigma

# Helper function to load a number as int or float as appropriate
def load_number(number):
    if int(number)/float(number) == 1:
        return int(number)
    else:
        return float(number)

# Helper function to make new directory
def create_dir(name, iterate_name = True):

    # set default directory name
    newDir = name
    # If directory doesn't exist, create
    if not os.path.exists(name):
        os.makedirs(name)

    # Otherwise, if iterate_name is set to true, iterate version number
    # to create new directory
    elif iterate_name:
        # Set initial version number
        version = 2
        # set base name to add version number to
        base_name = name + "_v"
        # while directory exists, iterate version number
        while os.path.exists(base_name + str(version)):
            version += 1
        # overwrite directory name
        newDir = base_name + str(version)
        # make new directory
        os.makedirs(newDir)

    return newDir

# helper function to time order data FLAGGGEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
def order_data(raw_channel_data):
    times = raw_channel_data[0]
    amplitude = raw_channel_data[1]
    indices = argsort(times)
    new_channel_data = [[],[]]
    new_channel_data[0] = [times[x] for x in indices]
    new_channel_data[1] = [amplitude[x] for x in indices]
    return new_channel_data

# same as above but hypothetically gaurds against duplicates in the data set

## test this funciton ##############################
def order_data_and_guard(raw_channel_data):
    times = raw_channel_data[0]
    amplitude = raw_channel_data[1]
    times2 = list(set(times))
    indices1 = [times.index(x) for x in times2]
    times = [times[x] for x in indices1]
    amplitude = [amplitude[x] for x in indices1]
    indices = argsort(times)
    new_channel_data = [[],[]]
    new_channel_data[0] = [times[x] for x in indices]
    new_channel_data[1] = [amplitude[x] for x in indices]
    return new_channel_data
#####################################################

# helper function to load from blrms2frame output frames
# add thing in end to sort entries by time.
def load_blrms_frame(frame_dir, filename, channel):
    print("make sure to finish this one!")#debug note!
    # get time from filename
    slice1 = filename[filename.find("-") + 1:]
    slice2 = slice1[slice1.find("-") + 1:]
    #float(x) if '.' in x else int(x)
    # need to change this if non integer numbers are common maybe.
    #frame_start = float(slice2[:slice2.find("-")]) # maybe int instead?
    frame_start = int(slice2[:slice2.find("-")])
    file_path = frame_dir + "/" + filename
    raw_channel = fromframefile(file_path, channel)
    processed_channel = [[],[]]
#    print("start_time:")
#    print(frame_start)#debug
    for num in range(len(raw_channel[0])):
        processed_channel[1].append(raw_channel[1][num])
        time = raw_channel[0][num] + frame_start
        processed_channel[0].append(time)
    processed_channel = order_data_and_guard(processed_channel)
    return processed_channel

# load channels
def load_channels(channel_list, file_name_list, frame_dir, channel_data = {}):

    for channel in channel_list:
        temp_data = [[],[]]
        for file in file_name_list:
            channel_temp = load_blrms_frame(frame_dir, file, channel)
            for num in range(len(channel_temp[0])):
                temp_data[0].append(channel_temp[0][num])
                temp_data[1].append(channel_temp[1][num])
        channel_data[channel] = order_data_and_guard(temp_data)
    return channel_data

# load frames
def load_frames(config_dict):
    channel_data = {}
    # loop through and find each output directory
    for base_name in config_dict:
        # get handle to frame directory
        frame_dir = config_dict[base_name]["frame directory"]
        # get list of all files in directory
        content = os.listdir(frame_dir)
        # initialize list to hold frame file names from directory
        frameFiles = []
        # loop through and add only frames to list
        for item in content:
            if ".gwf in item":
                frameFiles.append(item)
        # load channel data from frames
        channel_data = load_channels(config_dict[base_name]["band names"],
                                     frameFiles, frame_dir, channel_data)
    return channel_data

def fill_bins(x_data, y_data, number_bins = 20):
    x_min = min(x_data)
    x_max = max(x_data)
    x_delta = (x_max - x_min) / number_bins
    x_bins = [[x for x in x_data if x_min + num * x_delta <= x < x_min + (num + 1) * x_delta] for num in range(number_bins)]
    x_bins[-1].append(x_max)
    y_bins = [[y_data[ind] for ind, x in enumerate(x_data) if x_min + num * x_delta <= x < x_min + (num + 1) * x_delta] for num in range(number_bins)]
    y_bins[-1].append(y_data[x_data.index(x_max)])
    # Think more about edge cases, data could have multiple copies of the same value
    print('x_bins length: ' + str(len(x_bins)) + '   y_bins length: ' + str(len(y_bins)))
    return x_bins, y_bins

# grab overlapping times
def get_overlap(data_1, data_2):
    # find if any times happen multiple times? maybe skip this

    # find overlapping times
    ind_dict_1 = dict((k,i) for i,k in enumerate(data_1[0]))
    ind_dict_2 = dict((k,i) for i,k in enumerate(data_2[0]))
    shared = set(data_1[0]).intersection(data_2[0])
    indices_1 = [ind_dict_1[x] for x in shared]
    indices_2 = [ind_dict_2[x] for x in shared]
    new_data_1 = [[data_1[0][x] for x in indices_1],
                  [data_1[1][x] for x in indices_1]]
    new_data_2 = [[data_2[0][x] for x in indices_2],
                  [data_2[1][x] for x in indices_2]]
    return new_data_1, new_data_2

# Helper function to find the minimal bound line
def find_minimal_bound(x_data, y_data, anChan, secChan, x_min, x_max, y_min, y_max, number_bins = 20):
    print("is running: minimal bound thingy")
    # standard number of sections should be 20
    #print("x_data")#debug
    #print(x_data)#debug
    '''indices = argsort(x_data)
    #print(indices)#debug
    sortedX = []
    sortedY = []
    for index in indices:
        sortedX.append(x_data[index])
        sortedY.append(y_data[index])
    #print(len(sortedX))#debug
    deltaVal = (sortedX[-1] - sortedX[0]) / number_sections
    minList = [[],[]]
    currentIndex = 0
    for step in range(number_sections):
        quit_val = False
        max_test = max(y_data) + 1
        step_val = max_test
        pos_val = 0
        start_step_index = currentIndex
        while not quit_val:
            if currentIndex > len(x_data) -1 or \
               sortedX[currentIndex] > (step + 1) * deltaVal + sortedX[0]:
                quit_val = True
            else:
                if sortedY[currentIndex] < step_val:
                    step_val = sortedY[currentIndex]
                    pos_val = sortedX[currentIndex]
                currentIndex += 1
        if currentIndex != start_step_index:
            minList[0].append(pos_val)
            minList[1].append(step_val)'''

    x_bins, y_bins = fill_bins(x_data, y_data, number_bins)
    #x_max = max(x_data)
    #x_min = min(x_data)
    domain = x_max - x_min
    interval_width = domain / float(number_bins)
    new_x = [x_min + (i + 1) * interval_width - interval_width / 2 for i in range(number_bins) if len(y_bins[i]) > 0]
    new_y = [min(binny) for binny in y_bins if len(binny) > 0]

    #new_x = [new_x[ind] for ind, val in enumerate(y_bins) if len(val) > 0]
    #new_y = [new_y[ind] for ind, val in enumerate(y_bins) if len(val) > 0]

    sigmas = [avg(y_bins[i]) / sqrt(len(y_bins[i])) for i in range(len(y_bins)) if len(y_bins[i]) > 0]
    #scaled_sigmas = [avg(y_bins[i]) / sqrt(len(y_bins[i])) for i in range(len(y_bins)) if len(y_bins[i]) > 0]
    scaled_sigmas = array(sigmas)

    new_data = [new_x, new_y, scaled_sigmas]

    slope, intercept, reduced_X2, p_value, sigmaSlope, sigmaIntercept = chiSqrdFit(new_x, new_y, sigmas)
    std_err = [sigmaSlope, sigmaIntercept] 
    p_value = testRank(slope, intercept, x_data, y_data, x_min, x_max, y_min, y_max)
    print("Minimal Bound p = " + str(p_value))

    '''slope, intercept, r_value, p_value, std_err = \
           stats.linregress(minList[0],minList[1])'''
    min_linear_fit = [anChan, secChan, slope, intercept, reduced_X2, p_value, \
                      std_err, new_data]
    return min_linear_fit

def find_sigma_bound(x_data, y_data, anChan, secChan, sigmas = 2, number_bins = 20):
    x_bins, y_bins = fill_bins(x_data, y_data, number_bins)
    test = [len(x) for x in x_bins]
    print("bin populations: " + str(test))
    print('X Data = ' + str(x_data))
    print('Y Data = ' + str(y_data))
    x_max = max(x_data)
    x_min = min(x_data)
    domain = x_max - x_min
    interval_width = domain / float(number_bins)
    new_x = [x_min + (i + 1) * interval_width - interval_width / 2 for i in range(number_bins)]
    new_y = [avg(binny) - sigmas * stdDev(binny) for binny in y_bins]# if binny]
    #print("first bin: " + str(y_bins[0]))
    #print("means: " + str([mean(binny) for binny in y_bins]))
    #print("new_x: " + str(new_x))
    #print("new_y: " + str(new_y))
    #print('Magic')
    new_x = [new_x[ind] for ind, val in enumerate(y_bins) if len(val) > 1]
    new_y = [new_y[ind] for ind, val in enumerate(y_bins) if len(val) > 1]
    sigma_list = [(avg(y_bins[i]) - sigmas * stdDev(y_bins[i])) / len(y_bins[i])**.5 for i in range(len(y_bins)) if len(y_bins[i]) > 1]
    print("new_x: " + str(new_x))
    print("new_y: " + str(new_y))
    print("sigma List = " + str(sigma_list))
    new_data = [new_x, new_y]

    slope, intercept, reduced_X2, p_value, sigmaSlope, sigmaIntercept = chiSqrdFit(new_x, new_y, sigma_list)
    std_err = [sigmaSlope, sigmaIntercept]

    '''slope, intercept, r_value, p_value, std_err = \
           stats.linregress(new_x, new_y)'''
    sigma_fit = [anChan, secChan, slope, intercept, reduced_X2, p_value, \
                      std_err, new_data]
    return sigma_fit

# make linear plot or minimum correlated line plot
'''
def lin_min_plot(pl_name, plot_dir, x_name, y_name, x_data, y_data, slope,
                 intercept, command_args, bound_x = None, bound_y = None):
    # Get stuff from command line arguments
 #   pltXmin = None
#    pltXmax = None
#    pltYmin = None
#    pltYmax = None
#    if command_args.yMaxGen.upper() == "L":
#        pltYmin =
    # check proposed plot name for decimal, if found replace with 'dot'
    if '.' in pl_name:
        pl_name = pl_name.replace('.', 'dot')
    sqrdX = []
    sqrdY = []
    for num in range(len(x_data)):
        sqrdX.append(x_data[num]*x_data[num])
        sqrdY.append(y_data[num]*y_data[num])
    xMax = max(sqrdX)
    xMin = min(sqrdX)
    line = [[xMin, xMax],[xMin * slope + intercept, xMax * slope + intercept]]
    if bound_x and bound_y:
        plt.plot(sqrdX, sqrdY, 'bs', line[0], line[1], 'r', bound_x, \
                 bound_y, 'y*')
    else:
        plt.plot(sqrdX, sqrdY, 'b.', line[0], line[1], 'g--')
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    if command_args.yMaxGen:
        if command_args.yMaxGen.upper() == "L":
            plt.ylim([xMin * slope + intercept, xMax * slope + intercept])
        else:
            plt.ylim([0, float(command_args.yMaxGen)])
    if command_args.xMaxGen:
        plt.xlim([0, float(command_args.xMaxGen)])
    #if bound_x:
        #plt.ylim([min(bound_y),max(bound_y)])
        #plt.xlim([min(bound_x),max(bound_x)])
    #plt.ylim([0,1e-42])#xMin * slope + intercept, xMax * slope + intercept])#debug
    plt.savefig(plot_dir + "/" + pl_name)
    plt.clf()
'''

def vsPlot(analysisChannel, secondaryChannel, channelDataDictionary, plotDir, command_args, startTime = None, endTime = None, mainPlot = True):
    anChan = channelDataDictionary[analysisChannel]
    secChan = channelDataDictionary[secondaryChannel]
    if not anChan[0] == secChan[0]:
        anChan, secChan = get_overlap(anChan, secChan)

    x_data = secChan[1]
    y_data = anChan[1]

    if startTime and endTime:
        t_data = anChan[0]
        indices = [index for index, val in enumerate(t_data) if startTime < val < endTime]
        t_data = [t_data[x] for x in indices]
        x_data = [x_data[x] for x in indices]
        y_data = [y_data[x] for x in indices]
    elif startTime or endTime:
        print("Missing either a startTime or endTime input, please enter both or none. Treating data as if neither have been entered.")

    # Get command line arguments
    lowData = None
    highData = None
    numPoints = len(x_data)
    if command_args.lowData:
        lowData = float(command_args.lowData)/100
    if command_args.highData:
        highData = float(command_args.highData)/100
    if lowData and highData:
        if lowData <= 1 - highData:
            print("Data limits don't make sense.")
            print("lowData")
            print(lowData)
            print("highData")
            print(highData)
    if lowData:
        max_point = int(numPoints*lowData)
        x_data = x_data[:max_point]
        y_data = y_data[:max_point]
    if highData:
        min_point = int(numPoints*(1-highData))
        x_data = x_data[min_point:]
        y_data = y_data[min_point:]

    x_sqrd = [x**2 for x in x_data]
    y_sqrd = [y**2 for y in y_data]
    #for calculating fit line
    xMax = max(x_sqrd)
    xMin = min(x_sqrd)
    yDataMax = max(y_sqrd)
    yDataMin = min(y_sqrd)

    analysisChannel = edit_label(analysisChannel)
    secondaryChannel = edit_label(secondaryChannel)
    anChanName = analysisChannel.replace(".","dot")
    secChanName = secondaryChannel.replace(".","dot")
    plotName = plotDir + "/" + secChanName + " vs " + anChanName

    # calculate fit line
    minimalBoundData = find_minimal_bound(x_sqrd, y_sqrd, analysisChannel,                      secondaryChannel, xMin, xMax, yDataMin, yDataMax)
    #print('first assignment p = ' + str(minimalBoundData[5]))
    '''if (mainPlot and minimalBoundData[5] >= .85):
        return [minimalBoundData[5], plotName + ".png"]'''

    minimalBoundLine = [[xMin, xMax],[xMin * minimalBoundData[2] \
        + minimalBoundData[3], xMax * minimalBoundData[2] + minimalBoundData[3]]]
    # Sigma Fit Line
#    sigmaLine2 = find_sigma_bound(x_data, y_data, analysisChannel, secondaryChannel, 2)
#    sigmaLine3 = find_sigma_bound(x_data, y_data, analysisChannel, secondaryChannel, 3)
    #sigmaBound2 = find_sigma_bound(x_sqrd, y_sqrd, analysisChannel, secondaryChannel, 2)
    #sigmaLine2 = [[xMin, xMax],[xMin * sigmaBound2[2] + sigmaBound2[3], xMax * sigmaBound2[2] + sigmaBound2[3]]]
    #sigmaBound3 = find_sigma_bound(x_sqrd, y_sqrd, analysisChannel, secondaryChannel, 3)
    #sigmaLine3 = [[xMin, xMax],[xMin * sigmaBound3[2] + sigmaBound3[3], xMax * sigmaBound3[2] + sigmaBound3[3]]]
    # get normal fit line info
    #print(len(x_sqrd))
    #print(len(y_sqrd))
    slope, intercept, r_value, p_value, std_err = \
                   stats.linregress(x_sqrd, y_sqrd)
    line = [[xMin, xMax],[xMin * slope + intercept, xMax * slope + intercept]]
    # make plot
    # raw data
    plt.plot(x_sqrd, y_sqrd, 'c.')#, label="test")
    # simple line
####
####precision to round labels in graph to (# of significant figures)
    precision = 3
####
####
    #plt.axvline(x=200,color='g')
    #plt.axhline(y=1.2 * 10**-43,color='g')
    #LIIIIIINNNEEEEEAAAARRRRRRR FIT
    #plt.plot(line[0], line[1], 'g', label="Linear Fit\nslope = " + rounding(slope,precision) + "\np = " + rounding(p_value,precision) + "\nr^2 = " + str(rounding(r_value**2,precision)),alpha=0.5)
    # minimal bound line
    plt.plot(minimalBoundLine[0], minimalBoundLine[1], 'r', label="Minimal Bound\nslope = " + rounding(minimalBoundData[2],precision) + ' +- ' + rounding(minimalBoundData[6][0],precision) + "\np = " + rounding(minimalBoundData[5],precision) + "\nReduced X^2 = " + rounding(minimalBoundData[4],precision),alpha=0.5)
    # minimal bound data
    #plt.plot(minimalBoundData[7][0], minimalBoundData[7][1], 'r.')#, label = "test")
    #ERROR BARS
    #plt.errorbar(minimalBoundData[7][0], minimalBoundData[7][1], yerr = minimalBoundData[7][2], fmt = 'ro', mfc = 'red', mec = 'red', ecolor = 'red')#, label = "test")
    # 2 Sigma Line
    #print("sigma debug")
    #print(sigmaLine2[0])
    #print(sigmaLine2[1])
    #print("sigma debug more stuff")
    #print(sigmaBound2)
    #plt.plot(sigmaLine2[0], sigmaLine2[1], 'k', label="2 Sigma Bound\nslope = " + rounding(sigmaBound2[2],precision) + "\np = " + rounding(sigmaBound2[5],precision) + "\nReduced X^2 = " + rounding(sigmaBound2[4],precision),alpha=0.5)
    # 2 Sigma Data
    #plt.plot(sigmaBound2[7][0], sigmaBound2[7][1], 'k.')#, label = "test")
    # 3 Sigma Line
    #plt.plot(sigmaLine3[0], sigmaLine3[1], 'b', label="3 Sigma Bound\nslope = " + rounding(sigmaBound3[2],precision) + "\np = " + rounding(sigmaBound3[5],precision) + "\nReduced X^2 = " + rounding(sigmaBound3[4],precision),alpha=0.5)
    # 3 Sigma Data
    #plt.plot(sigmaBound3[7][0], sigmaBound3[7][1], 'b.')#, label = "test")

    if command_args.yMaxGen:
        if command_args.yMaxGen.upper() == "L":
            stdX = stdDev(x_sqrd)
            stdY = stdDev(y_sqrd)
            x_avg = avg(x_sqrd)
            y_avg = avg(y_sqrd)
            x_low = max(xMin, x_avg - 2 * stdX, 0)
            x_high = min(xMax, x_avg + 3 * stdX)
            #x_low = xMin
            #x_high = xMax

#yMinLinear = xMin * slope + intercept
            yMinMinimal = x_low * minimalBoundData[2] + minimalBoundData[3]
            ####yMinSigma2 = xMin + sigmaBound2[2] + sigmaBound2[3]
            ####yMinSigma3 = xMin + sigmaBound3[2] + sigmaBound3[3]
            #yMaxLinear = xMax * slope + intercept
            yMaxMinimal = x_high * minimalBoundData[2] + minimalBoundData[3]
            ####yMaxSigma2 = xMax * sigmaBound2[2] + sigmaBound2[3]
            ####yMaxSigma3 = xMax * sigmaBound3[2] + sigmaBound3[3]
            yMin = min(yMinMinimal, yMaxMinimal)
            yMin = max(0, yMin)
            yMax = max(yMinMinimal, yMaxMinimal)
            deltaY = yMax - yMin
            
            if yDataMax > yMax:
                yMax = min(3 * deltaY + yMin, yDataMax)
            yMax = max(yMax, y_avg)         

            plt.xlim([x_low, x_high])

            plt.ylim([yMin, yMax])
            #plt.ylim([yMin, 50 * 10**-43])

            #plt.ylim([xMin * slope + intercept, xMax * slope + intercept])
        else:
            plt.ylim([0, float(command_args.yMaxGen)])
    if command_args.xMaxGen:
        plt.xlim([0, float(command_args.xMaxGen)])

    # labels, title, gridlines

    plt.xlabel(secondaryChannel)
    plt.ylabel(analysisChannel)
    plt.title(secondaryChannel + " vs " + analysisChannel, fontsize = 14, y = 1.05)
    plt.grid(b=True, which='major',color='0.75',linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor',color='0.85',linestyle='--')
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits = (-2,2))
    plt.ticklabel_format(axis = 'x', style = 'sci', scilimits = (-2,2))
    legend = plt.legend(prop={'size':6})#, framealpha = 0.5)
    legend.get_frame().set_alpha(0.5)
    print(legend)
    # save and close figure
    '''anChanName = analysisChannel.replace(".","dot")
    secChanName = secondaryChannel.replace(".","dot")
    plotName = plotDir + "/" + secChanName + " vs " + anChanName'''
    plt.savefig(plotName, bbox_inches='tight')
#    plt.savefig(plotName, bb_inches = 'tight')
    plt.clf()
    plotName += ".png"
    #print("testing yeah")
    #print('VS Plot p = ' + str(minimalBoundData[5]))
    return [minimalBoundData[5], plotName]

# make a normal time plots
'''
def time_plot(pl_name, plot_dir, y_name, data, sqrd = False):
    x_data = data[0]
    y_data = data[1]
    if '.' in pl_name:
        pl_name = pl_name.replace('.','dot')
    if sqrd:
        dataX = []
        dataY = []
        for num in range(len(x_data)):
            dataX.append(x_data[num]*x_data[num])
            dataY.append(y_data[num]*y_data[num])
    else:
        dataX = x_data
        dataY = y_data
    plt.plot(dataX, dataY, 'bo')
    plt.xlabel("time GPS (s)")
    plt.ylabel(y_name)
    plt.yscale("log")
    plt.savefig(plot_dir + "/" + pl_name)
    plt.clf()
'''

'''
def timePlot(analysisChannel, secondaryChannel, channelDataDictionary, plotDir):
    t_data = channelDataDictionary[secondaryChannel][0]
    x_data = channelDataDictionary[secondaryChannel][1]
    y_data = channelDataDictionary[analysisChannel][1]
    plt.figure()
    plt.subplot(211)
    secondaryChannel = edit_label(secondaryChannel)
    analysisChannel = edit_label(analysisChannel)
    plt.title(secondaryChannel + " and\n" + analysisChannel, fontsize = 16, y = 1.03)
    plt.plot((t_data - t_data[0])/3600, y_data, 'b.-')
    plt.xlim(xmin = -.1)
    plt.grid(b=True, which='major',color='0.75',linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor',color='0.85',linestyle='--')
#    plt.xlabel(secondaryChannel)


    plt.ylabel(analysisChannel, fontsize = 10)
#    plt.title(secondaryChannel + " vs " + analysisChannel)
    plt.subplot(212)
    plt.plot((t_data - t_data[0])/3600, x_data, 'r.-')
    plt.xlim(xmin = -.1)
#    plt.xlabel(secondaryChannel)


    plt.ylabel(secondaryChannel, fontsize = 10)
    plt.xlabel("time (hr) ( +" + str(t_data[0]) + " [GPS])")
    plt.grid(b=True, which='major',color='0.75',linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor',color='0.85',linestyle='--')
    anChanName = analysisChannel.replace(".","dot")
    secChanName = secondaryChannel.replace(".","dot")
    plotName = plotDir + "/" + secChanName + " vs " + anChanName
    plt.savefig(plotName)#plotDir + "/" + secondaryChannel + " vs " + analysisChannel)
    plt.clf()
'''

def timeSqrdPlot(analysisChannel, secondaryChannel, channelDataDictionary, plotDir, startTime = None, endTime = None):
    t_data = channelDataDictionary[secondaryChannel][0]
    x_data = channelDataDictionary[secondaryChannel][1]
    y_data = channelDataDictionary[analysisChannel][1]

    t_zero = t_data[0]
    
    if not startTime and not endTime:
        x_length = len(x_data)
        y_length = len(y_data)
        min_length = min(x_length, y_length)
        if y_length > x_length:
            diff = y_length - x_length
            y_data = y_data[:-diff]
        elif x_length > y_length:
            diff = x_length - y_length
            x_data = x_data[:-diff]
    elif startTime and endTime:
        indices = [index for index, val in enumerate(t_data) if startTime < val < endTime]
        t_data = [t_data[x] for x in indices]
        x_data = [x_data[x] for x in indices]
        y_data = [y_data[x] for x in indices]
    elif startTime or endTime:
        print("Missing either a startTime or endTime input, please enter both or none. Treating data as if neither have been entered.")

    x_sqrd = [x**2 for x in x_data]
    y_sqrd = [y**2 for y in y_data]
    analysisChannel = edit_label(analysisChannel)
    secondaryChannel = edit_label(secondaryChannel)
    plt.figure(figsize=(9,6))
    if startTime and endTime:
    #    plt.figure(figsize=(9,7))
        title = secondaryChannel + " and " + analysisChannel + "\nFor " + str(startTime) + " - " + str(endTime) + " GPS"
    else:
      #  plt.figure()
        title = secondaryChannel + " and " + analysisChannel
####plt.subplot(211)
#    plt.title(title, y=1.1, fontsize = 18)
    plt.title(title, fontsize=18, y = 1.05)
#    plt.plot((t_data-t_data[0])/3600, y_sqrd, 'b-')
#    plt.xlim(xmin=-.1)
    #plt.plot((t_data-t_data[0])/60, y_sqrd, 'b-')
    #plt.xlim(xmin=-1/3)
####plt.plot((t_data - t_zero)/60, y_sqrd, 'b-')
    plt.xlim(xmin = (t_data[0] - t_zero)/60 -1/3)

#    plt.ylim([0,0.001])
    plt.grid(b=True, which='major',color='0.75',linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor',color='0.85',linestyle='--')
#    plt.xlabel(secondaryChannel)

    #edits channel labels - Paul
#    analysisChannel = edit_label(analysisChannel)

    ax1=plt.subplot(111)
    plt.ticklabel_format(axis = 'y', style = 'sci', scilimits = (-2,2))

####plt.ylabel(analysisChannel + "\n (amplitude squared)", fontsize = 12)
    ax1.set_ylabel(analysisChannel + "\n (amplitude squared)", fontsize = 12, color='red')
    for t1 in ax1.get_yticklabels():
        t1.set_color('r')
#    plt.title(secondaryChannel + " vs " + analysisChannel)
    #plt.subplot(212)
#   plt.plot((t_data-t_data[0])/3600, x_sqrd, 'r-')
#    plt.xlim(xmin=-.1)
    ax1.plot((t_data - t_zero)/60, x_sqrd, 'r-')
    #plt.xlim(xmin = (t_data[0] - t_zero)/60 -1/3)
#    plt.xlabel(secondaryChannel)
#    plt.ylim([0,0.001])

    #added editing of labels
#    secondaryChannel = edit_label(secondaryChannel)

####plt.ylabel(secondaryChannel + "\n (amplitude squared)", fontsize = 12)
    plt.xlabel("time (min) + (" + str(t_zero) + " [GPS])")
#    plt.xlabel("time [hr]  " + "(+ " + str(t_data[0]) + " GPS)", fontsize = 14)
    plt.grid(b=True, which='major',color='0.75',linestyle='-')
    plt.minorticks_on()
    plt.grid(b=True, which='minor',color='0.85',linestyle='--')
    #ax2=plt.subplot(212)
    ax2=ax1.twinx()
    ax2.plot((t_data - t_zero)/60, y_sqrd, 'b-')
    ax2.ticklabel_format(axis = 'y', style = 'sci', scilimits = (-2,2))
    ax2.set_ylabel(secondaryChannel + "\n (amplitude squared)", fontsize = 12, color='blue')
    for t1 in ax2.get_yticklabels():
        t1.set_color('b')
    anChanName = analysisChannel.replace(".","dot")
    secChanName = secondaryChannel.replace(".","dot")
    plotName = plotDir + "/" + secChanName + " vs " + anChanName
    plt.savefig(plotName, bbox_inches='tight')#plotDir + "/" + secondaryChannel + " vs " + analysisChannel)
#    plt.savefig(plotName, bb_inches = 'tight')#plotDir + "/" + secondaryChannel + " vs " + analysisChannel)
    plt.clf()

# run analysis
def make_plots(channel_data, band_name_dict, command_args):
    # get list of channels
    anChanList = band_name_dict["analysischannels"]
    secChanList = band_name_dict["secondarychannels"]

    # create directories for each plot type
    # get plot directory from command line arguments
    plot_dir = command_args.plotDir
    # create dir for each type of data analysis
    linPlotDir = create_dir(plot_dir + "/basicLinearPlots")
    minBoundPlotDir = create_dir(plot_dir + "/minimalLinearBoundPlots")
    #timePlotDir = create_dir(plot_dir + "/timePlots")
    timeSqrdPlotDir = create_dir(plot_dir + "/timeSqrdPlots")

    # create plots

    # for plot ranking
    pStats = []

    for anChan in anChanList:
        for secChan in secChanList:
            tempStat = vsPlot(anChan, secChan, channel_data, minBoundPlotDir, command_args)
            pStats.append(tempStat)
            #timePlot(anChan, secChan, channel_data, timePlotDir)
            timeSqrdPlot(anChan, secChan, channel_data, timeSqrdPlotDir)
####                sqrdAnChan, sqrdSecChan = get_overlap(sqrdAnChan, sqrdSecChan)
###            if highData:
###                min_point = int(numberDataPoints*(1-highData))
###            if lowData:
###                max_point = int(numberDataPoints*lowData)# - 1
    rankStats = [[x[0] for x in pStats if not "nan" == str(x[0]) and x[0] < .85],[x[1] for x in pStats if not "nan" == str(x[0]) and x[0] < .85]]
    rankStats = order_data_and_guard(rankStats)
#    print(rankStats)
    with open(minBoundPlotDir + "/rankings.txt", "w") as outfile:
        for num in range(len(rankStats[0])):
            outfile.write("p = " + str(rankStats[0][num]) + ", " + rankStats[1][num] + "\n")
#    return pStats

# analysis function to handle multiple time section option for plots
def sections_make_plots(channel_data, band_name_dict, command_args):
    if not command_args.timeIntervals:
        make_plots(channel_data, band_name_dict, command_args)
    else:
        # temporary. better method needed
        max_time = 0
        min_time = 1000000000000000000000
        # order channel data
        for channel in channel_data:
            channel_data[channel] = order_data_and_guard(channel_data[channel])
            # find min time
            temp_time = channel_data[channel][0][0]
            if temp_time < min_time:
                min_time = temp_time
            # find max time
            temp_time = channel_data[channel][0][-1]
            if temp_time > max_time:
                max_time = temp_time
        # find number of sections
        time_length = command_args.timeIntervals
        if str(time_length).lower() == "w":
            time_length = 7*24*3600
        number_sections = ceil((max_time - min_time + 1)/time_length) # +1 to account for if the difference between last time and first time exactly matches a multiple of time_length
        # split data by time sections
        # initial values
        front = 0
        output_channel_data = {}
        for num in range(number_sections):
            start = min_time + num*time_length
            end = min_time + (num+1)*time_length
            output_channel_data[str(start) + "-" + str(end)] = {}

            temp_channel_data = {}
            # cycle through channels
            for channel in channel_data:#[channel]:
                temp_data = [[],[]]
                temp_data[0] = [channel_data[channel][0][x] for x in range(len(channel_data[channel][0])) if start <= channel_data[channe][0][x] < end]
                temp_data[1] = [channel_data[channel][1][x] for x in range(len(channel_data[channel][0])) if start <= channel_data[channe][0][x] < end]
                command_args_copy = command_args
                command_args_copy.plotDir += "/" + str(start) + "-" + str(end)
                make_plots(channel_data, band_name_dict, command_args_copy)

# run analysis
def make_plots_v2(channel_data, band_name_dict, times_list, command_args, i, n):
    # get list of channels
    anChanList = band_name_dict["analysischannels"]
    secChanList = band_name_dict["secondarychannels"]

    # create plots

    # for plot ranking
    pStats = []

    anChan = anChanList[i]
    secChan = secChanList[n]

    plot_dir = command_args.plotDir
    minBoundPlotDir = plot_dir + "/minimalLinearBoundPlots/"
    timeSqrdPlotDir = plot_dir + "/timeSqrdPlots/"

    multipleTimePlotDirs = None
    multipleVsPlotDirs = None
    if len(times_list) > 1:
        print(len(times_list))#debug
        #keys = [x[0] + " - " + x[1] + " GPS" for x in times_list]
        #multipleTimePlotDirs = dict((x[0] + " - " + x[1] + " GPS", x[0] + "-" + x[1]) for x in times_list)
        multipleTimePlotDirs = dict((create_dir(timeSqrdPlotDir + "/" + x[0] + "-" + x[1], False), x) for x in times_list)
        multipleVsPlotDirs = dict((create_dir(minBoundPlotDir + "/" + x[0] + "-" + x[1], False), x) for x in times_list)

    tempStat = vsPlot(anChan, secChan, channel_data, minBoundPlotDir, command_args)
    if tempStat[0] <= .85 or 1 > 0:
        pStats.append(tempStat)
    #            timePlot(anChan, secChan, channel_data, timePlotDir)
        timeSqrdPlot(anChan, secChan, channel_data, timeSqrdPlotDir)
        if multipleTimePlotDirs:
            for key in multipleTimePlotDirs:
                startTime = load_number(multipleTimePlotDirs[key][0])
                endTime = load_number(multipleTimePlotDirs[key][1])
                #timeSqrdPlot_v2(anChan, secChan, channel_data, key, startTime, endTime)
                timeSqrdPlot(anChan, secChan, channel_data, key, startTime, endTime)
        if multipleVsPlotDirs:
            for key in multipleVsPlotDirs:
                startTime = load_number(multipleVsPlotDirs[key][0])
                endTime = load_number(multipleVsPlotDirs[key][1])
                #vsPlot_v2(anChan, secChan, channel_data, key, command_args, startTime, endTime)
                vsPlot(anChan, secChan, channel_data, key, command_args, startTime, endTime, False)
    ####                sqrdAnChan, sqrdSecChan = get_overlap(sqrdAnChan, sqrdSecChan)
    ###            if highData:
    ###                min_point = int(numberDataPoints*(1-highData))
    ###            if lowData:
    ###                max_point = int(numberDataPoints*lowData)# - 1
        rankStats = [[x[0] for x in pStats if not "nan" == str(x[0])],[x[1] for x in pStats if not "nan" == str(x[0])]]
        rankStats = order_data_and_guard(rankStats)
        #print(rankStats)
        with open(minBoundPlotDir + "/rankings_" + str(i) + "_" + str(n) + ".txt", "w") as outfile:
            for num in range(len(rankStats[0])):
                outfile.write("p = " + str(rankStats[0][num]) + ", " + rankStats[1][num] + "\n")
        with open(timeSqrdPlotDir + "/subplots_" + str(i) + "_" + str(n) + ".txt", "w") as outfile:
            #directories = [key for key in multipleTimePlotDirs]
            if multipleTimePlotDirs:
                for key in multipleTimePlotDirs:
                    descriptor = str(multipleTimePlotDirs[key][0]) + "-" + \
                                 str(multipleTimePlotDirs[key][1]) + " GPS"
                    outfile.write(key + ", " + descriptor + "\n")
    '''
        with open(minBoundPlotDir + "/subplots.txt", "w") as outfile:
            #directories = [key for key in multipleTimePlotDirs]
            for key in multipleVsPlotDirs:
                descriptor = str(multipleVsPlotDirs[key][0]) + "-" + \
                             str(multipleVsPlotDirs[key][1]) + " GPS"
                outfile.write(key + ", " + descriptor + "\n")
    '''
    #    return pStats

# analysis function to handle multiple time section option for plots
def sections_make_plots_v2(channel_data, band_name_dict, times_list, command_args, i, n):
    if not command_args.timeIntervals:
        make_plots_v2(channel_data, band_name_dict, times_list, command_args, i, n)
    else:
        ####
        # Note: this section may no longer be compatible
        ####

        # temporary. better method needed
        max_time = 0
        min_time = 1000000000000000000000
        # order channel data
        for channel in channel_data:
            channel_data[channel] = order_data_and_guard(channel_data[channel])
            # find min time
            temp_time = channel_data[channel][0][0]
            if temp_time < min_time:
                min_time = temp_time
            # find max time
            temp_time = channel_data[channel][0][-1]
            if temp_time > max_time:
                max_time = temp_time
        # find number of sections
        time_length = command_args.timeIntervals
        if str(time_length).lower() == "w":
            time_length = 7*24*3600
        number_sections = ceil((max_time - min_time + 1)/time_length) # +1 to account for if the difference between last time and first time exactly matches a multiple of time_length
        # split data by time sections
        # initial values
        front = 0
        output_channel_data = {}
        for num in range(number_sections):
            start = min_time + num*time_length
            end = min_time + (num+1)*time_length
            output_channel_data[str(start) + "-" + str(end)] = {}

            temp_channel_data = {}
            # cycle through channels
            for channel in channel_data:#[channel]:
                temp_data = [[],[]]
                temp_data[0] = [channel_data[channel][0][x] for x in range(len(channel_data[channel][0])) if start <= channel_data[channe][0][x] < end]
                temp_data[1] = [channel_data[channel][1][x] for x in range(len(channel_data[channel][0])) if start <= channel_data[channe][0][x] < end]
                command_args_copy = command_args
                command_args_copy.plotDir += "/" + str(start) + "-" + str(end)
                make_plots(channel_data, band_name_dict, command_args_copy)

#function that edits a label; takes last two '_' to ' - ' and adds ' Hz' to end
def edit_label(inlabel):
    lb1 = inlabel[::-1]
    lb2 = lb1.replace('_',' - ',2)
    outlabel = lb2[::-1]
    outlabel = outlabel + ' Hz'
    return outlabel

def rounding(x, prec):
    if x==0:
        return str(0.0)
    elif str(abs(x)) == "nan":
        return "NaN"
    elif abs(x) == float("inf"):
        return "inf"
    else:
        y = abs(x)
        print('y = ' + str(y))
        z = -int(floor(log10(y)))
        return str(round(x,z+prec-1))
