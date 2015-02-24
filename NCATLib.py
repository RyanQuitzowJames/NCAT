from __future__ import print_function
from sys import hexversion
import os, datetime, subprocess, pickle

def load_number(number):
    if int(number)/float(number) == 1:
        return int(number)
    else:
        return float(number)

# Helper function to copy input parameter files
def copy_input_file(filePath, outputDirectory):
    if "/" in filePath:
        reversePath = filePath[::-1]
        reverseName = reversePath[:reversePath.index("/")]
        fileName = reverseName[::-1]
    else:
        fileName = filePath
    if outputDirectory[-1] == "/":
        outputPath = outputDirectory + fileName
    else:
        outputPath = outputDirectory + "/" + fileName
    with open(filePath, "r") as infile:
        text = [line for line in infile]
    with open(outputPath,"w") as outfile:
        outfile.write("".join(line for line in text))

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

# Helper function to create dated directory
def dated_dir(name, date = None, iterate_name = True):

    # If date empty, get time and date
    if not date:
        date = datetime.datetime.now()
    # create dated name
    dated_name = name + "-" + str(date.year) + "_" + str(date.month) + \
                 "_" + str(date.day)
    # create directory
    newDir = create_dir(dated_name, iterate_name)

    return newDir

# Helper function to read in data from a file
# Possibly add a check here for whitespace characters that aren't spaces.  Think
# about how robust I want this to be and whether I want to check for this
# elsewhere instead.
def read_text_file(file_name, delimeter, strip_trailing = True):
    '''
    file = open(file_name, "r")
    data = []
    for line in file:
        temp_list = []
        list_index = -1
        for letter in line:
            list_index += 1
            if letter == delimeter:
                temp_list.append(list_index)
        if temp_list:
            if temp_list[-1] != list_index:
                if line[temp_list[-1]+1] != ' ':
                    if line.find('\n') == -1:
                        list_index += 1
                    temp_list.append(list_index)
        # for single word lines
        elif list_index > -1 and line[0] != '\n':
            if line[-1] == '\n':
                temp_list.append(list_index)
            else:
                temp_list.append(list_index + 1)
        current_index = -1
        row_data = []
        for number in temp_list:
            if number != current_index + 1:
                row_data.append(line[current_index+1:number])
            current_index = number
        if row_data:
            data.append(row_data)
    file.close()
    '''

    # is the empty variable needed?
    data = None
    with open(file_name, "r") as infile:
        if strip_trailing:
            data = [filter(None,line.strip().split(delimeter)) for line in infile if not line.isspace()]#(x.isspace or not x)]
        else:
            data = [line.strip().split(delimeter) for line in infile if not line.isspace()]#(x.isspace or not x)]
    return data

# Helper function to take take input properly in python 2 or 3
def version_input(input_string):
    output = ""
    if hexversion >= 0x3000000:
        output = input(input_string)
    else:
        output = raw_input(input_string)
    return output

# Helper function to ask a yes or no question
def ask_yes_no(question):
    answered = False
    while not answered:
        answer = version_input(question)
        if answer.lower() == 'y' or answer.lower() == 'n':
            answered = True
        else:
            print("\nSorry, '" + answer + "' is not a valid answer.  Please \
answer either 'y' or 'n'.\n")
    return answer.lower()

# Helper function to concatenate a list of strings
# Intended to be used to figure out which section of the config file we are
# looking at.  This is here so I don't have to worry about whether people have
# used spaces or not in the config file.
def concat_line(word_list):
    '''
    line = ''
    for word in word_list:
        line += word
    '''
    line = "".join(x for x in word_list)
    return line.lower()

# Helper function to sort data from list into lists
def sort_list (file_data, categories, catDictionary):
    category = ""
    # create for loop to go through the list
    for line in file_data:
        # Check if line is a category
        isCat = False
        for cat in categories:
            if concat_line(line) == cat:
                category = cat
                isCat = True
                catChanged = True
        # Input data into proper category
        if not isCat:
            if category == 'generalinformation':
                catDictionary[category].append(line)
            else:
                if line[0].lower() != 'band' and line[0].lower() != 'weight':
                    catChanged = False
                if catChanged:
                    catDictionary[category + 'Gen'].append(line)
                else:
                    catDictionary[category].append(line)
    return catDictionary

# Helper function to determine if a string is a number
def number_check(num):
    try:
        temp = float(num)
        if temp == float('NaN') or temp == float('inf'):
            return False
        else:
            return True
    except ValueError:
        return False

# Helper function to check to make sure non category line is entered correctly.
def check_line(line, category):

    # entry checks
    # interval or stride
    quit_val = False
    if (line[0] == "interval" or line[0] == "stride"):

        # it's in gen info
        if category != 'generalinformation':
            print("\n'" + line[0].title() + "' should only be a member of \
'General Information'.\nPlease edit file appropriately.")
            quit_val = True

        # it's followed by one number
        elif len(line) != 2 or number_check(line[1]) == False:
            print("\n'" + line[0].title() + "' must only be followed by one \
non-infinite real number.\nPlease edit file appropriately.")
            quit_val = True

    # channel
    elif line[0].lower() == 'channel':

        # followed by three strings
        if len(line) != 4:
            print("\n'Channel' must only be followed by channel name, frame \
type and observatory.  See line '" + " ".join(word for word in line) + " in '" \
+ category.title() + "' section.\nPlease edit file appropriately.")
            quit_val = True

    # band
    elif line[0].lower() == 'band':

        # followed by two numbers
        if len(line) != 3 or number_check(line[1]) == False or number_check(line[2]) == False:
            print("\n'Band' must only be followed by two non-infinite real \
numbers.  See line '" + " ".join(word for word in line) + "' in '" \
                  + category.title() + "' section.\nPlease edit file \
appropriately.")
            quit_val = True

    # weight
    elif line[0].lower() == 'weight':

        # followed by one string, and then followed by 2 numbers
        if len(line) != 4 or number_check(line[2]) == False or number_check(line[3]) == False:
            print("\n'Weight' must only be followed by a path to a file \
(including file name) and then two non-infinite real numbers.  See line '" \
                  + " ".join(word for word in line) + "' in '" \
                  + category.title() + "' section.\nPlease edit file \
appropriately.")
            quit_val = True

    elif not ((line[0] == "interval" or line[0] == "stride") and category == 'generalinformation'):
        print("\n'" + line[0].upper() + "' is not a valid entry.  Only 'Band' and \
'Weight are valid options, with the addtional options of 'stride and 'interval' \
in the 'General Information' category.\nPlease edit file appropriately.")
        quit_val = True
    return quit_val

# Helper function to write line to output string to be written to
# configuration files later.
def write_config_string(line, category, output_text, channel = None):

    # interval or stride
    if (line[0] == "interval" or line[0] == "stride") and not channel:
        output_text += line[0].upper() + " " + line[1] + "\n"
        if line[0] == 'interval':
            output_text += "\n"

    # band
    elif line[0].lower() == 'band' and channel != None:
        output_text += "\nBAND " + "_".join(n for n in line[1:]) + " " \
                       + channel + "".join(" " + n for n in line[1:])

    # weight
    elif line[0].lower() == 'weight' and channel != None:
        wf_name = line[1][:line[1].find('.')].lower().replace('weight', 'w')
        output_text += "\nBAND " + "_".join(n for n in line[2:]) + "_" \
                       + wf_name + " " + channel.upper() \
                       + "".join(" " + n for n in line[2:]) + " " + wf_name

    return output_text

# Helper function to handle config file dictionary keys
def config_key(channel_line):
    observatory = channel_line[3]
    frame_type = channel_line[2]
    key = observatory + "-" + frame_type
    #key = "-".join(k for k in channel_line[2:4])
    return key

# Helper function to write line to designate weight file
def record_weight_file(line, output_text):
    if line[0].lower() == 'weight':
        if line[1] not in output_text:
            filename = line[1]
            output_text += "WEIGHT " \
                + filename[:filename.find('.')].lower().replace('weight','w') \
                + " " + filename + "\n"
    return output_text

# Helper function to keep track of which channels that will be created
# by the blrms executable will be analysis channels or secondary
# channels, as well as track which frame file directories they will
# be found in.
def track_output_channel(line, channel, band_name_list):

    # entry checks
    # band
    if line[0].lower() == 'band':

        # output to specified band list for plot analysis later or to
        # track which config file they are related to.
        band_name = channel + "_" + line[1] + "_" + line[2]
        band_name_list.append(band_name)

    # weight
    elif line[0].lower() == 'weight':

        # output to specified band list for plot analysis later or to
        # track which config file they are related to.
        band_name = channel + "_" + line[2] + "_" + line[3]
        band_name_list.append(band_name)

    return band_name_list

# Helper function to record band lines, make sure weighting functions are
# entered, and track whether an output channel is an analysis channel or
# secondary channel, as well as which config file it's related to
def record_BAND(band, category, current_channel, band_name_list, config_file_dict, key):

    # create handle to band string entry
    band_strings = config_file_dict[key]["bands string"]

    # create handle to weighting function string entry
    weighting_functions = config_file_dict[key]["weighting functions string"]

    # create handle to band name entry
    band_names = config_file_dict[key]["band names"]

    # record BAND line
    band_strings = write_config_string(band, category, \
                                                 band_strings, current_channel)

    # check if a weighted band, and record weight function if not already
    # recorded
    weighting_functions = record_weight_file(band, weighting_functions)

    # feed output name of channel and band into the band_name_dict
    band_name_list = track_output_channel(band, current_channel, \
                                              band_name_list)

    # feed output name of channel and band into the config file entry
    band_names = track_output_channel(band, current_channel, band_names)

    # set band line string
    config_file_dict[key]["bands string"] = band_strings

    # set weighting function string
    config_file_dict[key]["weighting functions string"] = weighting_functions

    return config_file_dict, band_name_list

# Helper function to grab channels from given category and put all appropriate
# bands of that channel in output text.
def build_bands(category, catDict, quit_program, config_file_dict, \
                band_name_dict = {'analysischannels':[], \
                                  'secondarychannels':[]}):

    # initialize internal variables
    current_channel = ''
    key = ''
    output_list = []
    band_name_list = band_name_dict[category]

    # loop through list of bands and channels
    for line in catDict[category]:
        if quit_program:
            break

        # if line contains channel information
        if line[0].lower() == "channel":

            # set current channel
            current_channel = line[1].upper()

            # build key for given config file that will contain channel
            # possibly change this to build config file name
            key = config_key(line)

            # if config file entry does not already exist, initialize
            # entry for config file
            if key not in config_file_dict:
                # create empty dictionary for entry
                config_file_dict[key] = {}
                # enter frame type information for config file in new dictionary
                config_file_dict[key]["frame type"] = line[2]
                # enter observatory information for config file in new dictionary
                config_file_dict[key]["observatory"] = line[3]
                # create empty weighting function entry
                config_file_dict[key]["weighting functions string"] = ""
                # create empty band entry
                config_file_dict[key]["bands string"] = ""
                # create empty list of all band names that will be created by this
                # config file
                config_file_dict[key]["band names"] = []
                # create empty frame directory entry (to track where the frame
                # created by the blrms executable will be located
                config_file_dict[key]["frame directory"] = ""
                # create empty config file name entry
                config_file_dict[key]["config file path"] = ""
                # create empty dictionary to hold times and corresponding lists
                # of frame files to input into the blrms executable
                config_file_dict[key]["frame file lists"] = {}

            # for current channel, cycle through universal bands and
            # weighting functions
            for band in catDict['generalinformation']:
                if quit_program:
                    break

                # record and track band and weighting function
                config_file_dict, band_name_list = record_BAND(band, \
                                               'generalinformation', \
                                               current_channel, band_name_list, \
                                               config_file_dict, key)

            # for current channel, cycle through information shared by analysis
            # or secondary channels, depending on what kind of channel this is
            for band in catDict[category + 'Gen']:
                if quit_program:
                    break

                # record and track band and weighting function
                config_file_dict, band_name_list = record_BAND(band, \
                                               category + 'Gen', \
                                               current_channel, band_name_list, \
                                               config_file_dict, key)

        # for current channel, cycle through channel specific information
        elif "band" or "weight" in line[0].lower():

            # record and track band and weighting function
            config_file_dict, band_name_list = record_BAND(line, category, \
                                               current_channel, band_name_list, \
                                               config_file_dict, key)

        else:
            quit_program = True
            # better error message needed
            print("quiting because: ")
            print(line)

    return config_file_dict, quit_program, band_name_dict

# Helper function to make sure only one instance of entry exists
def single_instance_check(line, instance, quit_val):
    if instance:
        quit_val = True
        print("\nOnly one instance of '" + line[0].title() + "' should be \
in input file.\nPlease edit file appropriately.")
    return quit_val

# Helper funciton to find frames of specified type during specified time
def create_frame_file_list(frame_type, start_time, end_time, observatory, \
                           quit_program):
    # search for file location for a given frame type during specified times
    data_find = ['ligo_data_find','-s', start_time, '-e', end_time, '-o',
                 observatory, '--url-type', 'file', '--lal-cache', '--type',
                 frame_type]
    frame_locations_raw = subprocess.Popen(data_find, stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()#[0]
    if frame_locations_raw[1]:
        print(frame_locations_raw[1])
        print("The following shell command caused the above message:")
        #print(data_find)
        print(" ".join(data_find))
        quit_program = True
    frame_file_list = []
    all_found = False
    rest_of_loc = frame_locations_raw[0]
    while not all_found:
        if quit_program:
            break
        str_pos = rest_of_loc.find("localhost")
        if str_pos != -1:
            start_pos = rest_of_loc.find("localhost") + len("localhost")
            end_pos = rest_of_loc.find("\n")
            if end_pos != -1:
                frame_file_list.append(rest_of_loc[start_pos:end_pos])
                rest_of_loc = rest_of_loc[end_pos+1:]
            else:
                frame_file_list.append(rest_of_loc[start_pos:])
                all_found = True
        else:
            all_found = True
    # create frame list
    return frame_file_list, quit_program

# Helper function to create a file from the list of frame file locations
def create_list_file(frame_list,file_name,quit_program):
    # create string to write to file
    output_string = ""
    # create list to hold list of files in archive directory
    archived = []
    # for loop to go through list
    for line in frame_list:
        if quit_program:
            break
        output_string += line + "\n"
        if "archive" in line:
            archived.append(line)

    # check data for possibly archived data
    if archived:
        continue_program = ask_yes_no("Some frame files during the time \
specified may have to be loaded from tape.  Continue programe? ('y' or 'n'): ")
        display_files = ask_yes_no("Display frame files in 'archive' \
directory? ('y' or 'n'): ")

        if continue_program == 'n':
            quit_program = True

        if display_files == 'y':
            for line in archived:
                print(line)
            version_input("\n\nPress 'enter' to end program. ")

    # create file
    if not quit_program:
        file = open(file_name, 'w')
        file.write(output_string)
        file.close()
    return quit_program

# create job
def create_dag_job(job_number, condor_sub_loc, vars_entries, arg_list, filePointer, retry = 2):
    # create job entry strings
    jobNum = str(job_number)
    varEntries = [jobNum] + vars_entries
    argList = ["jobNumber"] + arg_list
    job_line = "JOB " + jobNum + " " + condor_sub_loc + "\n"
    retry_line = "RETRY " + jobNum + " " + str(int(retry)) + "\n"
    vars_line = "VARS " + jobNum
    for num in range(len(varEntries)):
        vars_line += ' ' + argList[num] + '="' + varEntries[num] + '"'
    vars_line += '\n\n'

    # enter string in dag file
    job_string = job_line + retry_line + vars_line
    filePointer.write(job_string)
    # iterate job number
    job_number +=1

    return job_number

# Helper function to write dag job entry
def blrms_dag_job(job_number, condor_sub_loc, frame_list_dict, frame_list,
                  conf_path, output_dir, file, test_interval = None):
    # possible arguments
    argList = ["confFile", "startTime", "endTime", "inlistFile", "outputDir"]
    # create variable entry list
    start = frame_list_dict[frame_list][0]
    if not test_interval:
        end = frame_list_dict[frame_list][1]
    else:
        end = str(int(start) + test_interval)
    list_file_name = frame_list
    vars_entries = [conf_path, start, end, frame_list, output_dir]

    # create job entry
    job_number = create_dag_job(job_number, condor_sub_loc, vars_entries,
                                argList, file)

    return job_number

# Helper function to write list of dag job entries
def write_blrms_jobs(job_number, conf_file_dict, job_tracker, condor_sub_loc,
                     file, output_dir = None):
    if output_dir:
        test_job = True
    else:
        test_job = False
    start_job = job_number
    for base_conf_name in conf_file_dict:
        conf_path = conf_file_dict[base_conf_name]["config file path"]
        # create handle to handle frame files
        frame_list_dict = conf_file_dict[base_conf_name]["frame file lists"]
        # determine if test jobs
        if not test_job:
            # get output directory
            output_dir = conf_file_dict[base_conf_name]["frame directory"]
            # run through frame lists
            for frame_list in frame_list_dict:
                job_number = blrms_dag_job(job_number, condor_sub_loc,
                                           frame_list_dict, frame_list,
                                           conf_path, output_dir, file)
        else:
            # grab any entry for test job
            frame_list = list(frame_list_dict.keys())[0]
            # determine test job length
            test_interval = 60
            # write test job
            job_number = blrms_dag_job(job_number, condor_sub_loc,
                                       frame_list_dict, frame_list, conf_path,
                                       output_dir, file, test_interval)
    end_job = job_number - 1
    # record range of job numbers just written
    job_tracker += [[start_job, end_job]]
    return job_tracker, job_number

# Helper function to write list of dag job entries
def write_analysis_dag_jobs(job_number, job_tracker, condor_sub_loc, filePointer, band_name_dict, arg_list):
#def write_analysis_dag_job(job_number, pickled_conf_name, job_tracker,
#                       condor_sub_loc, file, output_dir):
    start_job = job_number
    #arg_list = ["analysisChannel","secondaryChannel"]

    anChanList = band_name_dict["analysischannels"]
    secChanList = band_name_dict["secondarychannels"]
    print('len(anChanList) = ' + str(len(anChanList)))
    print('len(secChanList) = ' + str(len(secChanList)))
    for i in range(len(anChanList)):
        for n in range(len(secChanList)):
            vars_entries = [str(i), str(n)]
            job_number = create_dag_job(job_number, condor_sub_loc, vars_entries,
                                arg_list, filePointer)

    end_job = job_number - 1
    # record range of job numbers just written
    job_tracker += [[start_job, end_job]]
    return job_tracker, job_number

def write_displayPage_dag_job(job_number, job_tracker, condor_sub_loc, file):
#def write_analysis_dag_job(job_number, pickled_conf_name, job_tracker,
#                       condor_sub_loc, file, output_dir):
    start_job = job_number
    arg_list = []#["pickledFile", "outputDir"]
    vars_entries = []#pickled_conf_name, output_dir]

    job_number = create_dag_job(job_number, condor_sub_loc, vars_entries,
                                arg_list, file)

    end_job = job_number - 1
    # record range of job numbers just written
    job_tracker += [[start_job, end_job]]
    return job_tracker, job_number

# Helper function to write condor sub file
def write_sub_file(filename_base, executable, support_dir, args):
    condor_sub_filename = support_dir + "/" + filename_base + ".sub"
    condor_sub_file = open(condor_sub_filename, "w")
    # Separate strings
    universe = "universe = vanilla"
    executable_line = "executable = " + executable
    log = "log = " + support_dir + "/dagLogs/" + filename_base + ".log"
    error = "error = " + support_dir + "/jobLogs/" + filename_base + "$(jobNumber).err"
    output = "output = " + support_dir + "/jobLogs/" + filename_base + "$(jobNumber).out"
    arguments = 'arguments ='
    if args:
        arguments += ' " ' + args + ' "'
    notifications = "notification = error"
    queue = "queue 1"
    string_list = [universe, executable_line, log, error, output, arguments,
                   notifications, queue]
    output_string = "\n".join(line for line in string_list)
    condor_sub_file.write(output_string)
    condor_sub_file.close()
    return condor_sub_filename

# Helper function to enter job hierarchy in dagfile
def job_heirarchy(job_tracker, dagfile):
    #list_of_orderings = []
    output_string = ""
    for num in range(len(job_tracker) - 1):
        pair_1 = job_tracker[num]
        pair_2 = job_tracker[num+1]
        parent = "PARENT "
        child = "CHILD"
        for num in range(pair_1[0], pair_1[1]+1):
            parent += str(num) + " "
        for num in range(pair_2[0], pair_2[1]+1):
            child += " " + str(num)
        #order_string = parent + child
        output_string += parent + child + "\n\n"
        #list_of_orderings.append(order_string)
    #output_string = "\n\n".join(string for string in list_of_orderings)
    dagfile.write(output_string)

# Helper function to handle text input
def get_input_from_file_V_C_line(quit_program, filename):
    # Set function internal variables
    category = ''
    categories = ['generalinformation','analysischannels','secondarychannels']
    genInfo = []
    anChanGen = []
    secChanGen = []
    anChan = []
    secChan = []
    catDictionary = {'generalinformation':genInfo, 'analysischannels':anChan,
                     'secondarychannels':secChan, 'analysischannelsGen':anChanGen,
                     'secondarychannelsGen':secChanGen}
    if not quit_program:
        # Read file
        raw_data = read_text_file(filename, ' ')
        # Check if first line with information is a category defining line.
        # If not, send error message to user and quit.
        first_line = concat_line(raw_data[0])
        for cat in categories:
            if first_line == cat:
                category = cat
        if category == '':
            quit_program = True
            print("""
Please list a category (General Information, Analysis Channels, or
Secondary Channels) as the top line in file.

If the first line in file is such a line, please check that file is
plain text.  Otherwise, please contact this program's administrator.
""")
        # Put input into dictionary
        if not quit_program:
            catDictionary = sort_list(raw_data, categories, catDictionary)
    return catDictionary, quit_program

# Helper function to load text file with start and end times
def load_time_file_V_C_line(quit_val, lower_limit, upper_limit, filename,
                            endbuffer = None):
    return_data = []
    if not quit_val:
        # Read file
        raw_data = read_text_file(filename, ' ')
        # Check that file includes only start and end times
        low = lower_limit
        for line in raw_data:
            if len(line) == 2:
                if not number_check(line[0]) or not number_check(line[1]):
                    quit_val = True
                    print("Not all entries in " + filename + " are numbers.")
                    break
                if float(line[0]) < lower_limit:
                    quit_val = True
                    print(line[0] + " in " + filename + " is before the allowed\
 earliest GPS time of " + lower_limit + ".")
                    break
                elif float(line[1]) > upper_limit:
                    quit_val = True
                    print(line[1] + " in " + filename + " is after the allowed \
latest GPS time of " + upper_limit + ".")
                    break
                elif float(line[0]) < low:
                    quit_val = True
                    print(line[0] + " overlaps with or occurs before an earlier\
 time segment")
                    break
                elif float(line[1]) < float(line[0]):
                    quit_val = True
                    print("Error: End of segment occurs before beginning \
of segment: " + line[0] + " - " + line[1])
                    break
                low = float(line[1])
                if endbuffer:
                    if int(line[1]) - int(line[0]) > int(endbuffer):
                        return_data.append([line[0],str(int(line[1])-int(endbuffer))])
                else:
                    return_data.append(line)
            else:
                quit_val = True
                print("Not all lines in " + filename + " have 2 entries (a \
start time and an end time).")
                break
    return return_data, quit_val

# Helper function to handle input and write to config file
def create_config_file(catDictionary, conf_dir, quit_program):

    # Set function intern1al variables
    stride_line = []
    interval_line = []
    genInfo = catDictionary['generalinformation']
    base_output = ""
    conf_names = {}

    # Look through data in entries in category dictionary and check
    # to make sure each entry makes sense.
    for category in catDictionary:
        if quit_program:
            break
        for line in catDictionary[category]:
            if quit_program:
                break
            quit_program = check_line(line, category)

    # Order gen info so stride shows up first, followed by interval
    for line in genInfo:
        if line[0].lower() == 'stride':
            # first check if another instance exists
            quit_program = single_instance_check(line, stride_line, quit_program)
            # record instance
            stride_line = line

        elif line[0] == 'interval':
            # first check if another instance exists
            quit_program = single_instance_check(line, interval_line, quit_program)
            interval_line = line

    # move interval line to front of General Information
    if interval_line:
        genInfo.remove(interval_line)
        genInfo.insert(0, interval_line)

    # move stride line to front of General Information
    if stride_line:
        genInfo.remove(stride_line)
        genInfo.insert(0, stride_line)

    # grab stride and interval from General Information and put in base output text
    for line in genInfo:
        if quit_program:
            break
        base_output = write_config_string(line, \
                                    'generalinformation', base_output)

    # Calculate min frequency and filter bands accordingly
    minFreq = 1/float(stride_line[1])
    for category in catDictionary:
        catDictionary[category] = [line for line in catDictionary[category] if not len(line) == 3 or load_number(float(line[2])) >= minFreq]
        min_shift_list = [index for index, line in enumerate(catDictionary[category]) if len(line) == 3 and load_number(float(line[2])) <= minFreq]

        if min_shift_list:
            for index in min_shift_list:
                catDictionary[category][index][1] = str(minFreq)
            print("Some bands in " + category + " have been modified such that \
            their lower bound is " + str(minFreq) + " Hz.")

    # grab bands from all categories and put in config file dictionary
    # Write bands from
    conf_names, quit_program, band_name_dict = \
                   build_bands('analysischannels', catDictionary, quit_program, \
                               conf_names)
    conf_names, quit_program, band_name_dict = \
                   build_bands('secondarychannels', catDictionary, quit_program, \
                               conf_names, band_name_dict)

    if not quit_program:

    # cycle through configuration file entries
        for key in conf_names:

            # create output file names
            # add date to this at some point
            filename = conf_dir + "/Config" + "_" + key + ".conf"
            outString = base_output + conf_names[key]["weighting functions string"] + \
                        conf_names[key]["bands string"]

            # write to file
            outfile = open(filename, "w")
            outfile.write(outString)
            outfile.close()
            #confPath = conf_dir + "/" + filename
            conf_names[key]["config file path"] = filename#confPath

    return conf_names, quit_program, band_name_dict

# Create frame file lists
def create_frame_file_list_files(conf_names, quit_program, times, list_dir):
    # for each config file find corresponding frame files
    for conf in conf_names:
        if quit_program:
            break
        # create handle to dictionary of file paths for config file
        for time in times:
            if quit_program:
                break
            # create handles to frame information and times
            frame_type = conf_names[conf]["frame type"]
            observatory = conf_names[conf]["observatory"]
            start = time[0]
            end = time[1]
            # create name of file that will list frame files
            list_file_name = list_dir + "/" + observatory + "-" + frame_type + \
                             "-" + time[0] + "_" + time[1] + ".list"
            # find frames
            frames, quit_program = create_frame_file_list(frame_type,start,end, observatory, quit_program)
            # write to file
            quit_program = create_list_file(frames, list_file_name, quit_program)
            # record path to list file
            #list_file_path = list_dir + "/" + list_file_name
            conf_names[conf]["frame file lists"][list_file_name] = time
    return conf_names, quit_program

# Create sub-directories to hold frame file output by frame type
def create_frame_output_directories(conf_file_dict, base_frame_dir, quit_val):
    if not quit_val:
        # loop through entries for each input frame file type and observatory
        for dirName in conf_file_dict:
            # create new directory
            dirPath = create_dir(base_frame_dir + "/" + dirName, False)
            # keep track of directory
            conf_file_dict[dirName]["frame directory"] = dirPath

# Pickle config file dictionary (it tracks where everything is)
def pickle_object(config_file_dict, support_dir, pickled_file_name,
                  quit_program):
    if not quit_program:
        # open file in byte write mode
        filename = support_dir + "/" + pickled_file_name + ".txt"
        file = open(filename, "wb")
        pickle.dump(config_file_dict, file)
        file.close()
        return filename

# Create shell script to allow dag to execute python analysis code
def create_shell_scripts(pickled_conf_dict, pickled_band_names, pickled_times,
                        support_dir, plot_dir, shell_path, python_script,
                        input_struct, quit_program):
    if not quit_program:
        bandFile = open(pickled_band_names, "rb")
        band_name_dict = pickle.load(bandFile)
        anChanList = band_name_dict["analysischannels"]
        secChanList = band_name_dict["secondarychannels"]
        for i in range(len(anChanList)):
            for n in range(len(secChanList)):     
                filename = support_dir + "/blrmsPyAnalysis_" + str(i) + "_" + str(n) + ".sh"
                file = open(filename, "w")
                # find bash path - can automate this later, maybe manual for now
                file.write(shell_path + "\n\n")
                file.write("python " + python_script + " -f " + pickled_conf_dict +
                           " -b " + pickled_band_names + " -t " + pickled_times + " -d "
                           + plot_dir + " -q " + str(i) + " -w " + str(n))
                if input_struct.lowData:
                    file.write(" -l " + input_struct.lowData)
                if input_struct.highData:
                    file.write(" -u " + input_struct.highData)
                if input_struct.xMaxGen:
                    file.write(" -x " + input_struct.xMaxGen)
                if input_struct.yMaxGen:
                    file.write(" -y " + input_struct.yMaxGen)
                if input_struct.timeIntervals:
                    file.write(" -i " + input_struct.timeIntervals)
                file.close()
        filename = support_dir + "/makeDisplayPage.sh"
        file = open(filename, "w")
        # find bash path - can automate this later, maybe manual for now
        file.write(shell_path + "\n\n")
        file.write("python " + "/home/vincent.roma/public_html/upconversionCode/prototypeNCAT6/makeDisplayPage.py -i " + str(len(anChanList)) + " -n " + str(len(secChanList)) + " -d " + plot_dir)
        file.close()
        os.chmod(filename, 0764)
        
        filename = support_dir + "/blrmsPyAnalysis"
        return filename

# Create shell script to allow dag to execute python analysis code
def create_shell_script(pickled_conf_dict, pickled_band_names, pickled_times,
                        support_dir, plot_dir, shell_path, python_script,
                        input_struct, quit_program):
    if not quit_program:
        bandFile = open(pickled_band_names, "rb")
        band_name_dict = pickle.load(bandFile)
        anChanList = band_name_dict["analysischannels"]
        secChanList = band_name_dict["secondarychannels"]
        filename = support_dir + "/blrmsPyAnalysis.sh"
        with open(filename, "w") as shellFile:
            # find bash path - can automate this later, maybe manual for now
            shellFile.write(shell_path + "\n\n")
            shellFile.write("python " + python_script + " -f " + pickled_conf_dict +
                           " -b " + pickled_band_names + " -t " + pickled_times + " -d "
                           + plot_dir + " -q $1 -w $2")
            if input_struct.lowData:
                shellFile.write(" -l " + input_struct.lowData)
            if input_struct.highData:
                shellFile.write(" -u " + input_struct.highData)
            if input_struct.xMaxGen:
                shellFile.write(" -x " + input_struct.xMaxGen)
            if input_struct.yMaxGen:
                shellFile.write(" -y " + input_struct.yMaxGen)
            if input_struct.timeIntervals:
                shellFile.write(" -i " + input_struct.timeIntervals)
        filename = support_dir + "/makeDisplayPage.sh"
        file = open(filename, "w")
        # find bash path - can automate this later, maybe manual for now
        file.write(shell_path + "\n\n")#shellpath
        file.write("python " + "/home/vincent.roma/public_html/upconversionCode/prototypeNCAT6/makeDisplayPage.py -i " + str(len(anChanList)) + " -n " + str(len(secChanList)) + " -d " + plot_dir)
        file.close()
        filename2 = filename

        filename = support_dir + "/blrmsPyAnalysis.sh"
        return [filename, filename2]

# create dag submission files
def create_dag(config_file_dict, pickled_conf_name, support_dir, plot_dir,
               blrmsExc, shellScript, pickled_bands, quit_program):
    if not quit_program:
        # create condor blrms executable submit file
        args = "-conf $(confFile) -start $(startTime) -end $(endTime) -inlist \
$(inlistFile) -dir $(outputDir)"
        blrms_sub_filename = write_sub_file("blrmsExc", blrmsExc,
                                            support_dir, args)
        # create condor shell script that executes python plot analysis submit file
        #shell_arg = "-file $(pickledFile) -dir $(outputDir)"
        shellArguments = ["analysisChannel", "secondaryChannel"]
        shell_arg = " " + " ".join("$(" + x + ")" for x in shellArguments) + " "
        python_arg = None

        # load band tracker file
        bandFile = open(pickled_bands, "rb")
        band_name_dict = pickle.load(bandFile)
        bandFile.close()
        anChanList = band_name_dict["analysischannels"]
        secChanList = band_name_dict["secondarychannels"]
        #for i in range(len(anChanList)):
         #   for n in range(len(secChanList)):
          #      tempShellScript = shellScript + "_" + str(i) + "_" + str(n)+".sh"
                #write_sub_file("blrmsAnalysis_" + str(i) + "_" + str(n), tempShellScript,
                 #                            support_dir, shell_arg)
        write_sub_file("blrmsAnalysis", shellScript, support_dir, shell_arg)
        python_sub_filename = support_dir + "/blrmsAnalysis.sub"
        #python_makePage_filename = write_sub_file("makeDisplayPage", "/home/vincent.roma/public_html/upconversionCode/prototypeNCAT6/makeDisplayPage.py", support_dir, python_arg)
        python_makePage_filename = write_sub_file("makeDisplayPage", support_dir + "/makeDisplayPage.sh", support_dir, python_arg)
        # create dag file
        filename = support_dir + "/blrmsJobs.dag"
        dagfile = open(filename, "w")
        job_number = 0
        job_tracker = []
            # create test blrms jobs
        job_tracker, job_number = write_blrms_jobs(job_number, config_file_dict,
                                       job_tracker, blrms_sub_filename, dagfile,
                                       support_dir)
            # create rest of jobs
        job_tracker, job_number = write_blrms_jobs(job_number, config_file_dict,
                                       job_tracker, blrms_sub_filename, dagfile)
            # create shell executable with pickled object input
        job_tracker, job_number = write_analysis_dag_jobs(job_number,
                                                         job_tracker,
                                                         python_sub_filename,
                                                         dagfile, band_name_dict, shellArguments)
        job_tracker, job_number = write_displayPage_dag_job(job_number, job_tracker,
                                                            python_makePage_filename,
                                                            dagfile)

#        job_tracker, job_number = write_analysis_dag_job(job_number,
#                                                         pickled_conf_name,
#                                                         job_tracker,
#                                                         python_sub_filename,
#                                                         dagfile, plot_dir)
            # create job hierarchy
        job_heirarchy(job_tracker, dagfile)

        dagfile.close()
        return filename

def run_dag(quit_val, dagfile, maxJobs = None):
    if not quit_val:
        dagCommand = ['condor_submit_dag']
        if maxJobs:
            dagCommand += ["-maxjobs",maxJobs]
        dagCommand += [dagfile]
        condor_output = subprocess.Popen(dagCommand, stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()#[0]
        if condor_output[1]:
            print(condor_output[1])
            print("The following shell command caused the above message:")
            print(dagCommand)
            #quit_program = True
    #condor_worked = frame_locations_raw[0]
