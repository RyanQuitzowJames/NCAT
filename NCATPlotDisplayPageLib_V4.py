import os

def tabify(string, indent):
    tabs = ""
    for num in range(indent):
        tabs += "    "
    output = tabs + string
    return output

def glueFileLocation(directory, filename):
    output = None
    if directory[-1] == "/":
        if filename[0] == "/":
            output = directory + filename[1:]
        else:
            output = directory + filename
    else:
        if filename[0] == "/":
            output = directory + filename
        else:
            output = directory + "/" + filename
    return output

def wrap(string, element, indent = 0):
    tabs = tabify("", indent)
#    for num in range(indent):
#        tabs += "    "
#    output = tabs + "<" + element + ">" + "\n" + string + "\n" + tabs + "</" \
#             + element + ">"
#    print(type(element))
#    print(type(tabs))
#    print(type(string))
    output =  tabs + "<" + element + ">" + "\n" + string + "\n" + tabs + "</" \
             + element + ">"
#    output = tabify(output, indent)
    return output

def wrap_v2(string, element, indent = 0):#, newlines = 3):
    tabs = tabify("", indent)
#    output = attachNewLines(string, newlines)
#    output =  tabs + "<" + element + ">" + "\n" + output + "\n" + tabs + "</" \
    output =  tabs + "<" + element + ">" + string + "</" \
             + element + ">"
    return output

def wrap_v3(string, element, indent = 0, modifiers = None):
    tabs = tabify("", indent)
    if modifiers:
        output =  tabs + "<" + element + " " + modifiers + ">" + string + "</" \
             + element + ">"
    else:
        output =  tabs + "<" + element + ">" + string + "</" \
             + element + ">"
    return output

def wrapBody(string):
    output = wrap(string,"body",1)
    #output = "\t<body>\n"+string+"\n\t</body>"
    return output

def wrapHead(string):
    output = wrap(string,"head",1)
    #output = "\t<body>\n"+string+"\n\t</body>"
    return output

def wrapHTML(string):
    output = "<!DOCTYPE html>\n" + wrap(string,"html")
    return output

def attachNewlines(string, newline):
    output = string
    if newline == 1:
        output = output + "\n"
    elif newline == 2:
        output = "\n" + output
    elif newline == 3:
        output = "\n" + output + "\n"
    return output

def wrapImage(imageurl, indent = 0, newline = 0, height = 200, width = 200):
#    tabs = ""
#    for num in range(indent):
#        tabs += "    "
#    output = tabs + '<img src="' + imageurl + '"'
    output = '<img src="' + imageurl + '"'
    if height:
        output += ' height="' + str(height) + '"'
    if width:
        output += ' width="' + str(width) + '"'
    output += ' />'
    output = tabify(output, indent)
    output = attachNewlines(output, newline)
#    if newline == 1:
#        output = output + "\n"
#    elif newline == 2:
#        output = "\n" + output
#    elif newline == 3:
#        output = "\n" + output + "\n"
    return output

def wrapLinkedImage(imageurl, indent = 0, newline = 0, height = 200, width = 200):
    output = wrapImage(imageurl, height = height, width = width)
    output = '<a href="' + imageurl + '">' + output + '</a>'
    output = tabify(output, indent)
    output = attachNewlines(output, newline)
    return output

def titleSection(title, string, indent):
    titleText = attachNewlines(wrap(tabify(title,indent+1),"p",indent),1)
    output = titleText + string
    return output

def smallSection(directory, fileList, indent = 2, newline = 0, title = None):
    fileLocations = [glueFileLocation(directory, file) for file in fileList]
    htmlTextList = [wrapLinkedImage(file) for file in fileLocations]
    htmlText = "".join(file for file in htmlTextList)
    output = tabify(htmlText, indent + 1)
    output = wrap(output,"p",indent)
    if title:
        output = titleSection(title, output, indent)
#        titleText = attachNewlines(wrap(tabify(title,indent+1),"p",indent),1)
#        output = titleText + output
    output = attachNewlines(output, newline)
    return output

def largeSection(directory, fileList, title = None, newline = 0, indent = 2, subDirList = None):
    if not subDirList:
        subDirList = os.listdir(directory)
    dirList = [glueFileLocation(directory, subDir) for subDir in subDirList]
    htmlTextList = [smallSection(dirList[num], fileList, indent, 2, subDirList[num]) for num in range(len(dirList))]
    htmlTextList[0] = htmlTextList[0][htmlTextList[0].find("\n")+1:]
    output = "".join(line for line in htmlTextList)
    if title:
        output = titleSection(title, output, indent)
#        titleText = attachNewlines(wrap(tabify(title,indent+1),"p",indent),1)
#        output = titleText + output
    output = attachNewlines(output, newline)
    return output

def makeSpecificPageString(directoryList, fileList, strideDir):
    focusFileList = [glueFileLocation(strideDir, file) for file in fileList]
    pageText = largeSection(directoryList[0],focusFileList, "Before Oct 14", 2)
    afterSection = largeSection(directoryList[1],focusFileList, "CR", 2)
    afterSection += largeSection(directoryList[2],focusFileList, "SB", 2)
    pageText += titleSection("After Oct 14", afterSection, 2)
    pageText = wrapBody(pageText)
    pageText = wrapHTML(pageText)
    return pageText

def makeSpecificPage(directoryList, fileList, strideDir, filename):#, saveDir = None):
    pageText = makeSpecificPageString(directoryList, fileList, strideDir)
    #if saveDir:
        
    with open(filename,"w") as outfile:
        outfile.write(pageText)

def wrap_table(string, indent):
    tabs = tabify("", indent)
    output = tabs +'<table id="theTable" border="1" style="border-collapse: collapse;">\n' \
             + string + "\n" + tabs + '</table>'
    return output

def make_row(listEl, indent = 0, baseString = None):
    if baseString:
        newList = [glueFileLocation(baseString, el) for el in listEl]
    else:
        newList = listEl
    output = "\n".join(tabify(wrap_v2(x,'td'),indent+1) for x in newList)
#    output = "".join(x for x in output)
    output = wrap(output,'tr',indent)
    return output

def make_row_v2(string,indent):
    output = wrap(string,'tr',indent)
    return output

def make_data_cells(listEl, indent = 0):#, baseString = None):
#    if baseString:
#        newList = [glueFileLocation(baseString, el) for el in listEl]
#    else:
#        newList = listEl
#    newList = listEl
#    output = "\n".join(tabify(wrap_v2(x,'td'),indent) for x in newList)
    output = "\n".join(tabify(wrap_v2(x,'td'),indent) for x in listEl)
#    output = "".join(x for x in output)
#    output = wrap(output,'tr',indent)
    return output

def make_tall_data_cell(string, indent = 0, rows = 2):
    output = tabify(wrap_v3(string,'td', modifiers = 'rowspan="'+ str(rows) + '"'), indent)
    return output

def small_table_section(rowTitle, directory, fileList, indent = 3, newline = 0, title = None):
    fileLocations = [glueFileLocation(directory, file) for file in fileList]
    htmlTextList = [rowTitle] + [wrapLinkedImage(file) for file in fileLocations]
    #htmlText = "".join(file for file in htmlTextList)
    output = make_row(htmlTextList,indent)
    #output = tabify(htmlText, indent + 1)
    #output = wrap(output,"p",indent)
    #if title:
#        output = titleSection(title, output, indent)
#        titleText = attachNewlines(wrap(tabify(title,indent+1),"p",indent),1)
#        output = titleText + output
    output = attachNewlines(output, newline)
    return output

def small_table_section_v2(rowTitle, directory, fileList, plotTypeDict, indent = 3, newline = 0, title = None):

    # make title cell for double (or more) row
    output = make_tall_data_cell(rowTitle, indent+1, len(plotTypeDict)) + "\n"

    first_row = True

    # iterate through and create rows
    for key in plotTypeDict:
        # create path string
        path = glueFileLocation(directory, plotTypeDict[key])
        fileLocations = [wrapLinkedImage(glueFileLocation(path, file)) for file in fileList]
        # create row data
        temp_output = make_data_cells(fileLocations,indent+1)
        # add single row title to row
        temp_output = make_data_cells([key], indent+1) + "\n" + temp_output
        # handle first row exception (title of double or more rows)
        if first_row:
            output = make_row_v2(output + temp_output, indent) + "\n"
            first_row = False
        else:
            output += make_row_v2(temp_output, indent) + "\n"
    # slice string to cut off trailing new line character
    output = output[:-1]
    output = attachNewlines(output, newline)

    return output

def small_table_section_v3(directory, fileList, plotTypeList, plotTypeDict, indent = 3, newline = 0, title = None):

    # make title cell for double (or more) row
#    output = make_tall_data_cell(rowTitle, indent+1, len(plotTypeDict)) + "\n"

#    first_row = True

# each row will have every type of plot. that means the top row should be a single
# blank followed by a cell for each plot type. currently two plot typee displayed.

#    titleCells = ["", "Channel vs channel plots", "Time vs rms^2 plots"]
#    titleCells = [""] + [x for x in plotTypeDict.keys()]
    titleCells = [""] + [x for x in plotTypeList]

    output = make_data_cells(titleCells, indent+1) + "\n"
    output = make_row_v2(output, indent) + "\n"

    # iterate through and create rows
    ranksLength = len(fileList)
    print('ranksLength = ' + str(ranksLength))
    i = 0
    while (i < 15 and i < ranksLength):
        name = fileList[i]
        name = name[::-1]
        name = name[:name.find("/")]
        name = name[::-1]
        # get row name
        temp_output = [name[:name.index("vs")] + "<br/>and<br/>" + name[name.index("vs") + 2:name.index(".")]]
        # iterate through plot types
        for key in plotTypeList:
            path = glueFileLocation(directory, plotTypeDict[key])
            fileLocation = wrapLinkedImage(glueFileLocation(path, name))
            print('fileLocation = ' + fileLocation)
            temp_output.append(fileLocation)
#        fileLocations = [wrapLinkedImage(glueFileLocation(path, file)) for file in fileList]
        # create row data
        temp_output = make_data_cells(temp_output,indent+1)
        output += make_row_v2(temp_output, indent) + "\n"
        i = i + 1

    # slice string to cut off trailing new line character
    output = output[:-1]
    output = attachNewlines(output, newline)

    # make table
    output = wrap_table(output, indent - 1)
    # Make display buttons
    tabs = tabify("", indent)
    buttonCode = '\n' + tabs + '<button onclick="buttonGo(15)" type="button">Display 15 More Rows</button>\n' + tabs + '<button onclick="buttonGo(50)" type="button">Display 50 More Rows</button>\n'
    fileList_commas = ','.join(fileList)
    
    scriptCode = tabs + '<script> \n' + tabs + 'var i = ' + str(i) + ";\n" + tabs + "var ranks_commas = '" + fileList_commas + "'; \n" + tabs + 'var ranks_list = ranks_commas.split(",");\n' + tabs + 'function buttonGo(x) {\n' + tabs + tabs + "var table = document.getElementById('theTable');\n" + tabs + tabs + "var currentLimit = i + x;\n" + tabs + tabs + "while (i < ranks_list.length && i < currentLimit) {\n" + tabs + tabs + tabs + "var fullName = ranks_list[i];\n" + tabs + tabs + tabs + 'var cutIndex = fullName.search("dPlots");\n' + tabs + tabs + tabs + 'fullName = fullName.slice(cutIndex + 7, -1) + "g";\n' + tabs + tabs + tabs + "var name = fullName;\n" + tabs + tabs + tabs + "var fullNameLength = fullName.length;\n" + tabs + tabs + tabs + "name = name.slice(0, fullNameLength - 3);\n" + tabs + tabs + tabs + 'name = name.replace(" vs ", "<br>and<br>");\n' + tabs + tabs + tabs + "var tr = document.createElement('tr');\n" + tabs + tabs + tabs + "var td0 = tr.appendChild(document.createElement('td'));\n" + tabs + tabs + tabs + "var td1 = tr.appendChild(document.createElement('td'));\n" + tabs + tabs + tabs + "var td2 = tr.appendChild(document.createElement('td'));\n" + tabs + tabs + tabs + "td0.innerHTML = name;\n" + tabs + tabs + tabs + '''td1.innerHTML = '<a href="plots/minimalLinearBoundPlots/' + fullName + '"><img src="plots/minimalLinearBoundPlots/' + fullName + '" height="200" width="200" /></a>';\n ''' + tabs + tabs + tabs + '''td2.innerHTML = '<a href="plots/timeSqrdPlots/' + fullName + '"><img src="plots/timeSqrdPlots/' + fullName + '" height="200" width="200" /></a>';\n ''' + tabs + tabs + tabs + "table.appendChild(tr);\n" + tabs + tabs + tabs + 'i = i + 1;\n' + tabs + tabs + '}\n' + tabs + '}\n' + tabs + '</script>\n'
    output = output + buttonCode + scriptCode

    return output

def make_column_titles(fileList, indent = 3, newline = 0, title = None):
    
    # get strings for title row
    firstLine = ["",""] + [name[:name.index("vs")] + "<br/>and<br/>" + name[name.index("vs") + 2:name.index(".")] for name in fileList]
    firstLine = make_data_cells(firstLine, indent + 1)
    output = [make_row_v2(firstLine, indent)]
#    print(output)

    return output

def large_section_v2(directory, fileList, plotTypeDict, title = None, newline = 0, indent = 2, subDirList = None):
    if not subDirList:
        subDirList = os.listdir(directory)
    dirList = [glueFileLocation(directory, subDir) for subDir in subDirList]
#    firstLine = ["",""] + [name[:name.index("vs")] + "and" + name[name.index("vs") + 2:name.index(".")] for name in fileList]
    firstLine = make_column_titles(fileList)
#    print(" ".join(x for x in firstLine))
#    print(directory)
#    print(fileList)
#    print(plotTypeDict)
    htmlTextList = [small_table_section_v2(subDirList[num] + " GPS", dirList[num], fileList, plotTypeDict, indent + 1, 2, subDirList[num]) for num in range(len(dirList))]
    htmlTextList[0] = htmlTextList[0][htmlTextList[0].find("\n")+1:]
    print(firstLine)
    print(htmlTextList)
    htmlTextList = firstLine + htmlTextList
    output = "".join(line for line in htmlTextList)
    output = wrap_table(output, indent)
    if title:
        output = titleSection(title, output, indent)
#        titleText = attachNewlines(wrap(tabify(title,indent+1),"p",indent),1)
#        output = titleText + output
    output = attachNewlines(output, newline)
#    print(output)
    return output

def makeSpecificPageString_tableV(directoryList, fileList, strideDir, plotTypeDict):
    focusDict = dict((key, glueFileLocation(strideDir, plotTypeDict[key])) for key in plotTypeDict)
    #print(focusDict)
 #   focusFileList = [glueFileLocation(strideDir, file) for file in fileList]
    pageText = large_section_v2(directoryList[0],fileList, focusDict, "Before Oct 14", 2)#, subDirList = ["dir1","dir2"])
    afterSection = large_section_v2(directoryList[1],fileList, focusDict, "CR", 2)#, subDirList = ["dir1","dir2"])
    afterSection += large_section_v2(directoryList[2],fileList, focusDict, "SB", 2)#, subDirList = ["dir1","dir2"])
    pageText += titleSection("After Oct 14", afterSection, 2)
    pageText = wrapBody(pageText)
    pageText = wrapHTML(pageText)
    return pageText

def makeSpecificPage_tableV(directoryList, fileList, strideDir, filename, plotTypeDict):#, saveDir = None):
    pageText = makeSpecificPageString_tableV(directoryList, fileList, strideDir, plotTypeDict)
    #if saveDir:
        
    with open(filename,"w") as outfile:
        outfile.write(pageText)

# Make sure to place this script in the base directory of the analysis ## No longer valid?

def make_display_page(directory):
    #ranks = []
    # load ranking text file
    with open(directory + "/plots/minimalLinearBoundPlots/rankings.txt","r") as infile:
        print('directory = ' + directory)
        #print('infile = ' + infile.read())
        '''inString = infile.read()
        print('inString = ' + inString)
        inList = inString.split('png')
        print('len(inList) = ' + str(len(inList)))
        inList = inList[:-1]
        inList = [x.strip() for x in inList]'''

#ranks = []
    ##with open("plots/minimalLinearBoundPlots_v2/rankings.txt","r") as infile:
        ranks = [x.split(",") for x in infile if len(x) > 0 and x != '\n']
        print(ranks)
    print('len(ranksInitial) = ' + str(len(ranks)))
    for x in ranks:
        print('x =')
        print(x)
    ranks = [[x[0], x[1].strip()] for x in ranks]
    ranks2 = [x[1] for x in ranks]
    #print('len(ranks2) = ' + str(len(ranks2)))
    #ranks2_text = ' '.join(ranks2)

    try:
        with open(directory + "/plots/timeSqrdPlots/subplots.txt","r") as infile:
            #subplots = [[x.split(","), ] for x in infile if len(x) > 0]
            subplots = [x.split(",") for x in infile if len(x) > 0]
            print(subplots) #debug
            subplots = [[x[0], x[1].strip()] for x in subplots]
    except:
        subplots = None 

    plotTypeList = ["Channel vs channel plots", "Time vs rms^2 plots"]

    plotTypeDict = {}
    
    ####
    # eventually change this to reference the current version of timeSqrdPlots
    # and other such directories instead of the absolute directory as currently
    # defined. Basically make this smart so it looks up the directory instead
    # of having the directory hard coded.
    ####
    
    #plotTypeDict[plotTypeList[0]] = "minimalLinearBoundPlots_v2"
    plotTypeDict[plotTypeList[0]] = "minimalLinearBoundPlots"
    #plotTypeDict[plotTypeList[1]] = "timeSqrdPlots_v2"
    plotTypeDict[plotTypeList[1]] = "timeSqrdPlots"

    if subplots:
        plotListExtended = [["Channel vs channel plots\n" + x[1], "Time vs rms^2 plots\n" + x[1]] for x in subplots]
        plotListExtendedFlat = [item for sublist in plotListExtended for item in sublist]
        plotTypeList += plotListExtendedFlat
        print(plotTypeList) #debug
        for index, val in enumerate(plotListExtended):
            currentDir = subplots[index][0][subplots[index][0].find("timeSqrdPlots"):]
            vsDir = currentDir.replace("timeSqrdPlots","minimalLinearBoundPlots")
    #        print(val)#debug
    #        print(subplots[index][0][subplots[index][0].find("timeSqrdPlots"):])#debug
            #plotTypeDict[plotTypeList[index + 2]] = currentDir
            plotTypeDict[val[0]] = vsDir
            plotTypeDict[val[1]] = currentDir
        print(plotTypeDict) #debug

    header = """        <style type="text/css">
            img
            {
                width:auto;
                height:200px
            }
            th, td
            {
                padding: 15px
            }
        </style>"""

    header = wrapHead(header)

    #body_text = small_table_section_v3("plots", ranks2, plotTypeList, plotTypeDict)
    #body_text = body_text + "/n" + script_text
    #page_text = wrapHTML(header + "/n" + wrapBody(body_text))
    page_text = wrapHTML(header + "\n" + wrapBody(small_table_section_v3("plots", ranks2, plotTypeList, plotTypeDict)))#, times_path)))

    with open(directory + "/pageDisplayTest.html","w") as outfile:
        outfile.write(page_text)
