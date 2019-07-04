import os.path

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
#    dirName = r'\\srv-echange\echange\Sylvain\2019-04-16 - T2594Al - 300K\4\HYP1-T2594Al-300K-Vacc5kV-spot7-zoom6000x-gr600-slit0-2-t5ms-cw440nm'
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            print(fullPath)
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
    return allFiles