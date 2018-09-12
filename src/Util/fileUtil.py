def isSameFileContent(fileNames):
    # guard bad input
    try:
        if len(fileNames < 2):
            raise Exception('Need at least two files to compare')
    except TypeError:
    files = []
    for fileName in fileNames:
        files.append(open(fileName, 'r'))
    content = files[1]
        
    