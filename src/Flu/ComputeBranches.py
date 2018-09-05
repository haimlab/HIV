################################################
# Program to compute the different branches
# computed based on the distance matrix of 
# strains of flu (distance of <= .007)
# Note: The distance matrix is stored in a csv
# named "Distance Data" stored on the desktop.
################################################
#Author Kallin Khan, Department of Microbiology
#Updated: 4/3/17
################################################

import csv

def ComputeBranches(Reader, numberOfRows, offSet):
    num = numberOfRows - 1

    Matrix = [[0.0 for x in range(num)] for y in range(num)]
    Names = ['' for x in range(num)]
    Branch = [-1 for x in range(num)]

    count = 0
    rowNum = 0
    
    ## This loop reads each row of the csv, reads each element of each row,
    ## adds the name of each case to the array 'Names', and lastly
    ## adds each value into a matrix in Python.
    for row in Reader:                          
        colNum = 0
        t = (', '.join(row))
        mylist = t.split(",")

        for element in mylist:                  
            if colNum == 0:                    
                Names[rowNum] = element         
            elif colNum != 0:
                Matrix[rowNum][colNum-1] = float(element)
            colNum += 1
        rowNum += 1
    
        count += 1              ## These lines end the loop after all elements are read
        if count == num:        ## and ignore any text below the matrix.
            break


    ## This for loop assigns a branch number to each element as strains
    ## with a diffference of .007 or less will have the same in a branch.
    BranchNum = 0
    thresh = .007    
    for row in range(num-1,0,-1):  
        if Branch[row] == -1:      
            Branch[row] = BranchNum
            for column in range(0,num):
                if Matrix[row][column] <= thresh and Branch[column] == -1:
                    Branch[column] = BranchNum
            BranchNum += 1


    ## Prints out the branch number, strain number, and name of each strain.
    rowNumbers = [1]*(num)
    rowNum = offSet
    for index in range(0,BranchNum): 
        ind = 0
        for element in Branch: 
            if element == index:
                print(str(index) + " " + str(ind + 1) + " " + Names[ind])
                rowNumbers[ind] = rowNum
                rowNum += 1
            ind += 1
        rowNum += 2 ##add more here to create more spaces between branches
        print('  ')

    if rowNumbers[0] == 1:
        count = 0
        for o in rowNumbers:
            rowNumbers[count] += 1
            count += 1

    return(rowNumbers, BranchNum)
