#!/usr/bin/python3

# an itty bitty script to find the differences between scaffolds

from os import listdir
import argparse

parser = argparse.ArgumentParser(description='find the difference between genetic scaffolds and save ones above certain thresholds')
parser.add_argument('--input', help='filename for the input file')
parser.add_argument('--output', help='filename for the output file', default="out")
parser.add_argument('--threshold', help='only extract values strictly greater than the threshold')
parser.add_argument('--moderate', help='define what is a moderate difference between scaffolds')
parser.add_argument('--severe', help='define what is a severe difference between scaffolds')
parser.add_argument('--debug', help='print debug output', action="store_true")
args = parser.parse_args()

aboveTInd = []
linkageGroupInd = []
moderateScaffs = []
severeScaffs = []
index1 = 0
index2 = 0
index3 = 0
modGaps = 0
sevGaps = 0
skip2Lines = 2
lgHasDiff = False

# open the file, find all indices above a certain threshold, and save them
with open(args.input) as f1:
    content = f1.readlines()
    tempList = content[skip2Lines].split()

    # store the index of every linkage group above threshold
    for i in tempList:
        if i.isnumeric():
            if int(i) > int(args.threshold):
                aboveTInd.append(index1)
            index1 += 1

    # find lines where scaffolds above threshold start and calculate the difference between adjacent scaffolds 
    for i in content:

        # end loop if there are no more linkage groups above threshold
        if (index3 >= len(aboveTInd)):
            break

        # find groups above threshold 
        strToComp = "group lg" + str(aboveTInd[index3])
        if strToComp in i:
            j = skip2Lines
            tempSev = []
            tempMod = []

            # iterate through linkage group and compare adjacent scaffolds
            while True:
                scaff1 = content[index2 + j].split()
                scaff2 = content[index2 + j + 1].split()

                # end loop once terminus is reached
                if scaff2[0] == ";ENDOFGROUP":
                    break
                diff = float(scaff2[1]) - float(scaff1[1])
                if diff > float(args.severe):
                    severe = "There is a severe difference between " + scaff1[0] + " and " + scaff2[0] + " and that is: " + str(diff) + "\n"
                    tempSev.append(severe)
                    lgHasDiff = True
                    sevGaps += 1
                elif diff > float(args.moderate):
                    moderate = "There is a moderate difference between " + scaff1[0] + " and " + scaff2[0] + " and that is: " + str(diff) + "\n"
                    tempMod.append(moderate)
                    lgHasDiff = True
                     modGaps += 1
                j += 1

            # add index of linkage group which had a difference
            if aboveTInd[index3] not in linkageGroupInd and lgHasDiff == True:
                lgHasDiff = False
                linkageGroupInd.append(aboveTInd[index3])
                severeScaffs.append(tempSev)
                moderateScaffs.append(tempMod)
            index3 += 1
        index2 += 1

# close the file
f1.close()

counter = 0

# write severeScaffs and moderateScaffs to log
with open(args.output, "x") as toWrite:

    # write and print header for file
    header = "Out of " + str(len(tempList)) + " linkage groups, " + str(len(aboveTInd)) + " had more than " + str(args.threshold) + " scaffolds, out of those " + str(len(aboveTInd)) + " groups, " + str(len(linkageGroupInd)) + " had moderate to severe gaps. There are a total of " + str(sevGaps) + " severe gaps, and " + str(modGaps) + " moderate gaps across all linkage groups\n\n"
    toWrite.write(header)
    print(header + "See: " + str(args.output) + " for a more detailed log.")
    for i in linkageGroupInd:

        # write and print linkage group number and severe differences
        if len(severeScaffs[counter]) != 0:
            sevIntro = "The severe differences in linkage group: " + str(i) + " are:\n\n"
            toWrite.write(sevIntro)
            toWrite.writelines(severeScaffs[counter])
            toWrite.write("\n")

        # write and print linkage group number and moderate differences
        if len(moderateScaffs[counter]) != 0:
            modIntro = "The moderate differences in linkage group: " + str(i) + " are:\n\n"
            toWrite.write(modIntro)
            toWrite.writelines(moderateScaffs[counter])
            toWrite.write("\n")

        # increment counter
        counter += 1

# close the file 
toWrite.close()