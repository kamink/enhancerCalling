#!/usr/bin/python
#Written by Kamin Kahrizi, CHORI
import sys, getopt, csv, itertools

#removes enhancers wholly contained within a window of size windowSize centered at the tss of each gene
def excludePromoters(windowSize):
   global chrCol 
   global startCol
   global endCol 
   global tagCol
   global inputFile
   global outputFile
   global croppedMatrix
   global genomeMatrix
   global exclusionWindow 
   global excludedPromoters
   strandCol = 5
   for row in range(0,len(genomeMatrix)):
       if genomeMatrix[row][strandCol] == '+\n': 
          startSite = genomeMatrix[row][startCol]
       elif genomeMatrix[row][strandCol] == '-\n':
          startSite = genomeMatrix[row][endCol]
       exclusionWindow.append([genomeMatrix[row][chrCol],int(startSite)-windowSize/2,int(startSite)+windowSize/2])

   print 'Cropped Matrix before removing promoters: ', '\n'.join([ str(myElement) for myElement in croppedMatrix[0:50]])
#   print exclusionWindow[0]
#   print exclusionWindow[1]
#   print exclusionWindow[2]
#   print exclusionWindow[len(exclusionWindow)-1]
   promoter = 0 
   enhancer = 1 
   promoterChr = ''
   promoterStart = 0
   promoterEnd = 0
   while promoter < len(exclusionWindow):
       promoterChr = exclusionWindow[promoter][0]
       promoterStart = exclusionWindow[promoter][1]
       promoterEnd = exclusionWindow[promoter][2]
       enhancer = 0 
       while enhancer < len(croppedMatrix)-1:
           if(promoterChr == croppedMatrix[enhancer][0] and int(croppedMatrix[enhancer][1])>=promoterStart and int(croppedMatrix[enhancer][2])<=promoterEnd):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop(enhancer)) 
           enhancer+=1
           if(promoterChr == croppedMatrix[enhancer][0] and int(croppedMatrix[enhancer][1])>=promoterStart and int(croppedMatrix[enhancer][2])<=promoterEnd):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop(enhancer)) 
           enhancer+=1
       if enhancer < len(croppedMatrix):
           if(promoterChr == croppedMatrix[enhancer][0] and int(croppedMatrix[enhancer][1])>=promoterStart and int(croppedMatrix[enhancer][2])<=promoterEnd):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop(enhancer)) 
       promoter+=1
#merges enhancers within mergeDistance
def mergeEnhancers(mergeDistance):
   global chrCol 
   global startCol
   global endCol 
   global tagCol
   global inputFile
   global outputFile
   global croppedMatrix

   counter = 1
   distance = 0 
   thisChr, thisStartCol, thisEndCol, thisTagCol = '', '', '', ''
   nextChr, nextStartCol, nextEndCol, nextTagCol = '', '', '', ''

#   while counter < len(croppedMatrix)-1
#             

def main(argv):
   global firstDataLine 
   global lastDataLine 
   global firstDataCol 
   firstDatacol = 1
   global lastDataCol
   lastDataCol = 9
   global chrCol 
   chrCol = 0
   global startCol
   startCol = 1
   global endCol 
   endCol = 2  
   global tagCol
   tagCol = 5
   global inputFile
   inputFile = ''
   global genomeFile
   genomeFile = ''
   global array
   array = []
   global matrix
   matrix = []
   global croppedMatrix
   croppedMatrix = []  
   global genomeMatrix 
   genomeMatrix = []
   global exclusionWindow
   exclusionWindow = []
   global excludedPromoters 
   excludedPromoters = []
   try:
      opts, args = getopt.getopt(argv,"",["sourceFile=","excludeTSS="])
   except getopt.GetoptError:
      print 'Syntax Error. Syntax: \ncallSuperEnhancers.py --sourceFile <csvFile> --excludeTSS <genomeFile>)'
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('--sourceFile'):
         inputFile = arg
      elif opt in ('--excludeTSS'):
         genomeFile = arg
   print 'Input file is ', inputFile 
   print 'Genome file is ', genomeFile 

  #open source file and store as array in croppedMatrix (only chromosomal data)
   with open(inputFile, 'rb') as enhancerPeaksAndTags:
        counter = 0
        for row in enhancerPeaksAndTags.readlines():
            array.append(row)
            if row == 'chr\tstart\tend\tlength\tsummit\ttags\t-10*log10(pvalue)\tfold_enrichment\tFDR(%)\n':
                firstDataLine = counter
            counter += 1
        lastDataLine = len(array)-1

        for row in range(firstDataLine,lastDataLine+1):
            matrix.append(array[row].split('\t'))
            croppedMatrix.append([matrix[row-firstDataLine][chrCol],matrix[row-firstDataLine][startCol],matrix[row-firstDataLine][endCol],matrix[row-firstDataLine][tagCol]])

  #open genome file and store as array in genomeMatrix
   array = []
   matrix = []
   with open(genomeFile, 'rb') as genomeFile:
        for row in genomeFile.readlines():
            array.append(row)

        for row in range(0,len(array)):
            genomeMatrix.append(array[row].split('\t'))

   print genomeMatrix[0]
   print genomeMatrix[1]
   print genomeMatrix[2]
   print genomeMatrix[len(genomeMatrix)-1]
 
   excludePromoters(2000)

  #print output as bed  file
   reFormString = '' 
   index = 0
   reFormString = inputFile.replace('.xls','')

   with open(reFormString + '_tagsExtracted.bed', 'a') as outputEnhancers:
       enhancerPrinter = csv.writer(outputEnhancers, delimiter = '\t')
       for i in range(1,len(croppedMatrix)):
            enhancerPrinter.writerow(croppedMatrix[i])

if __name__ == "__main__":
   main(sys.argv[1:])
