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
   promoter = len(exclusionWindow)-1 
   enhancer = len(croppedMatrix)-1 
   originalEnhancerSize = len(croppedMatrix) 
   enhancerChr = ''
   enhancerStart = 0
   enhancerEnd = 0
   reachedChromosome = 0
   passedChromosome = 0
   switchDirection = 0
   while enhancer > 0:
       enhancerChr = croppedMatrix[enhancer][0]
       enhancerStart = int(croppedMatrix[enhancer][1])
       enhancerEnd = int(croppedMatrix[enhancer][2])
       promoter = len(exclusionWindow)-1 
       reachedChromosome = 0
       if enhancer < originalEnhancerSize/2:
           switchDirection = 1
       while promoter >= 1 and ~switchDirection:
           if(enhancerChr == exclusionWindow[promoter][0]):
               reachedChromosome = 1
           if(reachedChromosome and enhancerChr != exclusionWindow[promoter][0]):
               break
           if(enhancerChr == exclusionWindow[promoter][0] and enhancerStart>=exclusionWindow[promoter][1] and enhancerEnd<=exclusionWindow[promoter][2]):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop()) 
              break
           promoter-=1
           if(enhancerChr == exclusionWindow[promoter][0] and enhancerStart>=exclusionWindow[promoter][1] and enhancerEnd<=exclusionWindow[promoter][2]):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop()) 
              break
           promoter-=1
       if promoter == 0 and ~switchDirection:
           if(enhancerChr == exclusionWindow[promoter][0] and enhancerStart>=exclusionWindow[promoter][1] and enhancerEnd<=exclusionWindow[promoter][2]):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop()) 
       while promoter < len(exclusionWindow)-1 and switchDirection:
           if(enhancerChr == exclusionWindow[promoter][0]):
               reachedChromosome = 1
           if(reachedChromosome and enhancerChr != exclusionWindow[promoter][0]):
               break
           if(enhancerChr == exclusionWindow[promoter][0] and enhancerStart>=exclusionWindow[promoter][1] and enhancerEnd<=exclusionWindow[promoter][2]):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop()) 
              break
           promoter+=1
           if(enhancerChr == exclusionWindow[promoter][0] and enhancerStart>=exclusionWindow[promoter][1] and enhancerEnd<=exclusionWindow[promoter][2]):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop()) 
              break
           promoter+=1
       if promoter == len(exclusionWindow)-1 and switchDirection:
           if(enhancerChr == exclusionWindow[promoter][0] and enhancerStart>=exclusionWindow[promoter][1] and enhancerEnd<=exclusionWindow[promoter][2]):
              print 'Removing row ', enhancer
              excludedPromoters.append(croppedMatrix.pop()) 
       enhancer -= 1

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
   global outputLocation
   outputLocation = ''
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
      opts, args = getopt.getopt(argv,"",["sourceFile=","excludeTSS=","outputLocation="])
   except getopt.GetoptError:
      print 'Syntax Error. Syntax: \ncallSuperEnhancers.py --sourceFile <csvFile> --excludeTSS <genomeFile> --outputLocation <outputDirectory>)'
      sys.exit(2)
   for opt, arg in opts:
      if opt in ('--sourceFile'):
         inputFile = arg
      elif opt in ('--excludeTSS'):
         genomeFile = arg
      elif opt in ('--outputLocation'):
         outputLocation = arg
   print 'Input file is ', inputFile 
   print 'Genome file is ', genomeFile 
   print 'Output Location is ',outputLocation
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
   stringFragments = []
  #print output as bed  file
   reFormString = '' 
   reFormString = inputFile.replace('.xls','')
   stringFragments = reFormString.rpartition('/') 
   reFormString = stringFragments[2]
   with open(outputLocation + reFormString + '_tagsExtracted.bed', 'a') as outputEnhancers:
       enhancerPrinter = csv.writer(outputEnhancers, delimiter = '\t')
       for i in range(1,len(croppedMatrix)):
            enhancerPrinter.writerow(croppedMatrix[i])

   with open(outputLocation + reFormString + '_removedEnhancers.bed', 'a') as outputEnhancers:
       enhancerPrinter = csv.writer(outputEnhancers, delimiter = '\t')
       for i in range(1,len(excludedPromoters)):
            enhancerPrinter.writerow(excludedPromoters[i])

if __name__ == "__main__":
   main(sys.argv[1:])
