#!/usr/bin/python
#Written by Kamin Kahrizi, CHORI
import sys, getopt, csv, itertools

def mergeEnhancers():
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
   with open(inputFile, 'rb') as enhancerPeaksAndTags:
        enhancerReader = csv.reader(enhancerPeaksAndTags, dialect='excel')
        counter = 0
        for row in enhancerPeaksAndTags.readlines():
            array.append(row)
            if row == 'chr\tstart\tend\tlength\tsummit\ttags\t-10*log10(pvalue)\tfold_enrichment\tFDR(%)\n':
                firstDataLine = counter
            counter += 1
        lastDataLine = len(array)-1

        print array[firstDataLine] 
        print array[lastDataLine]
        
        for row in range(firstDataLine,lastDataLine+1):
            matrix.append(array[row].split('\t'))
            croppedMatrix.append([matrix[row-firstDataLine][chrCol],matrix[row-firstDataLine][startCol],matrix[row-firstDataLine][endCol],matrix[row-firstDataLine][tagCol]])

        print "Length of cropped matrix: ", len(croppedMatrix)
   
   print croppedMatrix[0]
   print croppedMatrix[1]
   print croppedMatrix[2]
   print croppedMatrix[len(croppedMatrix)-1]
   inputParse = inputFile.split('.')
   print 'Length of input parse ', len(inputParse)
   reFormString = '' 
   index = 0
   while index < (len(inputParse)-1):
       if(inputParse[index] == ''):
           reFormString = reFormString + '.'
       else:
           reFormString = reFormString + inputParse[index]
       index += 1

   with open(reFormString + '_tagsExtracted.csv', 'a') as outputEnhancers:
       enhancerPrinter = csv.writer(outputEnhancers, delimiter = '\t')
       for i in range(1,len(croppedMatrix)):
            enhancerPrinter.writerow(croppedMatrix[i])

if __name__ == "__main__":
   main(sys.argv[1:])
