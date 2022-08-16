# APADetect
Detection of Alternative Polyadenylation events by probe level analysis of Affymetrix Microarrays

APADetect Command Line Version Running Instructions

1. Requirements

In order to compile the java source files, Java Development Kit (JDK) is needed.

	- Java Platform (JDK) can be downloaded from:

	  http://www.oracle.com/technetwork/java/javase/downloads/index.html
	  Select JDK to download and install it in your PC.

	- During installation, note the installation directory. By default it is installed in:
          
          C:\Program Files\Java\jdk1.<version>

	        the "bin" directory in this directory contains the "java" runtime and "javac" compiler.
          You will run these programs from the command prompt for compiling and running APADetect.
          At the command prompt, you may need to type the whole path name such as:
          "C:\Program Files\Java\jdk1.7.0_45\bin\javac" if "C:\Program Files\Java\jdk1.7.0_45\bin" is not
          included in your PATH environment variable. Run command prompt from Accessories or by
          executing "cmd.exe" if typing "javac" only results in a response such as "Command not found"
	        it means that the "bin" directory is not included in your PATH variable. Type the whole
          path name as above or adding the "C:\Program Files\Java\jdk1.7.0_45\bin" directory at the end of the
	        PATH variable by selecting Computer->Properties->Advanced system settings->Advanced Tab-> Environment Variables->Edt (PATH).

2. Compilation

- Unzip the contents of the APADetect zip bundle. On Windows machines, open the command prompt by typing "cmd.exe" at the start menu.

- Go to the APADetect directory using the "cd" command at the command prompt.

- type "javac GetProbeIntensities.java"

- type "javac ProcessProbeIntensities.java"

- type "javac ReportDetailedIntensitiesBatch.java"

3. Running APADetect

To run APADetect, you need to sets of input files:
	1. The list of control samples (as a list of GSM ids line by line in a separate text file) and
     the list of treated samples (as a list of GSM ids line by line in a separate text file).

	2. The raw CEL files of the listed samples in ASCII format. To find out whether the CEL files you have are in ASCII format,
     try to open the CEL file using WordPad (or another text edfitor). If you see strange unreadable characters,
	   it means that the CEL file is NOT in ASCII format. Convert it to ASCII using the apt-cel-convert.exe,
	   provided in the zip bundle. Refer to apa-cel-convert help file for usage. Sample usage:
	   apt-cel-convert -f text -o text-cels GSM213332.CEL. Copy the converted CEL files in text-cels directory
	   to the same directory as the compiled java class files. If the CEL files for some samples are missing,
     the first step described below tries to download them from NCBI GEO ftp site automatically. But the
	   downloaded file is a gz (zip) file. You need to unzip the downloaded files and re-run the first step.

1. The first step is to run the GetProbeIntensties program. Run it by:
   java GetProbeIntensities splitProbeSetsU133A.txt splitProbeSetsU133Plus2.txt <control file> <treated file> <output>

   This program processes CEL files and produces a single intensity file as output.

   If your dataset contains CEL files from a single platform (U133A or U133Plus2.0 only) then use the same platfrom name as the first
   two argumenst to the program.
   Use different paramter names as in the example above only if your dataset contains samples from both platforms.

   To see the usage of the program run it without providing any parameters, like:
   java GetProbeIntensities

2. The second step is to process the intensity output file produced by the first step using the ProcessProbeIntensities program
   Run it by:
   java ProcessProbeIntensities <geneName unigene file> <probe intensity file, e.g., intensities.txt> <control file> <treated file> <output file, e.g., results.txt> <distal gt proximal ratio, e.g., 0.2>

   To see the usage of the program run it without providing any parameters, like:
   java ProcessProbeIntensities

   This program processes the intenstity output from the first step. Unigene->Gene name conversion file is provided in the zip bundle.
   The last paramter is a quality threshold which filters out genes with unexpected intensities when the distal probes have higher
   intensities that the proximal probes. If this parameter is t, then the gene is removed from the result if its distal probes have 
   t*proximal average intensity greater average intensity than the proximal average intensity. Typical value of t is 0.2.

   This program outputs 3 files. The output file specified in the parameters is a table with proximal/distal intensity ratios for genes.
   A second file contains log(proximal/distal ratio) values for every gene for every sample. That matrix can be further analyzed for statistically
   significantly different genes between two groups (i.e., control versus treated) by using program such as TIGR MEV.
   The third output file is a list of gene names, probeset IDs, and polyA site IDs which can be used by the third program as detailed next.

   The first and second outputs can be opened in Microsoft Excel directly since they are tab separated files.    

3. The third step is optional. You may run ReportDetailedIntensitiesBatch to get detailed sample level proximal/distal intensity statistics for each gene
   in a given gene list. If the gene list output by the second step is used, detailed output files are produced for all the genes in the study (about 2000    files may be generated
   by this program). You may execute this program for only a selected list of genes of interest.
   Run this program by:

   java ReportDetailedIntensitiesBatch <geneName unigene file> <probe intensity file> <gene list> <control samples> <treatment samples> <distal gt proximal ratio, e.g., 0.2>
