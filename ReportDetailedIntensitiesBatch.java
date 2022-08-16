
import java.net.*;
import java.io.*;
import java.util.*;
import java.text.*;

public class ReportDetailedIntensitiesBatch
{
	public static void main(String[] argv) throws FileNotFoundException, IOException
	{


		if (argv.length!=6)
		{
			System.err.println("Usage: java ReportDetailedIntensitiesBatch <geneName unigene file> <probe intensity file> <gene list> <control samples> <treatment samples> <distal gt proximal ratio, e.g., 0.2>");
			System.exit(1);
		}

		String geneNameFileName = argv[0];
		String intensityFileName = argv[1];
		String geneListFileName = argv[2];
		String controlFileName = argv[3];
		String treatedFileName = argv[4];
		float threshold = (new Float(argv[5])).floatValue();
		

		// create gene name dictionary

		Hashtable unigeneToGenename = new Hashtable(50000);

		BufferedReader dictionaryReader = new BufferedReader(new FileReader(geneNameFileName));

		String dictionaryLine = dictionaryReader.readLine();

		while (dictionaryLine!=null)
		{
			String names[] = split(dictionaryLine,'\t',false);
			if (!unigeneToGenename.containsKey(names[1]))
			{
				unigeneToGenename.put(names[1],names[0]);
			}
			else
			{
				//System.err.println("Duplicate unigene entry: "+names[1]);
			}
			dictionaryLine = dictionaryReader.readLine();
		}

		// read gene list file information and put it in a Hashtable
		BufferedReader geneListReader = new BufferedReader(new FileReader(geneListFileName));

		String geneLine = geneListReader.readLine();

		Hashtable geneList = new Hashtable(500);

		while (geneLine!=null)
		{

			String[] geneInfo = split(geneLine,'\t',false);

			String geneName = geneInfo[0];
			String probeSetName = geneInfo[1];
			String polyASiteName = geneInfo[2];
			geneList.put(probeSetName+"::"+polyASiteName,geneName+"_"+probeSetName+"_"+polyASiteName);
			geneLine = geneListReader.readLine();
		}

		// read control file list and treatment file list
		BufferedReader controlReader = new BufferedReader(new FileReader(controlFileName));

		String controlLine = controlReader.readLine();
		int numControlFile = 0;
		while (controlLine!=null)
		{
			numControlFile++;
			controlLine = controlReader.readLine();
		}

		controlReader.close();

		controlReader = new BufferedReader(new FileReader(controlFileName));

		controlLine = controlReader.readLine();
		String[] controlSamples = new String[numControlFile];
		numControlFile = 0;
		while (controlLine!=null)
		{
			String[] sampleall = split(controlLine,'\\',false);
			String sample = sampleall[sampleall.length-1];
			sampleall = split(sample,'.',false);
			sample = sampleall[0];
			controlSamples[numControlFile++] = sample;
			controlLine = controlReader.readLine();
		}


		// read control file list and treatment file list
		BufferedReader treatedReader = new BufferedReader(new FileReader(treatedFileName));

		String treatedLine = treatedReader.readLine();
		int numTreatedFile = 0;
		while (treatedLine!=null)
		{
			numTreatedFile++;
			treatedLine = treatedReader.readLine();
		}

		treatedReader.close();

		treatedReader = new BufferedReader(new FileReader(treatedFileName));

		treatedLine = treatedReader.readLine();
		String[] treatedSamples = new String[numTreatedFile];
		numTreatedFile = 0;
		while (treatedLine!=null)
		{
			String[] sampleall = split(treatedLine,'\\',false);
			String sample = sampleall[sampleall.length-1];
			sampleall = split(sample,'.',false);
			sample = sampleall[0];
			treatedSamples[numTreatedFile++] = sample;
			treatedLine = treatedReader.readLine();
		}


		// read probe intensity information and process it gene by gene
		BufferedReader intensityReader = new BufferedReader(new FileReader(intensityFileName));

		String intensityLine = intensityReader.readLine();

		while (intensityLine!=null)
		{
			String[] strs = split(intensityLine,'\t',true);

			if (!geneList.containsKey(strs[1]+"::"+strs[2]))
			{
				intensityLine = intensityReader.readLine();
				continue;
			}
			String fileName = (String)(geneList.get(strs[1]+"::"+strs[2]));
			FileOutputStream outFile=new FileOutputStream(fileName+"_detailed.txt");
			DataOutputStream outData=new DataOutputStream(outFile);

			outData.writeBytes("Sample:\t# of Proximal Probes\t# of Distal Probes\tProximal average:\tDistal average:\tProximal/Distal Ratio:\tLOG Proximal/Distal Ratio:\tProximal Median\tProximal Median Dev (avg)\tDistal Median\tDistal Median Dev (avg)\t\r\n"); // 31.1.13

			int numValid = 	(new Integer(strs[5])).intValue();
			int numInvalid = (new Integer(strs[6])).intValue();
			int numControlSamples = 0;
			int numTreatedSamples = 0;

			String stmp = strs[8].substring(1,strs[8].length()-1); // get rid of the opening and closing paranthesis // 29.1.13
		    String[] stmp2 = split(stmp,':',true); // 29.1.13
		    String[] stmp3 = split(stmp2[0],'$',true); // 29.1.13
		    String[] stmp4 = split(stmp2[1],'$',true); // 29.1.13
	     	numControlSamples = stmp3.length; // 29.1.13
	     	numTreatedSamples = stmp4.length; // 29.1.13

		    float[][] controlValidValues = new float[numControlSamples][numValid]; // 29.1.13
		    float[][] controlInvalidValues = new float[numControlSamples][numInvalid]; // 29.1.13
		    float[][] treatedValidValues = new float[numTreatedSamples][numValid]; // 29.1.13
		    float[][] treatedInvalidValues = new float[numTreatedSamples][numInvalid]; // 19.1.13

			boolean[][] controlValidOutlierProbes = new boolean[numControlSamples][numValid];  // 29.1.13
			boolean[][] controlInvalidOutlierProbes = new boolean[numControlSamples][numInvalid];  // 29.1.13
			boolean[][] treatedValidOutlierProbes = new boolean[numTreatedSamples][numValid];  // 29.1.13
			boolean[][] treatedInvalidOutlierProbes = new boolean[numTreatedSamples][numInvalid];  // 29.1.13


		    for (int i = 0; i<numValid; i++)
		    {
		    	String source = strs[8+i*2].substring(1,strs[8+i*2].length()-1); // get rid of the opening and closing paranthesis // 19.1.13
		    	String[] values = split(source,':',true); // 19.1.13
		     	// 19.1.13 String controlValid = values[0].substring(1,values[0].length());
		     	// 19.1.13 String treatedValid = values[1].substring(0,values[1].length()-1);
		     	String[] controlValid = split(values[0],'$',true); // 19.1.13
		     	String[] treatedValid = split(values[1],'$',true); // 19.1.13

		    	// 19.1.13 controlValidValues[i] = (new Float(controlValid)).floatValue();
		    	// 19.1.13 treatedValidValues[i] = (new Float(treatedValid)).floatValue();
		    	for (int j=0;j<controlValid.length;j++) // 19.1.13
		    		controlValidValues[j][i]=(new Float(controlValid[j])).floatValue();  // 29.1.13
		    	for (int j=0;j<treatedValid.length;j++)  // 19.1.13
		    		treatedValidValues[j][i]=(new Float(treatedValid[j])).floatValue();  // 29.1.13
		    }

		    for (int i = 0; i<numInvalid; i++)
		    {
		    	String source = strs[8+2*numValid+i*2].substring(1,strs[8+2*numValid+i*2].length()-1); // get rid of the opening and closing paranthesis // 19.1.13
		    	String[] values = split(source,':',true); // 19.1.13
		     	// 19.1.13 String controlInvalid = values[0].substring(1,values[0].length());
		     	// 19.1.13 String treatedInvalid = values[1].substring(0,values[1].length()-1);
		     	String[] controlInvalid = split(values[0],'$',true); // 19.1.13
		     	String[] treatedInvalid = split(values[1],'$',true); // 19.1.13

		    	// 19.1.13 controlInvalidValues[i] = (new Float(controlValid)).floatValue();
		    	// 19.1.13 treatedInvalidValues[i] = (new Float(treatedValid)).floatValue();
		    	for (int j=0;j<controlInvalid.length;j++) // 19.1.13
		    		controlInvalidValues[j][i]=(new Float(controlInvalid[j])).floatValue();  // 29.1.13
		    	for (int j=0;j<treatedInvalid.length;j++)  // 19.1.13
		    		treatedInvalidValues[j][i]=(new Float(treatedInvalid[j])).floatValue();  // 29.1.13
		    }

			// 31.1.13 store average of median deviations for each sample

			float[] controlProximalMedian = new float[numControlSamples]; // 31.1.13
			float[] controlDistalMedian = new float[numControlSamples]; // 31.1.13
			float[] treatedProximalMedian = new float[numTreatedSamples]; // 31.1.13
			float[] treatedDistalMedian = new float[numTreatedSamples]; // 31.1.13

			float[] controlProximalAverageMedDev = new float[numControlSamples]; // 31.1.13
			float[] controlDistalAverageMedDev = new float[numControlSamples]; // 31.1.13
			float[] treatedProximalAverageMedDev = new float[numTreatedSamples]; // 31.1.13
			float[] treatedDistalAverageMedDev = new float[numTreatedSamples]; // 31.1.13

			// 29.1.13 new filter. find probes that are inconsistent and ignore them. after ignoring check whether we have <2 probes in proximal or distal set. if so, mark that sample as outlier (i.e. invalid)
			// we again use Iglewicz and Hoaglin's method for this purpose // see below.
			for (int j=0;j<numControlSamples;j++) // 29.1.13
			{ // 29.1.13
				// proximal probes first
				// find the median of the probes
				float[] tmp = Arrays.copyOf(controlValidValues[j],controlValidValues[j].length);
				Arrays.sort(tmp);
				float controlValidMedian = tmp[controlValidValues[j].length/2];
				controlProximalMedian[j] = controlValidMedian;

				float[] controlValidDeviations = new float[controlValidValues[j].length];
				for (int i=0;i<numValid;i++)
				{
					controlValidDeviations[i]=(controlValidValues[j][i]>controlValidMedian)?(controlValidValues[j][i]-controlValidMedian):(controlValidMedian-controlValidValues[j][i]);
				}
				tmp = Arrays.copyOf(controlValidDeviations,controlValidDeviations.length);
				Arrays.sort(tmp);
				float controlValidMAD = tmp[controlValidDeviations.length/2];
				int cnt = 0;
				controlProximalAverageMedDev[j]=0;
				for (int i=0;i<numValid;i++)
				{
					float mi = 0.6745f*(controlValidValues[j][i]-controlValidMedian)/controlValidMAD;
					if (mi>3.5 || mi<-3.5) controlValidOutlierProbes[j][i] = true;
				 	else {
				 		controlValidOutlierProbes[j][i] = false;
				 		controlProximalAverageMedDev[j]+=(controlValidValues[j][i]>controlValidMedian)?(controlValidValues[j][i]-controlValidMedian):(controlValidMedian-controlValidValues[j][i]);
				 		cnt++;
				 	}
				}
				controlProximalAverageMedDev[j]/=cnt;

				// distal probes next
				// find the median of the probes
				tmp = Arrays.copyOf(controlInvalidValues[j],controlInvalidValues[j].length);
				Arrays.sort(tmp);
				float controlInvalidMedian = tmp[controlInvalidValues[j].length/2];
				controlDistalMedian[j] = controlInvalidMedian;

				float[] controlInvalidDeviations = new float[controlInvalidValues[j].length];
				for (int i=0;i<numInvalid;i++)
				{
					controlInvalidDeviations[i]=(controlInvalidValues[j][i]>controlInvalidMedian)?(controlInvalidValues[j][i]-controlInvalidMedian):(controlInvalidMedian-controlInvalidValues[j][i]);
				}
				tmp = Arrays.copyOf(controlInvalidDeviations,controlInvalidDeviations.length);
				Arrays.sort(tmp);
				float controlInvalidMAD = tmp[controlInvalidDeviations.length/2];
				cnt = 0;
				controlDistalAverageMedDev[j]=0;
				for (int i=0;i<numInvalid;i++)
				{
					float mi = 0.6745f*(controlInvalidValues[j][i]-controlInvalidMedian)/controlInvalidMAD;
					if (mi>3.5 || mi<-3.5) controlInvalidOutlierProbes[j][i] = true;
				 	else
				 	{
				 		controlInvalidOutlierProbes[j][i] = false;
				 		controlDistalAverageMedDev[j]+=(controlInvalidValues[j][i]>controlInvalidMedian)?(controlInvalidValues[j][i]-controlInvalidMedian):(controlInvalidMedian-controlInvalidValues[j][i]);
						cnt++;
					}
				}
				controlDistalAverageMedDev[j]/=cnt;
			} // 29.1.13 all of the for loop is added on 29.1.13



			for (int j=0;j<numTreatedSamples;j++) // 29.1.13
			{ // 29.1.13
				// proximal probes first
				// find the median of the probes
				float[] tmp = Arrays.copyOf(treatedValidValues[j],treatedValidValues[j].length);
				Arrays.sort(tmp);
				float treatedValidMedian = tmp[treatedValidValues[j].length/2];
				treatedProximalMedian[j] = treatedValidMedian;

				float[] treatedValidDeviations = new float[treatedValidValues[j].length];
				for (int i=0;i<numValid;i++)
				{
					treatedValidDeviations[i]=(treatedValidValues[j][i]>treatedValidMedian)?(treatedValidValues[j][i]-treatedValidMedian):(treatedValidMedian-treatedValidValues[j][i]);
				}
				tmp = Arrays.copyOf(treatedValidDeviations,treatedValidDeviations.length);
				Arrays.sort(tmp);
				float treatedValidMAD = tmp[treatedValidDeviations.length/2];
				int cnt = 0;
				treatedProximalAverageMedDev[j]=0;
				for (int i=0;i<numValid;i++)
				{
					float mi = 0.6745f*(treatedValidValues[j][i]-treatedValidMedian)/treatedValidMAD;
					if (mi>3.5 || mi<-3.5) treatedValidOutlierProbes[j][i] = true;
				 	else {
				 		treatedValidOutlierProbes[j][i] = false;
				 		treatedProximalAverageMedDev[j]+=(treatedValidValues[j][i]>treatedValidMedian)?(treatedValidValues[j][i]-treatedValidMedian):(treatedValidMedian-treatedValidValues[j][i]);
				 		cnt++;
				 	}
				}
				treatedProximalAverageMedDev[j]/=cnt;

				// distal probes next
				// find the median of the probes
				tmp = Arrays.copyOf(treatedInvalidValues[j],treatedInvalidValues[j].length);
				Arrays.sort(tmp);
				float treatedInvalidMedian = tmp[treatedInvalidValues[j].length/2];
				treatedDistalMedian[j] = treatedInvalidMedian;

				float[] treatedInvalidDeviations = new float[treatedInvalidValues[j].length];
				for (int i=0;i<numInvalid;i++)
				{
					treatedInvalidDeviations[i]=(treatedInvalidValues[j][i]>treatedInvalidMedian)?(treatedInvalidValues[j][i]-treatedInvalidMedian):(treatedInvalidMedian-treatedInvalidValues[j][i]);
				}
				tmp = Arrays.copyOf(treatedInvalidDeviations,treatedInvalidDeviations.length);
				Arrays.sort(tmp);
				float treatedInvalidMAD = tmp[treatedInvalidDeviations.length/2];
				cnt = 0;
				treatedDistalAverageMedDev[j] = 0;
				for (int i=0;i<numInvalid;i++)
				{
					float mi = 0.6745f*(treatedInvalidValues[j][i]-treatedInvalidMedian)/treatedInvalidMAD;
					if (mi>3.5 || mi<-3.5) treatedInvalidOutlierProbes[j][i] = true;
				 	else
				 	{
				 		treatedInvalidOutlierProbes[j][i] = false;
				 		treatedDistalAverageMedDev[j]+=(treatedInvalidValues[j][i]>treatedInvalidMedian)?(treatedInvalidValues[j][i]-treatedInvalidMedian):(treatedInvalidMedian-treatedInvalidValues[j][i]);
				 		cnt++;
				 	}
				}
				treatedDistalAverageMedDev[j]/=cnt;
			} // 29.1.13 all of the for loop is added on 29.1.13


			// 29.1.13 now find samples with <2 proximal/distal probes and mark them as ignore

			boolean controlOutliers[] = new boolean[numControlSamples]; // 29.1.13
			boolean treatedOutliers[] = new boolean[numTreatedSamples]; // 29.1.13


			for (int j=0;j<numControlSamples;j++) // 29.1.13
			{ // 29.1.13
				int cnt = 0;
				for (int i=0;i<numValid;i++)
					if (controlValidOutlierProbes[j][i]==false) cnt++;
				if (cnt<2) {
					controlOutliers[j] = true;
					continue;
				} else controlOutliers[j] = false;

				cnt = 0;
				for (int i=0;i<numInvalid;i++)
					if (controlInvalidOutlierProbes[j][i]==false) cnt++;
				if (cnt<2)
					controlOutliers[j] = true;
			} // 29.1.13 all of the for loop is added on 29.1.13


			for (int j=0;j<numTreatedSamples;j++) // 29.1.13
			{ // 29.1.13
				int cnt = 0;
				for (int i=0;i<numValid;i++)
					if (treatedValidOutlierProbes[j][i]==false) cnt++;
				if (cnt<2) {
					treatedOutliers[j] = true;
					continue;
				} else treatedOutliers[j] = false;

				cnt = 0;
				for (int i=0;i<numInvalid;i++)
					if (treatedInvalidOutlierProbes[j][i]==false) cnt++;
				if (cnt<2)
					treatedOutliers[j] = true;
			} // 29.1.13 all of the for loop is added on 29.1.13


			// 29.1.13 count number of usable control and treated samples
			int numUsableControlSamples = 0;
			int numUsableTreatedSamples = 0;

			for (int j=0;j<numControlSamples;j++)
			{
				if (controlOutliers[j]==false) numUsableControlSamples++;
			}
			for (int j=0;j<numTreatedSamples;j++)
			{
				if (treatedOutliers[j]==false) numUsableTreatedSamples++;
			}
			if (numUsableControlSamples<1 || numUsableTreatedSamples<1)
			{
				intensityLine = intensityReader.readLine();
				continue;
			}
			// 29.1.13 all of the above code until the comment is added on 29.1.13


			// 19.1.13 Now find the meanControlValid and meanControlInvalid for each sample and remove samples that are outliers whem computing the overall averages
			float meanControlValidSample[] = new float[numControlSamples]; // 19.1.13
			float meanTreatedValidSample[] = new float[numTreatedSamples]; // 19.1.13
			float meanControlInvalidSample[] = new float[numControlSamples]; // 19.1.13
			float meanTreatedInvalidSample[] = new float[numTreatedSamples]; // 19.1.13

			for (int j=0;j<numControlSamples;j++) // 19.1.13
			{ // 19.1.13
				meanControlValidSample[j] = 0; // 19.1.13
				meanControlInvalidSample[j] = 0; // 19.1.13
				int cnt = 0; // 29.1.13
				for (int i=0;i<numValid;i++) // 19.1.13
				{ // 19.1.13
					if (controlValidOutlierProbes[j][i]==false) { meanControlValidSample[j]+=controlValidValues[j][i]; cnt++; } // 29.1.13
				} // 19.1.13
				meanControlValidSample[j]/=cnt; // 29.1.13
				cnt = 0; // 29.1.13
				for (int i=0;i<numInvalid;i++) // 19.1.13
				{ // 19.1.13
					if (controlInvalidOutlierProbes[j][i]==false) { meanControlInvalidSample[j]+=controlInvalidValues[j][i]; cnt++; } // 29.1.13
				} // 19.1.13
				meanControlInvalidSample[j]/=cnt; // 29.1.13
			} // 19.1.13

			for (int j=0;j<numTreatedSamples;j++) // 19.1.13
			{ // 19.1.13
				meanTreatedValidSample[j] = 0; // 19.1.13
				meanTreatedInvalidSample[j] = 0; // 19.1.13
				int cnt = 0; // 29.1.13
				for (int i=0;i<numValid;i++) // 19.1.13
				{ // 19.1.13
					if (treatedValidOutlierProbes[j][i]==false) { meanTreatedValidSample[j]+=treatedValidValues[j][i]; cnt++; } // 29.1.13
				} // 19.1.13
				meanTreatedValidSample[j]/=cnt; // 29.1.13
				cnt = 0; // 29.1.13
				for (int i=0;i<numInvalid;i++) // 19.1.13
				{ // 19.1.13
					if (treatedInvalidOutlierProbes[j][i]==false) { meanTreatedInvalidSample[j]+=treatedInvalidValues[j][i]; cnt++; } // 29.1.13
				} // 19.1.13
				meanTreatedInvalidSample[j]/=cnt; // 29.1.13
			} // 19.1.13


			// 19.1.13 Now find samples (in the control group and treated group) that are outliers. For this we use Iglewicz and Hoaglin's method which recommends using the modified Z-score
            // M(i) = 0.6745*(x(i) - xtilde)/MAD with MAD denoting the median absolute deviation and xtilde denoting the median.
            // These authors recommend that modified Z-scores with an absolute value of greater than 3.5 be labeled as potential outliers.

			float  controlRatioMedian = 0; // 19.1.13
			float[]  controlRatios = new float[numUsableControlSamples]; // 29.1.13
			int usableCount = 0;
			for (int j=0;j<numControlSamples;j++) // 19.1.13
			{ // 19.1.13
				if (controlOutliers[j] == false) // 29.1.13
					controlRatios[usableCount++]=(meanControlValidSample[j]/meanControlInvalidSample[j]); // 29.1.13
			} // 19.1.13
			float[] tmp = Arrays.copyOf(controlRatios,controlRatios.length);
			Arrays.sort(tmp); // 19.1.13
			controlRatioMedian = tmp[controlRatios.length/2]; // 29.1.13

			float controlMAD = 0; // 19.1.13
			float[] controlRatioDeviations = new float[numUsableControlSamples]; // 29.1.13
			usableCount = 0; // 29.1.13
			for (int j=0;j<numControlSamples;j++) // 19.1.13
			{ // 19.1.13
				if (controlOutliers[j] == false) // 29.1.13
					controlRatioDeviations[usableCount]=(controlRatios[usableCount]>controlRatioMedian)?(controlRatios[usableCount]-controlRatioMedian):(controlRatioMedian-controlRatios[usableCount++]); // 29.1.13
			} // 19.1.13
			tmp = Arrays.copyOf(controlRatioDeviations,controlRatioDeviations.length);
			Arrays.sort(tmp); // 19.1.13
			controlMAD = tmp[controlRatioDeviations.length/2]; // 29.1.13

			float  treatedRatioMedian = 0; // 19.1.13
			float[]  treatedRatios = new float[numUsableTreatedSamples]; // 29.1.13
			usableCount = 0;
			for (int j=0;j<numTreatedSamples;j++) // 19.1.13
			{ // 19.1.13
				if (treatedOutliers[j] == false) // 29.1.13
					treatedRatios[usableCount++]=(meanTreatedValidSample[j]/meanTreatedInvalidSample[j]); // 29.1.13
			} // 19.1.13
			tmp = Arrays.copyOf(treatedRatios,treatedRatios.length);
			Arrays.sort(tmp); // 19.1.13
			treatedRatioMedian = tmp[treatedRatios.length/2]; // 29.1.13

			float treatedMAD = 0; // 19.1.13
			float[] treatedRatioDeviations = new float[numUsableTreatedSamples]; // 29.1.13
			usableCount = 0; // 29.1.13
			for (int j=0;j<numTreatedSamples;j++) // 19.1.13
			{ // 19.1.13
				if (treatedOutliers[j] == false) // 29.1.13
				treatedRatioDeviations[usableCount]=(treatedRatios[usableCount]>treatedRatioMedian)?(treatedRatios[usableCount]-treatedRatioMedian):(treatedRatioMedian-treatedRatios[usableCount++]); // 29.1.13
			} // 19.1.13
			tmp = Arrays.copyOf(treatedRatioDeviations,treatedRatioDeviations.length);
			Arrays.sort(tmp); // 19.1.13
			treatedMAD = tmp[treatedRatioDeviations.length/2]; // 19.1.13

			usableCount = 0; // 29.1.13
			for (int j=0;j<numControlSamples;j++) // 19.1.13
			{ // 19.1.13
				if (controlOutliers[j] == false) // 29.1.13
				{
					float mi = 0.6745f*(controlRatios[usableCount++]-controlRatioMedian)/controlMAD; // 29.1.13
					// if (mi>3.5 || mi<-3.5) controlOutliers[j] = true; // 29.1.13 // removed this filter 31.1.13
				}
			} // 19.1.13
			usableCount = 0; // 29.1.13
			for (int j=0;j<numTreatedSamples;j++) // 19.1.13
			{ // 19.1.13
				if (treatedOutliers[j] == false) // 29.1.13
				{
					float mi = 0.6745f*(treatedRatios[usableCount++]-treatedRatioMedian)/treatedMAD; // 29.1.13
					// if (mi>3.5 || mi<-3.5) treatedOutliers[j] = true; // 29.1.13 // removed this filter 31.1.13
				}
			} // 19.1.13



			for (int j=0;j<numControlSamples;j++)
			{
				if (controlOutliers[j]==true) continue;
			    float stdevControlValid = 0.0f;
			    float stdevControlInvalid = 0.0f;
			    int cnt = 0; // 29.1.13
			    for (int i=0;i<controlValidValues[j].length;i++)
			    {
			    	if (controlValidOutlierProbes[j][i]==false) // 29.1.13
			    	{ // 29.1.13
				    	stdevControlValid+=(controlValidValues[j][i]-meanControlValidSample[j])*(controlValidValues[j][i]-meanControlValidSample[j]);
			    	    cnt++; // 29.1.13
			    	} // 29.1.13
			    }
			    stdevControlValid/=cnt; // 29.1.13
		    	if (stdevControlValid!=0) stdevControlValid=(float)Math.sqrt(stdevControlValid);
		    	stdevControlValid/=meanControlValidSample[j]; // normalized stdev

				cnt = 0; // 29.1.13
			    for (int i=0;i<controlInvalidValues[j].length;i++)
			    {
			    	if (controlInvalidOutlierProbes[j][i]==false) // 29.1.13
			    	{ // 29.1.13
				    	stdevControlInvalid+=(controlInvalidValues[j][i]-meanControlInvalidSample[j])*(controlInvalidValues[j][i]-meanControlInvalidSample[j]);
				    	cnt++;
			    	} // 29.1.13
			    }
			    stdevControlInvalid/=cnt; // 29.1.13
		    	if (stdevControlInvalid!=0) stdevControlInvalid=(float)Math.sqrt(stdevControlInvalid);
		    	stdevControlInvalid/=meanControlInvalidSample[j]; // normalized stdev
		        if ( // removed this filter 01.02.13 stdevControlValid < 0.50 && stdevControlInvalid < 0.50 &&
		        	(meanControlInvalidSample[j]-meanControlValidSample[j]) < (threshold*meanControlValidSample[j])) // changed threshold to 20% 01.02.13 // check if the intensities are consistent within samples and also look for irregular levels such as long > short
						controlOutliers[j] = false;
				else
						controlOutliers[j] = true; // removed this filter 31.1.13 // added this filter again on 01.02.13
			}

			for (int j=0;j<numTreatedSamples;j++)
			{
				if (treatedOutliers[j]==true) continue;
			    float stdevTreatedValid = 0.0f;
			    float stdevTreatedInvalid = 0.0f;
			    int cnt = 0; // 29.1.13
			    for (int i=0;i<treatedValidValues[j].length;i++)
			    {
			    	if (treatedValidOutlierProbes[j][i]==false) // 29.1.13
			    	{
				    	stdevTreatedValid+=(treatedValidValues[j][i]-meanTreatedValidSample[j])*(treatedValidValues[j][i]-meanTreatedValidSample[j]);
				    	cnt++; // 29.1.13
			    	}
			    }
			    stdevTreatedValid/=cnt;
		    	if (stdevTreatedValid!=0) stdevTreatedValid=(float)Math.sqrt(stdevTreatedValid);
		    	stdevTreatedValid/=meanTreatedValidSample[j]; // normalized stdev

				cnt = 0; // 29.1.13
			    for (int i=0;i<treatedInvalidValues[j].length;i++)
			    {
			    	if (treatedInvalidOutlierProbes[j][i]==false) // 29.1.13
			    	{
				    	stdevTreatedInvalid+=(treatedInvalidValues[j][i]-meanTreatedInvalidSample[j])*(treatedInvalidValues[j][i]-meanTreatedInvalidSample[j]);
				    	cnt++; // 29.1.13
				    }
			    }
			    stdevTreatedInvalid/=cnt;
		    	if (stdevTreatedInvalid!=0) stdevTreatedInvalid=(float)Math.sqrt(stdevTreatedInvalid);
		    	stdevTreatedInvalid/=meanTreatedInvalidSample[j]; // normalized stdev
		        if (// removed this filter 01.02.13 stdevTreatedValid < 0.50 && stdevTreatedInvalid < 0.50 &&
		        	(meanTreatedInvalidSample[j]-meanTreatedValidSample[j]) < (threshold*meanTreatedValidSample[j])) // changed threshold to 20% 01.02.13 // check if the intensities are consistent within samples and also look for irregular levels such as long > short
						treatedOutliers[j] = false;
				else
						treatedOutliers[j] = true; // removed this filter 31.1.13 // added this filter again on 01.02.13
			}


			float controlratioAvg = 0;
			float treatedratioAvg = 0;
			int controlCnt = 0;
			int treatedCnt = 0;

			DecimalFormat df = new DecimalFormat("#0.00");

			outData.writeBytes("----Control Samples----\r\n");
			for (int j=0;j<numControlSamples;j++)
			{
				if (controlOutliers[j]==false)
				{

					int pcnt = 0;
					for (int i=0;i<numValid;i++)
						if (controlValidOutlierProbes[j][i]==false) pcnt++;

					int dcnt = 0;
					for (int i=0;i<numInvalid;i++)
						if (controlInvalidOutlierProbes[j][i]==false) dcnt++;

					outData.writeBytes(controlSamples[j]+"\t"+pcnt+"\t"+dcnt+"\t"+df.format(meanControlValidSample[j])+"\t"+df.format(meanControlInvalidSample[j])+"\t"+df.format(controlRatios[j])+"\t"+df.format(Math.log(controlRatios[j]))+"\t"+df.format(controlProximalMedian[j])+"\t"+df.format(controlProximalAverageMedDev[j])+"\t"+df.format(controlDistalMedian[j])+"\t"+df.format(controlDistalAverageMedDev[j])+"\r\n");
					controlratioAvg += meanControlValidSample[j]/meanControlInvalidSample[j]; // 29.1.13
					controlCnt++;
				}
			}
			if (controlCnt==0)
			{
				System.out.println("All control samples are identified as outliers for "+strs[1]+"::"+strs[2]+"!");
				System.exit(1);
			}
			controlratioAvg/=controlCnt;

			outData.writeBytes("Average:\t \t \t"+df.format(controlratioAvg)+"\r\n");
			outData.writeBytes("----Treated Samples----\r\n");
			for (int j=0;j<numTreatedSamples;j++)
			{
				if (treatedOutliers[j]==false)
				{
					int pcnt = 0;
					for (int i=0;i<numValid;i++)
						if (treatedValidOutlierProbes[j][i]==false) pcnt++;

					int dcnt = 0;
					for (int i=0;i<numInvalid;i++)
						if (treatedInvalidOutlierProbes[j][i]==false) dcnt++;

					outData.writeBytes(treatedSamples[j]+"\t"+pcnt+"\t"+dcnt+"\t"+df.format(meanTreatedValidSample[j])+"\t"+df.format(meanTreatedInvalidSample[j])+"\t"+df.format(treatedRatios[j])+"\t"+df.format(Math.log(treatedRatios[j]))+"\t"+df.format(treatedProximalMedian[j])+"\t"+df.format(treatedProximalAverageMedDev[j])+"\t"+df.format(treatedDistalMedian[j])+"\t"+df.format(treatedDistalAverageMedDev[j])+"\r\n");
					treatedratioAvg += meanTreatedValidSample[j]/meanTreatedInvalidSample[j]; // 29.1.13
					treatedCnt++;
				}
			}
			if (treatedCnt==0)
			{
				System.out.println("All treatment samples are identified as outliers for "+strs[1]+"::"+strs[2]+"!");
				System.exit(1);
			}
			treatedratioAvg/=treatedCnt;

			outData.writeBytes("Average:\t \t \t"+df.format(treatedratioAvg)+"\r\n");

			double ratio = Math.log((double)treatedratioAvg/controlratioAvg);

			outData.writeBytes("LOG Treated/Control Ratio: "+df.format(ratio)+"\r\n");

	    	outData.close();
			intensityLine = intensityReader.readLine();

		}

		intensityReader.close();

    }

   	public static String[] split(String str, char delim, boolean zeroSizeOK)
    {
                // begin split
                Vector strsVec = new Vector(0,1);
                String tmp = str;
                while (tmp.indexOf(delim)!=-1)
                {
                		if (zeroSizeOK || tmp.substring(0,tmp.indexOf(delim)).length()>0)
                        	strsVec.addElement(new String(tmp.substring(0,tmp.indexOf(delim))));
                        tmp = tmp.substring(tmp.indexOf(delim)+1,tmp.length());
                }
                if (zeroSizeOK || tmp.length()>0)
                	strsVec.addElement(new String(tmp));
                String[] strs = new String[strsVec.capacity()];
                for (int s = 0; s < strsVec.capacity(); s++)
                        strs[s] = (String)strsVec.elementAt(s);
                // end of split
                return strs;
    }

   	public static String[] split(String str, String delim, boolean zeroSizeOK)
    {
                // begin split
                Vector strsVec = new Vector(0,1);
                String tmp = str;
                while (tmp.indexOf(delim)!=-1)
                {
                		if (zeroSizeOK || tmp.substring(0,tmp.indexOf(delim)).length()>0)
                        	strsVec.addElement(new String(tmp.substring(0,tmp.indexOf(delim))));
                        tmp = tmp.substring(tmp.indexOf(delim)+delim.length(),tmp.length());
                }
                if (zeroSizeOK || tmp.length()>0)
	                strsVec.addElement(new String(tmp));
                String[] strs = new String[strsVec.capacity()];
                for (int s = 0; s < strsVec.capacity(); s++)
                        strs[s] = (String)strsVec.elementAt(s);
                // end of split
                return strs;
    }


}