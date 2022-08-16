
import java.net.*;
import java.io.*;
import java.util.*;
import java.text.*;

public class GetProbeIntensities
{
	public static void main(String[] argv) throws FileNotFoundException, IOException
	{


		if (argv.length!=5)
		{
			System.err.println("Usage: java GetProbeIntensities <split probe list output file for U133A> <split probe list output file for U133 Plus 2.0> <control file list> <treated file list> <output file>");
			System.exit(1);
		}

		String splitFileName = argv[0];
		String splitFileName2 = argv[1];
		String controlFileName = argv[2];
		String treatedFileName = argv[3];
		String outFileName = argv[4];

		// first read split probe information get coordinate pairs in a hashtable
		BufferedReader splitReader = new BufferedReader(new FileReader(splitFileName));

		Hashtable probeCoords = new Hashtable(20000);
		Hashtable probeCoordsInv = new Hashtable(20000);

		Hashtable split1 = new Hashtable(20000);

		String splitLine = splitReader.readLine();

		while (splitLine!=null)
		{
			String[] strs = split(splitLine,'\t',true);
			String id = strs[1]+"::"+strs[2];

			split1.put(id,strs);

			int numValid = 	(new Integer(strs[5])).intValue();
			int numInvalid = (new Integer(strs[6])).intValue();
		    for (int i = 0; i<numValid; i++)
		    {
		    	if (!probeCoords.containsKey(strs[8+i*2]))
		    	{
		    		probeCoords.put(strs[8+i*2],id+"#"+i);
		    		probeCoordsInv.put(id+"#"+i,strs[8+i*2]);

		    		//System.out.println(strs[8+i*2]);
		    	}
		    }

		    for (int i = 0; i<numInvalid; i++)
		    {
		    	if (!probeCoords.containsKey(strs[8+2*numValid+i*2]))
		    	{
		    		probeCoords.put(strs[8+2*numValid+i*2],id+"#"+(numValid+i));
		    		probeCoordsInv.put(id+"#"+(numValid+i),strs[8+2*numValid+i*2]);
		    	}
		    }

			splitLine = splitReader.readLine();
		}

		splitReader.close();

		System.out.println("Number of unique probes: "+probeCoords.size());


		// first read split probe information get coordinate pairs in a hashtable
		BufferedReader splitReader2 = new BufferedReader(new FileReader(splitFileName2));

		Hashtable probeCoords2 = new Hashtable(20000);

		Hashtable split2 = new Hashtable(20000);

		String splitLine2 = splitReader2.readLine();

		while (splitLine2!=null)
		{
			String[] strs = split(splitLine2,'\t',true);

			String id = strs[1]+"::"+strs[2];

			split2.put(id,strs);

			int numValid = 	(new Integer(strs[5])).intValue();
			int numInvalid = (new Integer(strs[6])).intValue();
		    for (int i = 0; i<numValid; i++)
		    {
		    	if (!probeCoords2.containsKey(strs[8+i*2]))
		    	{
		    		probeCoords2.put(strs[8+i*2],id+"#"+i);
		    		//System.out.println(strs[8+i*2]);
		    	}
		    }

		    for (int i = 0; i<numInvalid; i++)
		    {
		    	if (!probeCoords2.containsKey(strs[8+2*numValid+i*2]))
		    		probeCoords2.put(strs[8+2*numValid+i*2],id+"#"+(i+numValid));
		    }

			splitLine2 = splitReader2.readLine();
		}

		splitReader2.close();

		System.out.println("Number of unique probes: "+probeCoords2.size());


		// now for each probe we want to get the average probe intensity from the control sample and the average probe intensity from the treated samples
		// maintain two hashtables treatedAverage and controlAverage for this purpose

		// read control samples first

		BufferedReader controlListReader = new BufferedReader(new FileReader(controlFileName));

		Hashtable controlValues = new Hashtable(20000);

		String controlListLine = controlListReader.readLine();

		while (controlListLine!=null)
		{
			BufferedReader controlReader = null;
			try {
				controlReader = new BufferedReader(new FileReader(controlListLine+".CEL"));
			}
			catch (Exception e)
			{
				int k;
				System.out.println("Could not find "+controlListLine+".CEL Downloading from NCBI GEO. Unzip the downloaded files and run the program again");
				//URL link = new URL("http://www.ncbi.nlm.nih.gov/geo/download/?acc="+controlListLine+"&format=file&file="+controlListLine+"%2ECEL%2Egz");
				URL link = new URL("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/"+controlListLine.substring(0,controlListLine.length()-3)+"nnn/"+controlListLine+"suppl/"+controlListLine+"%2ECEL%2Egz");
				DataInputStream fin = null;
				try {
					        fin = new DataInputStream(link.openStream());
				}
				catch (Exception e2)
				{
					fin = null;
					System.out.println("Unable to download "+controlListLine+".CEL. Check your network connection.");
				}

		        if (fin!=null)
		        {
		        	int i;
		        	FileOutputStream outFile=new FileOutputStream(controlListLine+".CEL.gz");
					DataOutputStream outData=new DataOutputStream(outFile);
    			    while ((i = fin.read()) != -1){
	   			        outData.writeByte(i);
	    	    	}
	    	    	outData.close();
			        fin.close();
		        }
	
		        controlListLine = controlListReader.readLine();
		        continue;
			}
			System.out.println("Reading CEL file: "+controlListLine);
			String controlLine = controlReader.readLine();
			int type  = 1; // 1: for U133A, 2: for U133 Plus 2.0
			while (controlLine!=null)
			{
				String[] strs = split(controlLine,'\t',true);
				if (strs.length!=5)
				{
					if (controlLine.startsWith("[INTENSITY]"))
					{
						controlLine = controlReader.readLine();
						int numcells = new Integer(controlLine.substring(12,controlLine.length())).intValue();
						//System.out.println(numcells);
						if (numcells<1000000) type = 1; else type = 2;
					}
					controlLine = controlReader.readLine();
					continue;
				}
				String coords = "("+strs[0].trim()+":"+strs[1].trim()+")";

				if (type == 2) {
					if (probeCoords2.containsKey(coords))
					{
						String probeid = (String)(probeCoords2.get(coords));
						if (probeCoordsInv.containsKey(probeid))
						{
							//System.out.print(probeid+"before: "+coords+" after: ");
							coords=(String)(probeCoordsInv.get(probeid));
							//System.out.println(coords);
						}
						else
						{
							controlLine = controlReader.readLine();
							continue;
						}
					}
					else
					{
						controlLine = controlReader.readLine();
						continue;
					}
				}

				if (probeCoords.containsKey(coords))
				{
					if (controlValues.containsKey(coords))
					{
						Vector values = (Vector)controlValues.get(coords);
						values.addElement(new Float(strs[2].trim()));
						controlValues.remove(coords);
						controlValues.put(coords,values);
					}
					else
					{
						Vector values = new Vector(0,1);
						values.addElement(new Float(strs[2].trim()));
						controlValues.put(coords,values);
					}
				}
				controlLine = controlReader.readLine();
			}
			controlReader.close();
			controlListLine = controlListReader.readLine();
		}
		controlListReader.close();
		// Now go over the controlValues hashtable and process the vector of values to get the average intensities (Float)

		Hashtable controlAverage = new Hashtable(20000);

		for (Enumeration en = controlValues.keys(); en.hasMoreElements(); )
		{
			String id = (String)en.nextElement();

			Vector values = (Vector)controlValues.get(id);

			float total = 0.0f;

			for (int i=0;i<values.size();i++)
			{
				float f = ((Float)values.elementAt(i)).floatValue();
				total+=f;
			}
			total = total/values.size();

			Float value = new Float(total);

			controlAverage.put(id,value);
		}

		// now read treated samples

		BufferedReader treatedListReader = new BufferedReader(new FileReader(treatedFileName));

		Hashtable treatedValues = new Hashtable(20000);

		String treatedListLine = treatedListReader.readLine();

		while (treatedListLine!=null)
		{
			BufferedReader treatedReader = null;
			try {
				treatedReader = new BufferedReader(new FileReader(treatedListLine+".CEL"));
			}
			catch (Exception e)
			{
				System.out.println("Could not find "+treatedListLine+".CEL Downloading from NCBI GEO. Unzip the downloaded files and run the program again");
				URL link = new URL("http://www.ncbi.nlm.nih.gov/geo/download/?acc="+treatedListLine+"&format=file&file="+treatedListLine+"%2ECEL%2Egz");

				DataInputStream fin = null;
				try {
					        fin = new DataInputStream(link.openStream());
				}
				catch (Exception e2)
				{
					fin = null;
					System.out.println("Unable to download "+treatedListLine+".CEL. Check your network connection.");
				}

		        if (fin!=null)
		        {
		        	int i;
		        	FileOutputStream outFile=new FileOutputStream(treatedListLine+".CEL.gz");
					DataOutputStream outData=new DataOutputStream(outFile);
    			    while ((i = fin.read()) != -1){
	   			        outData.writeByte(i);
	    	    	}
	    	    	outData.close();
		        }
		        fin.close();
		        treatedListLine = treatedListReader.readLine();
		        continue;
			}


			System.out.println("Reading CEL file: "+treatedListLine);
			String treatedLine = treatedReader.readLine();
			int type = 1;
			int c=0;

			while (treatedLine!=null)
			{
				String[] strs = split(treatedLine,'\t',true);
				if (strs.length!=5)
				{
					if (treatedLine.startsWith("[INTENSITY]"))
					{
						treatedLine = treatedReader.readLine();
						int numcells = new Integer(treatedLine.substring(12,treatedLine.length())).intValue();
						//System.out.println(numcells);
						if (numcells<1000000) type = 1; else type = 2;
					}

					treatedLine = treatedReader.readLine();
					continue;
				}
				String coords = "("+strs[0].trim()+":"+strs[1].trim()+")";

				if (type == 2) {
					if (probeCoords2.containsKey(coords))
					{
						String probeid = (String)(probeCoords2.get(coords));
						if (probeCoordsInv.containsKey(probeid))
						{
							//System.out.print(probeid+"before: "+coords+" after: ");
							coords=(String)(probeCoordsInv.get(probeid));
							//System.out.println(coords);
						}
						else
						{
							treatedLine = treatedReader.readLine();
							continue;
						}
					}
					else
					{
						treatedLine = treatedReader.readLine();
						continue;
					}
				}

				if (probeCoords.containsKey(coords))
				{
					if (treatedValues.containsKey(coords))
					{
						Vector values = (Vector)treatedValues.get(coords);
						values.addElement(new Float(strs[2].trim()));
						treatedValues.remove(coords);
						treatedValues.put(coords,values);
						c++;
					}
					else
					{
						Vector values = new Vector(0,1);
						values.addElement(new Float(strs[2].trim()));
						treatedValues.put(coords,values);
						c++;
					}
				}
				treatedLine = treatedReader.readLine();
			}
			treatedReader.close();
			treatedListLine = treatedListReader.readLine();
		}
		treatedListReader.close();
		// Now go over the treatedValues hashtable and process the vector of values to get the average intensities (Float)

		Hashtable treatedAverage = new Hashtable(20000);

		for (Enumeration en = treatedValues.keys(); en.hasMoreElements(); )
		{
			String id = (String)en.nextElement();

			Vector values = (Vector)treatedValues.get(id);

			float total = 0.0f;

			for (int i=0;i<values.size();i++)
			{
				float f = ((Float)values.elementAt(i)).floatValue();
				total+=f;
			}
			total = total/values.size();

			Float value = new Float(total);

			treatedAverage.put(id,value);
		}

		// for debugging
		//System.out.println(controlAverage.size());
		//System.out.println(treatedAverage.size());

		// Now go over the split file once more and replace the probe coordinates with two sets of average probe instensities (control and treated)

		FileOutputStream outFile=new FileOutputStream(outFileName);
		DataOutputStream outData=new DataOutputStream(outFile);


		splitReader = new BufferedReader(new FileReader(splitFileName));

		splitLine = splitReader.readLine();

		while (splitLine!=null)
		{
			String[] strs = split(splitLine,'\t',true);

			if (!split2.containsKey(strs[1]+"::"+strs[2])){
				splitLine = splitReader.readLine();
				System.out.println("U133A contains something that is not found in U133 Plus 2.0 "+strs[1]+" "+strs[2]);
				continue;
			}

			String[] strs2 = (String[])(split2.get(strs[1]+"::"+strs[2]));

			if (strs.length!=strs2.length)
				System.out.println("Mismatching number of probes in U133A and U133Plus 2.0: "+strs[1]+"\t"+strs.length+"\t"+strs2.length);


			String id = strs[1]+"::"+strs[2];

			int numValid = 	(new Integer(strs[5])).intValue();
			int numInvalid = (new Integer(strs[6])).intValue();

		    outData.writeBytes(strs[0]+"\t"+strs[1]+"\t"+strs[2]+"\t"+strs[3]+"\t"+strs[4]+"\t"+strs[5]+"\t"+strs[6]);

		    for (int i = 0; i<numValid; i++)
		    {
		    	// 19.1.13 if (controlAverage.containsKey(strs[8+i*2]) && treatedAverage.containsKey(strs[8+i*2]))
		    	if (controlValues.containsKey(strs[8+i*2]) && treatedValues.containsKey(strs[8+i*2])) // 19.1.13
		    	{
		    		// 19.1.13 float caverage = ((Float)controlAverage.get(strs[8+i*2])).floatValue();
		    		// 19.1.13 float taverage = ((Float)treatedAverage.get(strs[8+i*2])).floatValue();
					Vector cvalues = (Vector)controlValues.get(strs[8+i*2]); // 19.1.13
					Vector tvalues = (Vector)treatedValues.get(strs[8+i*2]); // 19.1.13

					//if (cvalues.size()!=74) System.out.println("Number of control samples = "+cvalues.size()+"...."+id);
					//if (tvalues.size()!=520) System.out.println("Number of treated samples = "+tvalues.size()+"...."+id);

		    		DecimalFormat df = new DecimalFormat("#.00");
		    		// 19.1.13 outData.writeBytes("\t"+strs[7+i*2]+"\t"+"("+df.format(caverage)+":"+df.format(taverage)+")");
					outData.writeBytes("\t"+strs[7+i*2]+"\t"+"("); // 19.1.13

					outData.writeBytes(df.format(((Float)cvalues.elementAt(0)).floatValue())); // 19.1.13
					for (int j=1;j<cvalues.size();j++) // 19.1.13
					{ // 19.1.13
						float f = ((Float)cvalues.elementAt(j)).floatValue(); // 19.1.13
						outData.writeBytes("$"+df.format(f)); // 19.1.13
					} // 19.1.13

					outData.writeBytes(":"); // 19.1.13

					outData.writeBytes(df.format(((Float)tvalues.elementAt(0)).floatValue())); // 19.1.13
					for (int j=1;j<tvalues.size();j++) // 19.1.13
					{ // 19.1.13
						float f = ((Float)tvalues.elementAt(j)).floatValue(); // 19.1.13
						outData.writeBytes("$"+df.format(f)); // 19.1.13
					} // 19.1.13

					outData.writeBytes(")"); // 19.1.13
		    	}
		    	else
		    	{
		    		System.err.println("Probe "+id+" not found in the CEL files! Exiting program.");
		    		System.exit(1);
		    	}
		    }

		    for (int i = 0; i<numInvalid; i++)
		    {
		    	// 19.1.13 if (controlAverage.containsKey(strs[8+2*numValid+2*i]) && treatedAverage.containsKey(strs[8+2*numValid+2*i]))
		    	if (controlValues.containsKey(strs[8+2*numValid+2*i]) && treatedValues.containsKey(strs[8+2*numValid+2*i])) // 19.1.13
		    	{
		    		// 19.1.13 float caverage = ((Float)controlAverage.get(strs[8+2*numValid+2*i])).floatValue();
		    		// 19.1.13 float taverage = ((Float)treatedAverage.get(strs[8+2*numValid+2*i])).floatValue();
					Vector cvalues = (Vector)controlValues.get(strs[8+2*numValid+2*i]); // 19.1.13
					Vector tvalues = (Vector)treatedValues.get(strs[8+2*numValid+2*i]); // 19.1.13

					//if (cvalues.size()!=74) System.out.println("Number of control samples = "+cvalues.size()+"...."+id);
					//if (tvalues.size()!=520) System.out.println("Number of treated samples = "+tvalues.size()+"...."+id);

		    		DecimalFormat df = new DecimalFormat("#.00");
		    		// 19.1.13 outData.writeBytes("\t"+strs[7+2*numValid+2*i]+"\t"+"("+df.format(caverage)+":"+df.format(taverage)+")");
					outData.writeBytes("\t"+strs[7+2*numValid+2*i]+"\t"+"("); // 19.1.13

					outData.writeBytes(df.format(((Float)cvalues.elementAt(0)).floatValue())); // 19.1.13
					for (int j=1;j<cvalues.size();j++) // 19.1.13
					{ // 19.1.13
						float f = ((Float)cvalues.elementAt(j)).floatValue(); // 19.1.13
						outData.writeBytes("$"+df.format(f)); // 19.1.13
					} // 19.1.13

					outData.writeBytes(":"); // 19.1.13

					outData.writeBytes(df.format(((Float)tvalues.elementAt(0)).floatValue())); // 19.1.13
					for (int j=1;j<tvalues.size();j++) // 19.1.13
					{ // 19.1.13
						float f = ((Float)tvalues.elementAt(j)).floatValue(); // 19.1.13
						outData.writeBytes("$"+df.format(f)); // 19.1.13
					} // 19.1.13

					outData.writeBytes(")"); // 19.1.13
		    	}
		    	else
		    	{
		    		System.err.println("Probe "+id+" not found in the CEL files! Exiting program.");
		    		System.exit(1);
		    	}
		    }

		    outData.writeBytes("\r\n");

			splitLine = splitReader.readLine();
		}

		splitReader.close();

    	outData.close();

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

class polyAInfo
{
	public String siteName;
	public long position;
	public int numSupport;
	public String strand;
	public polyAInfo(String n, long p, int s, String st)
	{
		siteName = n;
		position = p;
		numSupport = s;
		strand = st;
	}
}

class SplitProbeInfo
{
	public int startLoc;
	public int x;
	public int y;

	public SplitProbeInfo(int a, int b, int c)
	{
		startLoc = a;
		x = b;
		y = c;
	}
}
