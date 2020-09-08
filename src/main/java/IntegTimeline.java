import java.awt.geom.Point2D; 
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.UUID;

import cern.colt.list.DoubleArrayList;
import jsky.coords.WorldCoords;
import nom.tam.fits.BinaryTableHDU;
import nom.tam.fits.Fits;
import nom.tam.fits.FitsException;
import nom.tam.fits.ImageHDU;
import nom.tam.util.BufferedDataInputStream;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

import gb.esac.integral.IntegralException;
import gb.esac.integral.PointingListWriter;
import gb.esac.integral.ScwList;
import gb.esac.integral.ScwListMaker;
import gb.esac.io.MyFile;

public class IntegTimeline {

    //  Arguments
    private static String gnrlIdxFilename = null;
    private static double ra = 266.4168; // SgrA*
    private static double dec = -29.0078;
    private static Point2D.Double radec = null;
    private static double radius = 12;

    //  Other class variables
    private static Logger logger  = Logger.getLogger(IntegTimeline.class);
    private static File loggerFile;
    private static DecimalFormat twoDigits = new DecimalFormat("0.00");

    // main
    public static void main (String[] args) throws Exception {
		configureLogger();
		handleArgs(args);
		File pointLisFile = getPointLisFile();
		ScwList scwlis = ScwListMaker.makeScwList(ra, dec, distInL, distInB, firstRev, lastRev, pointLisFile, normEffArea);
		getDataFromIndexFile();
		selectScws();
		writeTimelineAsQDP();
    }

    //  configureLogger
    private static void configureLogger() throws IOException {
		String loggerFilename= "logger.config";
		InputStream log = getFileFromJarAsStream(loggerFilename);
		UUID uuid = UUID.randomUUID();
		String tmpdir = System.getProperty("java.io.tmpdir");
		loggerFilename = new String(tmpdir+sep+"logger.config_"+uuid.toString());
		inputStreamToFile(log, loggerFilename);
		(new File(loggerFilename)).deleteOnExit();
	    PropertyConfigurator.configure(loggerFilename);
    }
    
    private static void inputStreamToFile(InputStream io, String fileName) throws IOException {       
		FileOutputStream fos = new FileOutputStream(fileName);
		byte[] buf = new byte[256];
		int read = 0;
		while ((read = io.read(buf)) > 0) {
		    fos.write(buf, 0, read);
		}
		fos.flush();
		fos.close();
    }
    
    // arguments
    private static void handleArgs(String[] args) {
		if (args.length < 3 || args.length > 4 ) {
		    logger.error("Usage: java -jar IntegTimeline.jar gnrlIdxFile ra dec (maxDist)");
		    logger.info("  gnrlIdxFile = full path to GNRL-SCWG-GRP-IDX.fits.gz");
		    logger.info("  ra = RA (deg)");
		    logger.info("  dec = Dec (deg)");
		    logger.info("  maxDist = max angular distance from (ra,dec)");
		    System.exit(-1);
		}
		else {
		    gnrlIdxFilename = args[0];
		    ra = (Double.valueOf(args[1])).doubleValue();
		    dec = (Double.valueOf(args[2])).doubleValue();
		    radec = new Point2D.Double(ra,dec);
		    if ( args.length == 4 ) {
			radius = (Double.valueOf(args[3])).doubleValue();
		    }
		}
		logger.info("Running IntegTimeline");
		logger.info("General index filename = "+gnrlIdxFilename);
		logger.info("RA = "+ra);
		logger.info("Dec = "+dec);
		logger.info("Max Dist = "+radius);
    }

    //  Data
    private static String[] scwNum;
    private static double[] pointing_ra;
    private static double[] pointing_dec;
    private static double[] tstart;
    private static double[] tstop;
    private static double[] telapse;
    private static int[] allRevs;
    private static int lastRevNum;
    private static int nscw;

    private static void getDataFromIndexFile() throws Exception {
		// Open or construct point.lis
		File pointlisFile = new File("point.lis");
		AsciiDataFileReader pointlisReader;
		if ( pointlisFile.exists() ) {
		    pointlisReader = new AsciiDataFileReader(pointlisFile);
		}
		else {
		    logger.warn("File point.lis does not exist");
		    constructPointlis(gnrlIdxFilename);
		    pointlisFile = new MyFile("point.lis");
		}

		//  Get the data from point.lis
		scwNum = (String[]) pointlisFile.getStrCol(0);
		pointing_ra = (double[]) pointlisFile.getDblCol(1);
		pointing_dec = (double[]) pointlisFile.getDblCol(2);
		tstart = (double[]) pointlisFile.getDblCol(3);
		tstop = (double[]) pointlisFile.getDblCol(4);
		telapse = (double[]) pointlisFile.getDblCol(5);
		nscw = scwNum.length;
	 	logger.info("File "+pointlisFile.getName()+" contains "+nscw+" rows");

		//  Make list of revs
		allRevs = new int[scwNum.length];
		for ( int i=0; i < scwNum.length; i++ ) {
		    String rev = scwNum[i].substring(0,4);
		    allRevs[i] = (new Integer(rev)).intValue();
		}
		lastRevNum = allRevs[allRevs.length-1];
    }

    private static double[] scwTimes;
    private static double[] scwDts;
    private static double[] scwExpo;
    private static double[] scwOffAxis;
    private static double[] seasonExpo;
    private static double[] seasonAvgRev;
    private static double totTelapse = 0;
    private static double totEffExpo = 0;
    private static double totExpoPerSeason = 0;
    private static int nSelected=0;
    private static double fracFC=0;
    private static double fracPC=0;

    private static void selectScws() throws FitsException, IOException {
	 	logger.info("Selecting science windows ...");
		logger.info("Observation seasons:");
		WorldCoords sourceCoords = new WorldCoords(radec);
		DoubleArrayList scwTimeIDs = new DoubleArrayList();
		DoubleArrayList scwTimeDt = new DoubleArrayList();
		DoubleArrayList scwEffExpo = new DoubleArrayList();
		DoubleArrayList scwOffAxisAngle = new DoubleArrayList();
		DoubleArrayList expoInSeason = new DoubleArrayList();
		DoubleArrayList revsInSeason = new DoubleArrayList();
		DoubleArrayList avgRevInSeason = new DoubleArrayList();
		int previousRevNum = allRevs[0];
		int lastSelectedRevNum = Integer.MAX_VALUE;
		boolean previousRevWasSelected = false;
		int nFC=0;
		int i=0;
		double k=1;
		float[][] normEffArea = getNormIsgriEffArea();

		while ( i < nscw ) {

		    int revNum = allRevs[i];
		    if ( revNum == previousRevNum ) {
			k++;  // increment position within the revolution
		    }
		    else {
			k=1;  // reset position to start of rev
		    }
		    
		    //  Calculate distance from pointing axis
		    Point2D.Double pointing_radec = new Point2D.Double(pointing_ra[i], pointing_dec[i]);
		    WorldCoords pointingCoords = new WorldCoords(pointing_radec);
		    double dist = sourceCoords.dist(pointingCoords);
		    int distInPix = (new Double(Math.rint(dist/4.35))).intValue(); // each pix has a size of 4.35'
		    double distInDeg = dist/60d; //  divide armin dist by 60 to get degrees
		    
		    //  Select scw within search radius
		    if ( distInDeg <= radius && distInPix < 199 ) {

			nSelected++;
			scwOffAxisAngle.add(distInDeg);
			if ( distInDeg <= 8.3 ) nFC++;

			//  Calculate total and effective exposure per scw and per season
			totTelapse += telapse[i];
			totEffExpo += telapse[i]*normEffArea[199][199 - distInPix];
			totExpoPerSeason += telapse[i]*normEffArea[199][199 - distInPix];
			scwEffExpo.add(telapse[i]*normEffArea[199][199 - distInPix]);

			//  Define the time identifier as the scw's fractional position within the revolution
			int nScwsInRev = getNumOfScwsInThisRev(revNum, allRevs);
			double scwTimeID = revNum + k/nScwsInRev;
			scwTimeIDs.add(scwTimeID);
			scwTimeDt.add(1d/nScwsInRev);
			
			//  Check if we are still in the same observing season (min gap = 30 revs)
			double distFromPrevious = revNum - lastSelectedRevNum;
			double distToEnd = lastRevNum - revNum;
			if ( distFromPrevious < 30 ) {
			    revsInSeason.add(previousRevNum);
			}
			else {
			    revsInSeason.trimToSize();
			    double[] revs = revsInSeason.elements();
			    double first = revs[0];
			    double last = revs[revs.length-1];
			    double middle = (first + last)/2;
			    logger.info((int) first+"  to  "+(int) last+" with "+twoDigits.format(totExpoPerSeason/1e3)+" ks");
			    avgRevInSeason.add(middle);
			    revsInSeason.clear();
			    expoInSeason.add(totExpoPerSeason);
			    totExpoPerSeason = 0;
			}
			lastSelectedRevNum = revNum;
		    }
		    previousRevNum = revNum;
		    i++;
		}
		fracFC = (double)nFC/(double)nSelected;
		fracPC = 1d - fracFC;

		//  Process last group of revs
		revsInSeason.trimToSize();
		double[] revs = revsInSeason.elements();
		double first = revs[0];
		double last = revs[revs.length-1];
		double middle = (first + last)/2;
		logger.info((int) first+"  to  "+(int) last+" with "+twoDigits.format(totExpoPerSeason/1e3)+" ks");
		avgRevInSeason.add(middle);
		revsInSeason.clear();
		expoInSeason.add(totExpoPerSeason);

		scwOffAxisAngle.trimToSize();
		scwEffExpo.trimToSize();
		scwTimeIDs.trimToSize();
		scwTimeDt.trimToSize();
		expoInSeason.trimToSize();
		avgRevInSeason.trimToSize();

		scwTimes = scwTimeIDs.elements();
		scwDts = scwTimeDt.elements();
		scwExpo = scwEffExpo.elements();
		scwOffAxis = scwOffAxisAngle.elements();
		seasonExpo = expoInSeason.elements();
		seasonAvgRev = avgRevInSeason.elements();
    }

    private static float[][] getNormIsgriEffArea() throws FitsException, IOException {

		//  Read map of ISGRI effective area and normalize it
		BufferedDataInputStream effAreaFileAsStream = new BufferedDataInputStream(getFileFromJarAsStream("eff_area.fits.gz"));
		Fits effAreaFits = new Fits(effAreaFileAsStream, true);
		ImageHDU effAreaHDU = (ImageHDU) effAreaFits.getHDU(0);
		float[][] effArea = (float[][]) effAreaHDU.getKernel();
		float max = -Float.MAX_VALUE;
		float min = Float.MAX_VALUE;

		//  Determine min and max values
		for ( int row=0; row < 400; row++ ) {
		    for ( int col=0; col < 400; col++ ) {
			max = Math.max(max, effArea[row][col]);
			min = Math.min(min, effArea[row][col]);
		    }
		}
		max = max - min;

		//  Normalize to min=0 and max=1
		float[][] normEffArea = new float[400][400];
		for ( int row=0; row < 400; row++ ) {
		    for ( int col=0; col < 400; col++ ) {
			normEffArea[row][col] = (effArea[row][col] - min)/max;
		    }
		}
		return normEffArea;
    }

    private static File obsTimelineFile;

    private static void writeTimelineAsQDP() throws IOException {
		DecimalFormat threeDigits = new DecimalFormat("0.000");

		String obsTimelineFilename = "timeline_"+ ra +"_"+ dec +"_"+radius+"deg.qdp";
		String obsTimelineFilenamePlot = "timeline_"+ ra +"_"+ dec +"_"+radius+"deg.ps";	
		obsTimelineFile = new File(obsTimelineFilename);

		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(obsTimelineFile)));
		String[] header = new String[] {
		    "! QDP data file",
		    "! Author: Guillaume Belanger - ESA/ESAC",
		    "!",
		    "DEV /XS",
		    "READ SERR 1",
		    "PLOT VERT",
		    "LAB F",
		    "TIME OFF",
		    "LW 3", 
		    "CS 1.0",
		    "LAB OT INTEGRAL Observation Timeline", 
		    "LAB T RA = "+ra+", Dec = "+dec+", Radius = "+radius,
		    //"LINE STEP ON",
		    "MA 1 ON",
		    "CO 2 ON 3",
		    //"MA 39 ON", "MA SIZE 5",
		    "LAB X Revolution Number",
		    "LAB Y2 Effective Exp (ks)",
		    "LAB Y3 Angular Dist (deg)",
		    "R Y2 -2.5 9.5",
		    "R Y3 -3 18",
		    "R X "+(-0.03*lastRevNum)+" "+(1.02*lastRevNum),
		    "VIEW 0.1 0.2 0.9 0.8",
	 	    "WIN 3",
		    "LAB 100 POS "+(-0.02*lastRevNum)+" 8.3 LINE 0 0.99 \"",
		    "LAB 100 LS 4",
		    "LAB 101 POS 0 9.5 \""+(int)Math.round(fracPC*100)+"%\" CS 0.55 JUST CEN",
		    "LAB 102 POS 0 7.1 \""+(int)Math.round(fracFC*100)+"%\" CS 0.55 JUST CEN",
	 	    "LOC 0 0.1 1 0.6",
	 	    "WIN 2",
	 	    "LOC 0 0.42 1 0.92",
		    "LAB 200 POS "+lastRevNum+" 8.25 \"Selected Observations = "+nSelected+"\" CS 0.55 JUST RIGHT",
		    "LAB 201 POS "+lastRevNum+" 7.5 \"Effective Exposure = "+twoDigits.format(totEffExpo/1e6)+" Ms\" CS 0.55 JUST RIGHT",
		    "LAB 202 POS 0 -1.25 \"(ks)\" CS 0.55 JUST CEN"
		};
		for ( int i=0; i < header.length; i++ ) {
		    pw.println(header[i]);
		}
		for ( int i=0; i < seasonExpo.length; i++ ) {
		    pw.println("LAB "+(i+1)+" POS "+twoDigits.format(seasonAvgRev[i])+" -1.25 \""+(int)Math.round(seasonExpo[i]/1e3)+"\" CS 0.55 JUST CEN");
		}
		pw.println("!");
		for ( int i=0; i < scwTimes.length; i++ ) {
		    pw.println(threeDigits.format(scwTimes[i]) +"\t"+ (scwDts[i]/2) +"\t"+ (scwExpo[i]/1e3)  +"\t"+ scwOffAxis[i]);
		}
		pw.println("HARD "+obsTimelineFilenamePlot+"/ps");
		pw.close();
		printStats();
    }

    private static void printStats() {
		logger.info(""+ nSelected +" scw selected");
		logger.info("Total TELAPSE = "+twoDigits.format(totTelapse/1e3)+" ks");
		logger.info("Total off-axis corrected exposure = "+twoDigits.format(totEffExpo/1e3)+" ks");
		logger.info("Total dead-time corrected (0.82) exposure = "+twoDigits.format(totEffExpo/1e3*0.82)+" ks");
		logger.info("Observation timeline written to: "+obsTimelineFile.getPath());
    }

    private static int getNumOfScwsInThisRev(int revToFind, int[] allRevs) {
		int i=0;
		while ( revToFind > allRevs[i] ) {
		    i++;
		}
		int nScws=0;
		while ( revToFind == allRevs[i] ) {
		    nScws++;
		    i++;
		}
		return nScws;
    }

    private static InputStream getFileFromJarAsStream(String name) {
		return ClassLoader.getSystemResourceAsStream(name);
    }

    private static void inputStreamToFile(InputStream io, String fileName) throws IOException {
		FileOutputStream fos = new FileOutputStream(fileName);
		byte[] buf = new byte[256];
		int read = 0;
		while ((read = io.read(buf)) > 0) {
		    fos.write(buf, 0, read);
		}
		fos.flush();
		fos.close();
    }

    private static void constructPointlis(File idxFile) throws Exception {
		DecimalFormat number = new DecimalFormat("0.000");

		//  Get data from GNRL-SCWG-GRP-IDX.fits.gz
		logger.info("Using general index file = "+idxFile.getCanonicalPath());
		logger.info("Reading the data ...");
		Fits f = openFits(idxFile);
		BinaryTableHDU hdu = (BinaryTableHDU) f.getHDU(1);
		String[] scwid = (String[]) hdu.getColumn("SWID");
		String[] scwType = (String[]) hdu.getColumn("SW_TYPE");
		double[] tStart = (double[]) hdu.getColumn("TSTART");
		double[] tStop = (double[]) hdu.getColumn("TSTOP");
		double[] telapse = (double[]) hdu.getColumn("TELAPSE");
		float[] ra_scx = (float[]) hdu.getColumn("RA_SCX");
		float[] dec_scx = (float[]) hdu.getColumn("DEC_SCX");
		byte[] ibismode = (byte[]) hdu.getColumn("IBISMODE");
		float[] ra_scz = (float[]) hdu.getColumn("RA_SCZ");
		float[] dec_scz = (float[]) hdu.getColumn("DEC_SCZ");
		
		//  Construct point.lis
		logger.info("Constructing list of pointings 'point.lis'");
		PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter("point.lis")));
		for ( int i=0; i < scwid.length; i++ ) {
		    if ( scwType[i].equals("POINTING") && ibismode[i] == 41 )
			pw.println(scwid[i] +"\t"+ ra_scx[i] +"\t"+ dec_scx[i] +"\t"+ tStart[i] +"\t"+ tStop[i] +"\t"+ telapse[i] +"\t"+ ra_scz[i] +"\t"+ dec_scz[i]);
		}
		pw.close();
    }

    private static void constructPointlis(String filename) throws Exception {
		constructPointlis(new File(filename));
    }

    private static Fits openFits(File file) throws Exception {
		boolean isGzipped = MyFile.isGzipped(file);
		BufferedDataInputStream dis = new BufferedDataInputStream(new FileInputStream(file));
		Fits fitsFile = new Fits(dis, isGzipped);
		return fitsFile;
    }

    private static Fits openFits(String filename) throws Exception {
		return openFits(new File(filename));
    }

}
