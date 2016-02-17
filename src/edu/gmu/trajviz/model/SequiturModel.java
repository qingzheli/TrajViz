package edu.gmu.trajviz.model;
/*
 * Author: Qingzhe Li
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.math.BigDecimal;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Observable;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.concurrent.TimeUnit;

import org.slf4j.LoggerFactory;

import edu.gmu.trajviz.gi.GrammarRuleRecord;
import edu.gmu.trajviz.gi.GrammarRules;
import edu.gmu.trajviz.gi.sequitur.SAXMotif;
import edu.gmu.trajviz.gi.sequitur.SAXRule;
import edu.gmu.trajviz.gi.sequitur.SequiturFactory;
import edu.gmu.trajviz.logic.*;
import ch.qos.logback.classic.Level;
import ch.qos.logback.classic.Logger;
import edu.gmu.trajviz.sax.datastructures.SAXRecords;
import edu.gmu.trajviz.timeseries.TSException;
import edu.gmu.trajviz.util.StackTrace;
public class SequiturModel extends Observable {
//	public static double MINLINK = 0.0;
//	public final static double (minLink*2) = 0.0;
	public final static double EVAL_ANOMALY_THRESHOLD = 0.2;
	public final static int EVAL_RESOLUTION = 100;
	
	final static Charset DEFAULT_CHARSET = StandardCharsets.UTF_8;
	public final static String EVALUATION_HEAD = "DataName,MinLink,AlphabetSize,MinBlocks,NCThreshold,RunningTime,AvgDistance,AvgeStdDev,MinInterDistance,SilhouetteCoefficient,TotalRules,TotalDataPoints, TotalSubTrajectories,CoveredPoints, ImmergableRuleCount\n";
//	public static final int ALPHABETSIZE = 50;
//	public static final int CONTINUALBLOCKTHRESHOLD = 10;
	//public static final int paaSize = 10;
	private ArrayList<Integer> status;
	private static final String SPACE = " ";
	private static final String CR = "\n";
	private static final int STEP = 2;
	private static final int DEFAULT_TIME_GAP = 6;//180;
//	private static final int DEFAULT_TIME_GAP = 180;
    private boolean[] isCovered;
    private int breakPoint; // the positions<breakPoint are normal, otherwise are abnormal.
    private int trueAnomalyCount;
    private int falsePositiveCount;
    private int trueNegativeCount;
    private int falseNegativeCount;
    private boolean[] ruleCovered;
    private int totalFP;
    private int totalTP;
    private int totalTN;
    private int totalFN;
//	private static final int NOISYELIMINATIONTHRESHOLD = 5;
	public static int alphabetSize;
	private double minLink;
	private static ArrayList<Integer> groundTruth;
	private int noiseThreshold;
	private static GrammarRules rules;
	public static HashMap<String, ArrayList<String>> allPostions;
	public static ArrayList<GrammarRules> allRules;
	public static TreeMap<String, GrammarRuleRecord> sortedRuleMap;
	public static ArrayList< ArrayList<HashSet<Integer>>> allClusters;
	private ArrayList<HashSet<Integer>> clusters;
	private Cluster cluster;
	private HashMap<String, Cluster> currentClusters;
	private ArrayList<Integer> filter;
	public static ArrayList<ArrayList<Integer>> allFilters;
	private HashMap<Integer,Integer> filterMap;
	HashMap<Integer, Integer> clusterMap;
	private ArrayList<Integer> mapTrimed2Original; 
	private ArrayList<Integer> mapToPreviousR0;
	private ArrayList<Integer> mapToOriginalTS;
	public static ArrayList<ArrayList<Integer>> allMapToOriginalTS;
	public static ArrayList<ArrayList<Integer>> allMapToPreviousR0;
	public static HashMap<String, ArrayList<RuleInterval>> finalIntervals;
	private boolean hasNewCluster = true;
	private int sortedCounter;
	private String dataFileName, fileNameOnly;
	private static double lat_center;
	private double latMax;
	private double latMin;
	private double lonMin;
	private double lonMax;
	public int trajCounter;
	private int coverCount;
	private int immergableRuleCount;
	private int totalSubTrajectory;
	private String[] r0;
	private String[] r0Recover;
	private String[] r0Ori; 
	public static ArrayList<String[]> allR0;
	private static double lon_center;
	//The outer arrayList includes all rules, the inner arrayList includes all route under the same rule
	private static ArrayList<ArrayList<Route>> routes;  
	private static ArrayList<Route> rawRoutes;  
	private static ArrayList<Route> anomalyRoutes;
	public ArrayList<Double> lat;
	public static ArrayList<Double> ncLat = new ArrayList<Double>();
	public static ArrayList<Double> ncLon = new ArrayList<Double>();
	//public ArrayList<Double> paaLat;
	//public ArrayList<Double> paaLon;
	public ArrayList<Double> lon;
	public ArrayList<Double> latOri;
	public ArrayList<Double> lonOri;
	private MotifChartData chartData;
	private ArrayList<ArrayList<RuleInterval>> ruleIntervals;
	private ArrayList<RuleInterval> rawAllIntervals;
	private ArrayList<RuleInterval> anomalyIntervals;
	private ArrayList<RuleInterval> anomalRuleIntervals;
	private Integer iteration;
	private boolean isLastIteration;
	private int realRuleSize;
//	private ArrayList<HashSet<Integer>> mapToOriginRules;
	private double runTime = -1;
	//private static GrammarRules filteredRules;
	@SuppressWarnings("rawtypes")
	public ArrayList<NumerosityReductionMapEntry> trimedTrack;
	// index is the rule# after filtering, Integer value is the actual rule number
//	private ArrayList<Integer> filteredRuleMap = new ArrayList<Integer>(); 
	private ArrayList<String> words;
	public Blocks blocks, eBlocks;
	private int minBlocks;

	public ArrayList<Double> allLatOri;

	public ArrayList<Double> allLonOri; 
	private static Logger consoleLogger;
	  private static Level LOGGING_LEVEL = Level.DEBUG;
	
	  static {
	    //consoleLogger = (Logger) LoggerFactory.getLogger(SequiturModel.class);
	    //consoleLogger.setLevel(LOGGING_LEVEL);
	}
	 /**
	   * The file name setter.
	   * 
	   * @param filename The file name to set.
	   */
	  /*
	  public static GrammarRules getFilteredGrammarRules(){
		  return filteredRules;
	  }
	  */
	  private synchronized void setDataFileName(String filename) {
	    this.dataFileName = filename;
	  }

	  /**
	   * The file name getter.
	   * 
	   * @return current filename.
	   */
	  public synchronized String getDataFileName() {
		return this.dataFileName;
	}
	  public synchronized void setFileNameOnly(String filename) {
		  fileNameOnly = filename;
		   
		  }
	  public synchronized void setDataSource(String filename) {

		    //consoleLogger.info("setting the file " + filename + " as current data source");

		    // action
		    this.setDataFileName(filename);

		    // notify the View
		    this.setChanged();
		    notifyObservers(new SequiturMessage(SequiturMessage.DATA_FNAME, this.getDataFileName()));

		    // this notification tells GUI which file was selected as the data source
		    this.log("set file " + filename + " as current data source");

		  }
	  /**
	   * Load the data which is supposedly in the file which is selected as the data source.
	   * 
	   * @param limitStr the limit of lines to read.
	   */
	  public synchronized void loadData(String limitStr) {
		  if((null == this.dataFileName)	|| this.dataFileName.isEmpty()){
			  //???this.log("unable to load data - no data source select yet");
			  return;
		  }
		  Path path = Paths.get(this.dataFileName);
		  if (!Files.exists(path)){
			  //???this.log("file"+ this.dataFileName + "doesn't exist.");
			  return;
		  }
		  // read the input
		  // init the data array
		  ArrayList<Double> data = new ArrayList<Double>();
		  ArrayList<Double> data1 = new ArrayList<Double>();
		  status = new ArrayList<Integer>();    // taxi loading status
		  ArrayList<Long> timeAsUnixEpoc = new ArrayList<Long>();   //number of seconds since Jan. 1 1970 midnight GMT, if the time is in milliseconds, it will need to be converted to seconds,otherwise it may over Integer's limite. 
		  
		  try{
			  long loadLimit = 0l;
			  if(!(null == limitStr)&&!(limitStr.isEmpty())){
				  loadLimit = Long.parseLong(limitStr);
			  }
				  BufferedReader reader = Files.newBufferedReader(path, DEFAULT_CHARSET);
				  String line = null;
				  long lineCounter = 0;
				  int trajectoryCounter = -1001;
				  while ((line = reader.readLine()) !=null){
					  String[] lineSplit = line.trim().split("\\s+|,");
					  double value = new BigDecimal(lineSplit[0]).doubleValue();
					  
					  double value1 = new BigDecimal(lineSplit[1]).doubleValue();
					  int value2 = Integer.parseInt(lineSplit[2]);
					  long value3 = Long.parseLong(lineSplit[3]);
					 /*
					  if(value2==1000)
						  breakPoint = status.size();
					  */
					//  if (value>=37.7254&&value<=37.8212&&value1>=-122.5432&&value1<=-122.3561)
					  {
					  if((lineCounter<=1)||(Math.abs(value3-timeAsUnixEpoc.get(timeAsUnixEpoc.size()-1))<=DEFAULT_TIME_GAP &&(value3-timeAsUnixEpoc.get(timeAsUnixEpoc.size()-1))!=0))
					  {
						  
					  
						  if((value<=90)&&(value>=-90))
						  {
							  data.add(value);
					  
							  data1.add(value1);
							  status.add(value2);
							 
							  timeAsUnixEpoc.add(value3);
						  }
						  else
						  {
							  data.add(value1);
							  data1.add(value);
							  status.add(value2);
							  timeAsUnixEpoc.add(value3);
						  }
					  }
					  else{
						  data.add((double)trajectoryCounter);  //adding dummy point to split two trajectories
						  data1.add((double)trajectoryCounter);
						  trajectoryCounter--;
						  status.add(-1);
						  timeAsUnixEpoc.add(value3);
						  //following is adding the first point of a new trajectories
						  if((value<=90)&&(value>=-90))
						  {
							  data.add(value);
					  
							  data1.add(value1);
							  status.add(value2);
							  timeAsUnixEpoc.add(value3);
						  }
						  else
						  {
							  data.add(value1);
							  data1.add(value);
							  status.add(value2);
							  timeAsUnixEpoc.add(value3);
						  }
					  }
				lineCounter++;
				  }
				if((loadLimit>0&&(lineCounter>=loadLimit))){
					break;
				}
			  }
				  data.add((double)trajectoryCounter);
				  data1.add((double)trajectoryCounter);
				  trajCounter = 0 - (trajectoryCounter+1000);
				  status.add(-1);
				  timeAsUnixEpoc.add(timeAsUnixEpoc.get(timeAsUnixEpoc.size()-1));
			  reader.close();
		  }
		  catch (Exception e){
			  String stackTrace = StackTrace.toString(e);
			  System.err.println(StackTrace.toString(e));
			  //???this.log("error while trying to read data from " + this.dataFileName + ":\n" + stackTrace);
		  }
		  
		
			  this.allLatOri = new ArrayList<Double>();
			  this.allLonOri = new ArrayList<Double>();
			  
		latMax = Double.valueOf(data.get(0));
		lonMax = Double.valueOf(data1.get(0));
		latMin = Double.valueOf(data.get(0));
		lonMin = Double.valueOf(data1.get(0));
		  for(int i = 0; i<data.size(); i++){
			  double temp_latitude = Double.valueOf(data.get(i));
			  double temp_longitude = Double.valueOf(data1.get(i));
			//  //???//1System.out.println("i = "+i+": "+temp_latitude+","+temp_longitude);
			  this.allLatOri.add(temp_latitude);
			  this.allLonOri.add(temp_longitude);
			  if((temp_latitude>=-90)&&temp_latitude>latMax)
				  
				  {
				   		
				  latMax = temp_latitude;
				  }
			  if((temp_latitude>=-90)&&temp_latitude<latMin)
				  {
				  
				  latMin = temp_latitude;
				  }
			  if(temp_longitude>=-180&&temp_longitude>lonMax)
				  lonMax = temp_longitude;
			  if(temp_longitude>=-180&&temp_longitude<lonMin)
				  lonMin = temp_longitude;
			//test loaded points
		//	  //???//1System.out.println(this.lat.get(i)+", "+this.lon.get(i)+","+data.get(i)+", "+data1.get(i));
		  }
		  data = new ArrayList<>();
		  data1 = new ArrayList<>();
		  lat_center = (latMax+latMin)/2;
		  lon_center = (lonMax+lonMin)/2;
		  //???//1System.out.println("lonMax:  "+lonMax+"       lonMin: "+lonMin);
		  //???//1System.out.println("latMax:  "+latMax+"       latMin: "+latMin);
		  //???//1System.out.println("Number of trajectories: "+trajCounter);
		  //???//consoleLogger.debug("loaded " + this.allLatOri.size() + " points and "+trajCounter+" Trajecoties... ");
		  //???this.log("loaded " + this.allLatOri.size() + " points from " + this.dataFileName);
		  
		  
		  
		  setChanged();
		  notifyObservers(new SequiturMessage(SequiturMessage.TIME_SERIES_MESSAGE, this.allLatOri,this.allLonOri));
		  
	}
	  public static double getLatitudeCenter(){
		  return lat_center;
	  }
	  public static double getLongitudeCenter(){
		  return lon_center;
	  }
	  @SuppressWarnings("rawtypes")
	
	  
	  
	  
	  public synchronized void processData(double minLink, int alphabetSize, int minBlocks, int noiseThreshold)throws IOException{
		  sortedCounter = 0;
		  this.isLastIteration = false;
		  this.minLink = minLink;
		  this.minBlocks = minBlocks;
		  this.noiseThreshold = noiseThreshold;
		  this.alphabetSize = alphabetSize;
		  
		  totalTP = 0;
		  totalFP = 0;
		  totalFN = 0;
		  totalTN = 0;
		  for(int ith = 0;ith<1000;ith++)
		  {
			  latOri = new ArrayList<Double>();
			  lonOri = new ArrayList<Double>();
			  for(int jth = ith*250*17; jth<(ith+1)*250*17; jth++)
			  {
				  latOri.add(allLatOri.get(jth));
				  lonOri.add(allLonOri.get(jth));
			  }
			  this.allRules = new ArrayList<GrammarRules>();
			  this.allFilters = new ArrayList<ArrayList<Integer>>();
			  this.allClusters = new ArrayList<ArrayList<HashSet<Integer>>>();
			  this.allR0 = new ArrayList<String[]>();
			  this.allMapToPreviousR0 = new ArrayList<ArrayList<Integer>>();
			  this.allMapToOriginalTS = new ArrayList<ArrayList<Integer>>();
			  this.rawRoutes = new ArrayList<Route>();
			  this.anomalyRoutes = new ArrayList<Route>();
			  this.lat = new ArrayList<Double>();
			  this.lon = new ArrayList<Double>();
			  this.groundTruth = new ArrayList<Integer>();
			  this.currentClusters = new HashMap<String,Cluster>();
			  sortedCounter = 0;
		  
		  
		  
		 // this.lat = new ArrayList<Double>();
		 // this.lon = new ArrayList<Double>();
		  Comparator<String> expandedRuleComparator = new Comparator<String>(){
			  @Override public int compare(String r1, String r2)
			  {
				  Integer iteration1 = 0;
				  Integer iteration2 = 0;
				  Integer rule1 = 0;
				  Integer rule2 = 0;
				  if (r1.charAt(0)=='I' )
					{
						if(r1.contains("r")){
							int rIndex = r1.indexOf("r");
							iteration1 = Integer.valueOf(r1.substring(1, rIndex));
							rule1 = Integer.valueOf(r1.substring(rIndex+1));
						//	//1System.out.println("r1: "+r1+" iteration: "+iteration1+" rule1: "+rule1);
							
					//		//1System.out.println(s+" = "+subRule );
						}
						else 
							throw new IllegalArgumentException(r1+" is not comparable with "+ r2);
					}
				  else
					 throw new IllegalArgumentException(r1+" is not comparable with "+ r2);
				  if (r2.charAt(0)=='I' )
					{
						if(r2.contains("r")){
							int rIndex = r2.indexOf("r");
							iteration2 = Integer.valueOf(r2.substring(1, rIndex));
							rule2 = Integer.valueOf(r2.substring(rIndex+1));
						//	//1System.out.println("r2: "+r2+" iteration2: "+iteration2+" rule2: "+rule2);
							
					//		//1System.out.println(s+" = "+subRule );
							
						}
						else 
							throw new IllegalArgumentException(r1+" is not comparable with "+ r2);
					}
				  else
					  new IllegalArgumentException(r1+" is not comparable with "+ r2);
			
				  if(allRules.get(iteration2).get(rule2).getActualRuleYield() > allRules.get(iteration1).get(rule1).getActualRuleYield())
						  return 1;
				 
				  else 
					  return -1;
				   
			  }
		  };

		  SequiturModel.sortedRuleMap = new TreeMap<String, GrammarRuleRecord>(expandedRuleComparator);
		  
		  hasNewCluster = true;
		  StringBuffer sb = new StringBuffer();
		  if (null == this.latOri ||null == this.lonOri|| this.latOri.size()==0 || this.lonOri.size()==0 ){
			  this.log("unable to \"Process data\" - no data were loaded...");
		  }
		  else{
			  //consoleLogger.info("setting up GI with params: ");
			  sb.append(" algorithm: Sequitur");
			  sb.append(" MinLink: ").append(minLink);
			  sb.append(" Alphabet size: ").append(alphabetSize);
			  sb.append(" Minimal Continuous Blocks: ").append(minBlocks);
			  sb.append(" Noise Cancellation Threshold: ").append(noiseThreshold);
			  //consoleLogger.info(sb.toString());
			 
			 
			  this.log(sb.toString());
		  }
		  rawAllIntervals = new ArrayList<RuleInterval>();
		  long beginTime = System.currentTimeMillis();
			 // beginTime  =  System.nanoTime();
			  //1System.out.println("begin time: "+beginTime);
		  buildModel();
		 
		  /*
		  for (int i = 0; i<rawAllIntervals.size();i++)
			  //1System.out.println("Trajectory "+i+":" + rawAllIntervals.get(i));
			  */
		  drawRawTrajectories();

		  
		  
          allMapToPreviousR0.add(mapToPreviousR0);

		//  runSequitur();
		
		  
          /*
           * 
           *    replace rules with rules' ids and clusters' ids
           * 
           */
		  
		  
          iteration = 0;
          
        /*  //1System.out.println("before:");
          for (int d = 0; d<words.size(); d++)
        	  //1System.out.print(words.get(d)+ " ");
          //1System.out.println();
          */
          
          /*
           * run the algorithm
           */
         
          while(hasNewCluster){
	  		  int lastIteration = iteration;
	  		  hasNewCluster = false;
	  		
        	  //1System.out.println("Iteration: "+iteration);
        	  
        //  if(hasNewCluster)
        	  
        	  runSequitur(iteration);
        	  iteration = iteration + 1;
        	 // this.minLink = this.minLink*2;
        	  //this.minLink = minLink*(iteration+1);
        	  this.isLastIteration = true;
        	 //  drawOnMap();
           	//1System.out.println("total anomalies: "+anomalyRoutes.size());

      }
          /*
          this.isLastIteration = true;
          runSequitur(iteration);
          drawOnMap();
          */
          drawOnMap();
         //	//1System.out.println("total anomalies: "+anomalyRoutes.size());
	  
	  //end while
		  
         
		  //1System.out.println("Sorted Map.size = "+ sortedRuleMap.size()+ "sortedCounter = "+sortedCounter);
		  /*
		  for (int i = 0 ; i<r0.length;i++)
			  //1System.out.println(i+ " : "+r0[i]);
		  while(sortedRuleMap.size()>0)
		  {
			 
			  Entry<String, GrammarRuleRecord> entry = sortedRuleMap.pollFirstEntry();
		//	  //1System.out.println(entry.getKey()+" : "+entry.getValue());
		  }
		  */
	//	  AnomalyDetection();
  
		  
		  
		  
		  this.log("processed data, painting on map");
		  //consoleLogger.info("process finished");
		  setChanged();

		
		 


		  //test
		/*  blocks.printBlockMap();
		  for (int i = 0; i<words.size(); i++)
		  {
			  //1System.out.print("  "+words.get(i));
		  }
		  */
		  /*r
		  //1System.out.println("trackMap:");
		  //1System.out.println(trackMap.toString());
		 // //1System.out.println(map2String(trackMap));
		  //1System.out.println("Postions:\t"+getTrimedPositions(trimedTrack).toString()+"\t");
		  //1System.out.println("TrimedStrs:\t"+getTrimedIds(trimedTrack));
		  */
		

		  //1System.out.println("running time: "+runTime);
		  //1System.out.println("finalInteravals: "+finalIntervals.size());
		  ArrayList<Integer> frequency = new ArrayList<Integer>();
		  	

		  /*
		  for (int i=0;i<ruleIntervals.size();i++)
			  frequency.add(ruleIntervals.get(i).size());
			  */
		  notifyObservers(new SequiturMessage(SequiturMessage.CHART_MESSAGE, this.chartData, ruleIntervals));///, mapToOriginRules));//, frequency ));
		  }
		  
		  System.out.println("FINAL Confusion Matrix:");
			System.out.println("Total True Anomaly:\t"+ totalTP+"\t"+ totalFN);
		  System.out.println("Total False Anomaly:\t"+ totalFP+"\t"+ totalTN);
		 double	errorRate = (double)totalFP/(totalTP+totalFN);
		 double accuracy  = (double)(totalTP+totalTN)/(totalTP+totalFP+totalFN+totalTN);
		 System.out.println("errorRate = "+errorRate);
		 System.out.println("accuracy = "+accuracy);
	  //evaluateResult();
	  }
	  
	  private void buildModel() {
		 routes = new ArrayList< ArrayList<Route>>();
		//  paaLat = new ArrayList<Double>();
		//  paaLon = new ArrayList<Double>();
		//  ArrayList<Double> latBuffer=new ArrayList<Double>();
		//  ArrayList<Double> lonBuffer=new ArrayList<Double>();
		  double avgLat;
		  double avgLon;
		  /*
		   * use the centroid(paaLat,paaLon) to represent the data
		   */
		//  //1System.out.println("paaSize: "+paaSize);
		 /*
		  if(paaSize==1)
		  {
			  paaLat = lat;
			  paaLon = lon;
			  
		  }
		  else
		  */
		   // //1System.out.println("Should not see this msg.");
		  
		  
		  
		  
		  

	//	  //1System.out.println("oriLon" + lon);
		 
		  
		  blocks = new Blocks(alphabetSize,latMin,latMax,lonMin,lonMax);
		  double latCut = blocks.latCut;
		  double lonCut = blocks.lonCut;
		  ncLat = new ArrayList<Double>();
		  ncLon = new ArrayList<Double>();
		  resample(latCut,lonCut);
		//  lat = latOri;
		 // lon = lonOri;
		  
		  isCovered= new boolean[lat.size()];
		  ruleCovered = new boolean[lat.size()];
		  for(int i=0;i<lat.size();i++){
			  isCovered[i] = true;
			 // //1System.out.println(lat.get(i)+" , "+lon.get(i));
			
			  
		  }
		  words = new ArrayList<String>();
		  // add all points into blocks.
		  Integer previousId=(Integer)(-1);
		
		//  HashMap<Integer,Integer> trackMap = new HashMap<Integer, Integer>();
		  trimedTrack = new ArrayList<NumerosityReductionMapEntry>();
		  mapTrimed2Original = new ArrayList<Integer>();  // The index is the position in trimmed array, and the content is position in original time series.
		  mapToOriginalTS = new ArrayList<Integer>();
		  int startPoint = 0;
		  int endPoint = 0;
		  for (int i = 0; i<ncLat.size();i++){
			  Location loc = new Location(ncLat.get(i),ncLon.get(i));
			//  blocks.addPoint2Block(loc); this should not work here because the point will change if it is a noisy point.
			  Integer id = new Integer(blocks.findBlockIdForPoint(loc));
			  
			  if(isNoise(id,i,noiseThreshold)){
				 // lat.set(i, lat.get(i-1));
				  ncLat.set(i, ncLat.get(i-1));
				 // lon.set(i, lon.get(i-1));
				  ncLon.set(i, ncLon.get(i-1));
				  id = previousId;
			  }
			  
			  if(id<-1000){
				  endPoint = i-1;
				  rawAllIntervals.add(new RuleInterval(startPoint, endPoint));
				  startPoint = i+1;
			  }
		
			 
			words.add(id.toString());
			
		//	  //1System.out.println("previousId, id:  "+previousId+",   "+id+"        i:   "+i+"   Lat,Lon: "+paaLat.get(i)+","+paaLon.get(i));
			  Integer trimedIndex = 0;
			  if (!id.equals(previousId))
			  {
				  
				  NumerosityReductionMapEntry<Integer, String> entry = new NumerosityReductionMapEntry<Integer, String>(new Integer(i),id.toString());
		//		  //1System.out.println("entry: "+i+","+id);
				  trimedTrack.add(entry);
				  mapTrimed2Original.add(i);
				  mapToOriginalTS.add(i);
				  //put the new <index,id> pair into a map 
				  //NumerosityReductionMapEntry entry = new NumerosityReductionMapEntry(i,id);
			//	  trackMap.put(i, id);
				  previousId = id;
			  }
			  
			  				  
		  }
		  //1System.out.print("mapTrimed2Original: ");
	//	  printArrayList(mapTrimed2Original);		  
		  /*
		   * Following is put the cleaned location data into block again
		   */
		  /*
		  for(int i=0; i<20;i++)
		  //1System.out.println("Orignal String: " + words.get(i));
		  */
		 // //1System.out.println("StringTrimedTrack:  "+trimedTrack);
	/*
		  for(int i=0; i<trimedTrack.size();i++){
			  //1System.out.println(i+" : "+trimedTrack.get(i).getValue()+" ");
		  }
		  */
		  //1System.out.println();
		  
		 		
	}

	/*
	 * resample latitude and longitude when two points skip blocks.
	 */
	  private void resample(double latCut, double lonCut) {
			int i = 1;
			double latPre = latOri.get(0);
			double lonPre = lonOri.get(0);
			lat.add(latPre);
			lon.add(lonPre);
			ncLat.add(latPre);
			ncLon.add(lonPre);
			groundTruth.add(status.get(0));
			boolean firstPoint = true;
			while(i<latOri.size())
			{
			
				
				
				if(latOri.get(i)<-180)
					{
						lat.add(latOri.get(i));
						lon.add(lonOri.get(i));
						ncLat.add(latOri.get(i));
						ncLon.add(lonOri.get(i));
						groundTruth.add(status.get(i));
						i++;
						firstPoint = true;
					}
				else{
					if(firstPoint){
						lat.add(latOri.get(i));
						lon.add(lonOri.get(i));
						ncLat.add(latOri.get(i));
						ncLon.add(lonOri.get(i));
						groundTruth.add(status.get(i));
						i++;
						firstPoint = false;
					}
					else{
					double latStep = Math.abs(latOri.get(i)-latOri.get(i-1));
					double lonStep = Math.abs(lonOri.get(i)-lonOri.get(i-1));
					
					if(latStep>latCut||lonStep>lonCut){
						int skip = Math.max((int)Math.round(latStep/latCut),(int)Math.round(lonStep/lonCut));
						double latstep = (latOri.get(i)-latOri.get(i-1))/skip;
						double lonstep = (lonOri.get(i)-lonOri.get(i-1))/skip;
						for (int j = 0; j<skip; j++){
							lat.add((latOri.get(i-1)+latstep*(j+1)));
							lon.add((lonOri.get(i-1)+lonstep*(j+1)));
							ncLat.add((latOri.get(i-1)+latstep*(j+1)));
							ncLon.add((lonOri.get(i-1)+lonstep*(j+1)));
							groundTruth.add(status.get(i-1));
							//  //???//1System.out.println(lat.get(i+j)+" , "+lon.get(i+j));

						}
						lat.add(latOri.get(i));
						lon.add(lonOri.get(i));
						ncLat.add(latOri.get(i));
						ncLon.add(lonOri.get(i));
						groundTruth.add(status.get(i));
						i++;
					}
					else
						{
						lat.add(latOri.get(i));
						lon.add(lonOri.get(i));
						ncLat.add(latOri.get(i));
						ncLon.add(lonOri.get(i));
						groundTruth.add(status.get(i));
						i++;
						}
					}
				}
				
			}
			
			
		}

	private void drawRawTrajectories() {
		
		
		
	  			
	  		for (int k=0;k<rawAllIntervals.size();k++)
	  				  {
	  					  Route singleRoute = new Route();
	  					  int startPos = rawAllIntervals.get(k).getStartPos();
	  						int endPos = rawAllIntervals.get(k).getEndPos();
	  						/*
	  						for(int index=startPos; index<=endPos;index++)
	  							isCovered[index]=true;
	  							*/
	  		//				//1System.out.println("startPos: "+startPos);
	  		//				//1System.out.println("endPos: " +endPos);
	  						
	  					//	//1System.out.print("track#: "+counter+":       ");
	  						for (int j = startPos; j<=endPos; j++){
	  							
	  							Location loca = new Location(lat.get(j),lon.get(j));
	  				
	  							singleRoute.addLocation(lat.get(j), lon.get(j));
	  								
	  						
	  							
	  						}
	  						rawRoutes.add(singleRoute);
	  						
	  				  }
	  		//		  //1System.out.println("position size: "+positions.size());
	  			//	  //1System.out.println("route size: "+route.size());
	  				
	  				  //  if(route.size()>2)
	
	  		  	
	}

	private void runSequitur(int iteration) {
			chartData = new MotifChartData(this.dataFileName, lat, lon, 1, alphabetSize); //PAA is always 1.
			  clusters = new ArrayList<HashSet<Integer>>();
			  filter = new ArrayList<Integer>();
			  clusterMap = new HashMap<Integer,Integer>();
			  mapToPreviousR0 = new ArrayList<Integer>();
			  rules = new GrammarRules();
				try{
				  SAXRecords saxFrequencyData = null;
				  saxFrequencyData = SequiturFactory.entries2SAXRecords(trimedTrack);
				  //1System.out.println("Input String Length: " + countSpaces(saxFrequencyData.getSAXString(SPACE)));
				  //consoleLogger.trace("String: " + saxFrequencyData.getSAXString(SPACE));
				//  //1System.out.println("String: "+ saxFrequencyData.getSAXString(SPACE));
				  //consoleLogger.debug("running sequitur...");
				  
				  SAXRule sequiturGrammar = SequiturFactory.runSequitur(saxFrequencyData.getSAXString(SPACE));
				//  //1System.out.println("sequiturGrammar: "+sequiturGrammar.toGrammarRulesData().getRuleRecord(1));
				  //consoleLogger.debug("collecting grammar rules data ...");
				 // GrammarRules rules1 = sequiturGrammar.toGRD();
				 // //1System.out.println("rules size: "+ rules1.size());			 
		          rules = sequiturGrammar.toGrammarRulesData();
		          rules.setParsedString();
		          realRuleSize = rules.size();
		     //     allRules.add(rules);
		          //1System.out.println("real rules size: "+ realRuleSize);
		          //debug
		          
		          
		          //consoleLogger.debug("mapping rule intervals on timeseries ...");
		          GrammarRuleRecord rule0 = rules.get(0);
		          
		          //String rule0 = rules.get(0).getRuleString();
		          int length3 = countSpaces(rule0.getRuleString());
		        		  
		          r0 = rule0.getRuleString().split(" ");
		          r0Recover = rule0.getRuleString().split(" ");
		          //1System.out.println("R0 = "+r0);
		          int length4 = r0.length;
		          if (length3!=length4)
	       		  throw new IndexOutOfBoundsException(length3+":"+length4);
		        /* print all rule details
		         */
		          setR0Occ();
		        //  SequiturFactory.updateRuleIntervals(rules, saxFrequencyData, lat.size());   //Both update intervals and intervals in R0
		          for(int i=0;i<rules.size();i++){
		        	  //1System.out.println("Rule number: "+rules.getRuleRecord(i).getRuleNumber()+" Fre in R0: "+rules.get(i).frequencyInR0()+" LEVEL: "+rules.get(i).getRuleLevel()+" "+rules.get(i)+" StringOccurence: "+rules.getRuleRecord(i).occurrencesToString()+"OccurenceInR0: "+rules.get(i).r0OccurrencesToString()+" Rule String: "+rules.getRuleRecord(i).getExpandedRuleString()+" Rule Positions: "+rules.getRuleRecord(i).getRuleIntervals());
		          }
		        allRules.add(rules);
		       /*  */
		         if(this.isLastIteration)
				  {
		        	 mergeTerminals();
				     clusterRules();
		             replaceBack();
				  }
		       //  if(this.alphabetSize<=100)
		         else
		          {
		        //	 mergeTerminals();
		        	 clusterRules();
		        	 //finalCluster();
		        	 //replaceBack();
		          }
		          
		          
		        
		        	  
		          
		          mapToPreviousR0();
		          allR0.add(r0);
		          
		        	  //1System.out.print("mapToOriginalTS: ");
		        	  ArrayList<Integer> previousMapToOriginalTS = mapToOriginalTS;
		        	 // printArrayList(previousMapToOriginalTS);
		        	  /*new ArrayList<Integer>();
		        	  */
		        	  for (int i = 0; i<mapToOriginalTS.size();i++)
		        	  	{
		        		 // previousMapToOriginalTS.add(mapToOriginalTS.get(i));
		        	  	  //1System.out.print( mapToOriginalTS.get(i) + " ");
		        	  	}
		             //1System.out.println();
		             
		              mapToOriginalTS = new ArrayList<Integer>();
		              for(int i = 0; i<r0.length;i++){
		            	          
		              //1System.out.print(r0[i]+" ");
		              }
		              //1System.out.println();
					  trimedTrack = new ArrayList<NumerosityReductionMapEntry>();
				  /*
				   * Replace Rules' Ids with Clusters' Ids
				   */
					  //1System.out.println("r0.length: "+r0.length);
		          for (int i = 0; i<r0.length;i++){
		        	  NumerosityReductionMapEntry<Integer, String> entry;
		        	  if(r0[i]==null)
		        		  {
		        		  	
		        		 // 	Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
        		  		//	mapToOriginalTS.add(pos);
        	  		//		//1System.out.println("BlockID: " +r0[i]+" : "+pos);//mapTrimed2Original.get(mapToPreviousR0.get(i)));

        		  		//	entry = new NumerosityReductionMapEntry<Integer, String>(pos, null);
        		  			//trimedTrack.add(entry);
		        		//  	continue;
		        		  }
		        	  else{
		        	//  //1System.out.println("r0_"+i+"="+r0[i] );
		        	  if (r0[i].charAt(0)=='R')
		        		  {
		        		  //	if(i==0)
		        		  	//	//1System.out.println("r0[i] = "+r0[i]);
		        		  	Integer ruleNumber = Integer.parseInt(r0[i].substring(1));
		        		  	String currentRule = "I"+iteration+"r"+ruleNumber;
		        		  	sortedRuleMap.put(currentRule, rules.get(ruleNumber));
		        		  //	//1System.out.println("sortedRuleMap.size() = " + sortedRuleMap.size()+" "+currentRule+" : "+rules.get(ruleNumber)+" "+sortedRuleMap);
		        		  	sortedCounter++;
		        		//  	int cursor = rules.get(ruleNumber).getCursor(); 
	
		        		  	if (clusterMap.containsKey(filterMap.get(ruleNumber))){
		        		  			hasNewCluster = true;
		        		  			String s = "I" + (iteration) + "C" + clusterMap.get(filterMap.get(ruleNumber));
		        		  			r0[i] = s;
		        		  			Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
		        		  			mapToOriginalTS.add(pos);
		        	  		//		//1System.out.println("BlockID: " +r0[i]+" : "+pos);//mapTrimed2Original.get(mapToPreviousR0.get(i)));
	
		        		  			entry = new NumerosityReductionMapEntry<Integer, String>(pos, s);
		        		  			trimedTrack.add(entry);
		        		  		
		        		  	}
		        		  	else{
		        		  			String s = "I" + (iteration) + "r"+ruleNumber;
		        		  			r0[i] = s;
		        		  			Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
		        		  			mapToOriginalTS.add(pos);
	
		        		  		//	//1System.out.println("RuleID: " +r0[i]+" : "+pos);
	
		        		  			entry = new NumerosityReductionMapEntry<Integer, String>(pos, s);
		        	  				trimedTrack.add(entry);
			
		        		  	
		        		  	}
		        		  		
		        			
		        		  	
		        		  }
		        	  else
		        		  
		        	  {
				  		
		        		   
		        		  	Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
				  			mapToOriginalTS.add(pos);
	
			  				entry = new NumerosityReductionMapEntry<Integer, String>(pos, r0[i]);
			  			//	//1System.out.println("BlockID: " +r0[i]+" : "+pos);//mapTrimed2Original.get(mapToPreviousR0.get(i)));
			  				trimedTrack.add(entry);
	
		        	  }
		          }
		          }
		          
		          /*
		          //1System.out.println("after:");
		          for (int d = 0; d<words.size(); d++)
		        	  //1System.out.print(words.get(d)+ " ");
		          //1System.out.println();
		          */
		                
		          
		          //1System.out.println();
		          allMapToPreviousR0.add(mapToPreviousR0);
		          
		         
		          
		          //consoleLogger.debug("done ...");
		          
		          
		          
		          
		          
		          chartData.setGrammarRules(rules);
		          //1System.out.println("chartData size: "+ chartData.getRulesNumber());
				
		
			  }
			  catch (TSException e){
				  this.log("error while processing data "+StackTrace.toString(e));
				  e.printStackTrace();
			  }
			//  allMapToOriginalTS.add(mapToOriginalTS);
			  		
		}

	
	 
	  
	  /*
		   * Generate All Motifs and record them on files respectively.
		   */
	/*
		private void drawOnMap() {
			 // Generate All Motifs and record them on files respectively.
			 // String header = "type,latitude,longitude";
	//		  //1System.out.println("Total rules:"+chartData.getRulesNumber());
			  
			//  ArrayList<SAXMotif> allMotifs = chartData.getAllMotifs();
			//  for (int i=1; i<chartData.getRulesNumber();i++){
			    // create merged rule interval data structure corresponding to "clusters" 
			    //ruleIntervals = new ArrayList<ArrayList<RuleInterval>>();
		//        mapToOriginRules = new ArrayList<HashSet<Integer>>();
		        
			    int totalRuleCount = 0;
			    immergableRuleCount = 0;
			    
			    for(int i=0;i<filter.size();i++){
			    	// getRulePositions() was modified to show the intervals only occurred in R0
			  //  	if(!clusterMap.containsKey(i)&&chartData.getRulePositionsByRuleNum(filter.get(i)).size()>=minBlocks) 
			    		{
			    		    ArrayList<RuleInterval> ri = chartData.getRulePositionsByRuleNum(filter.get(i));
			    		  
			    			{ruleIntervals.add(chartData.getRulePositionsByRuleNum(filter.get(i)));
			    			HashSet<Integer> set = new HashSet<Integer>();
			    			set.add(filter.get(i));
			    			mapToOriginRules.add(set);
			    			totalRuleCount++;
			    			immergableRuleCount++;
			    			}
			    			  if(ri.size()<=2)
			    		    	//1System.out.println("Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+ri);
			    		}
			    }
			    
			    	
			    for(int i = 0; i< clusters.size();i++){
			    	
			    	ArrayList<RuleInterval> mergedIntervals = new ArrayList<RuleInterval>();
			    	HashSet<Integer> set = new HashSet<Integer>();
			    	if(clusters.get(i).size()>0)
			    	{
			    		////1System.out.println("cluster "+i+" : {" +clusters.get(i)+"}");
			    		totalRuleCount = totalRuleCount+clusters.get(i).size();
			    	
			    	for(int r : clusters.get(i)){
			    	
			    		int rule = r; //filter.get(r);
			    		set.add(rule);
			    		if(mergedIntervals.size()==0)
			    			mergedIntervals.addAll(chartData.getRulePositionsByRuleNum(rule));
			    		else{
			    			ArrayList<RuleInterval> newIntervals = chartData.getRulePositionsByRuleNum(rule);
			    			for(int j = 0; j<newIntervals.size();j++){
			    				RuleInterval newComer = newIntervals.get(j);
			    				boolean hasMerged = false;
			    				for(int k = 0; k<mergedIntervals.size();k++){
			    					if(RuleInterval.isMergable(mergedIntervals.get(k),newComer)){
			    						//1System.out.println("I still need merge here.");
			    						//1System.out.println("1:" +mergedIntervals.get(k) +" 2:"+ newComer );
			    							
			    						RuleInterval newInterval = RuleInterval.merge(mergedIntervals.get(k), newComer);
			    						mergedIntervals.set(k, newInterval);
			    						hasMerged = true;
			    						//mergedIntervals.remove(k);
			    						//mergedIntervals.add(newInterval);
			    						//break;
			    						
			    					}
			    					
			    						
			    				}
			    				if(!hasMerged)
			    					mergedIntervals.add(newComer);
			    				
			    			}
			    		}
			    		}
			    	if(mergedIntervals.size()>2)
			    	{	
			    	ruleIntervals.add(mergedIntervals);
			    	mapToOriginRules.add(set);
			    	}
			    	else
			    		//1System.out.println("mergedIntervals.size = "+mergedIntervals.size());
			    	}
			    }
			    //1System.out.println("Immergable Rule  = "+ immergableRuleCount);
			    //1System.out.println("Total Rule Count = "+totalRuleCount);
			    boolean[] isCovered = new boolean[lat.size()];
			    coverCount = 0;
			    for (int i = 0; i<isCovered.length;i++)
			    	isCovered[i] = false;
			    totalSubTrajectory = 0;
				for (int i = 0; i<ruleIntervals.size();i++){
			  	totalSubTrajectory = totalSubTrajectory + ruleIntervals.get(i).size();
			  	
				  {//(countSpaces(chartData.getRule(i).getExpandedRuleString())>minBlocks){
				  ArrayList<RuleInterval> positions = ruleIntervals.get(i);//chartData.getRulePositionsByRuleNum(filteredRuleMap.get(i));
				  
				//  ArrayList<RuleInterval> positions = chartData.getRulePositionsByRuleNum(i);	  
			//	  //1System.out.println("rule" + i+" :  "+ positions);//.get(0).toString());
				  
				  
				  if(true)//(positions.size()>2)
					  //&&chartData.getRule(i).getMeanLength()>1)
				  {
					  
					  
						//  File fname = new File("./rules/motif_"+i+".csv");
						//  FileWriter motifPos = new FileWriter(fname);
					//Generating evaluation file
						  
						  int counter = 0;
						  ArrayList<Route> route = new ArrayList<Route>();
						 // Integer route0Id =-1;
						 // Integer route1Id =-1;
					  for (int k=0;k<positions.size();k++)
					  {
						  Route singleRoute = new Route();
				//		  motifPos.append(header+"\n");
						  int startPos = positions.get(k).getStartPos();
							int endPos = positions.get(k).getEndPos();
							
							for(int index=startPos; index<=endPos;index++)
								isCovered[index]=true;
			//				//1System.out.println("startPos: "+startPos);
			//				//1System.out.println("endPos: " +endPos);
							boolean firstPoint = true;
							
						//	//1System.out.print("track#: "+counter+":       ");
							for (int j = startPos; j<=endPos; j++){
								
					//			motifPos.append("T,"+lat.get(j)+","+lon.get(j));
							//	if(counter<2){
								Location loca = new Location(lat.get(j),lon.get(j));
							//	  blocks.addPoint2Block(loc);
						//		  Integer idss = new Integer(blocks.findBlockIdForPoint(loca));
						//		  //1System.out.print(idss+", ");
								singleRoute.addLocation(lat.get(j), lon.get(j));
									
							//	}
								if(firstPoint)
								{
						//			motifPos.append(",Track "+counter+",red\n");
									firstPoint = false;
									
								}
						//		else motifPos.append("\n");
								
							}
							route.add(singleRoute);
							
							counter++;
						//	//1System.out.println();
					  }
			//		  //1System.out.println("position size: "+positions.size());
				//	  //1System.out.println("route size: "+route.size());
					
					  //  if(route.size()>2)
					     routes.add(route);
					    
					    
					    
					    
					  //	motifPos.flush();
					  //	motifPos.close();
					  //	//1System.out.println(fname.getName());
				//	  	route.get(0).print();
					  	
					 
				  }
				  
				//  //1System.out.println("motif index: "+motif.getRuleIndex()+"   " +motif.toString());
				  //FileWriter motifPos = new FileWriter(new File("./motif_"+motif.getRuleIndex()+".csv"));
				 
			  	}
			  }	
				for (int i = 0;i<isCovered.length;i++){
					  if(isCovered[i]==true)
						  coverCount++;
				  }
				  //1System.out.println("Cover Count: "+ coverCount);
				  //1System.out.println("cover rate: " +(double)coverCount/isCovered.length);
			 
		}
	*/
	
		private void finalCluster() {
			//currentClusters = new HashMap<String, Cluster>();
			for(int i = 0; i<r0.length; i++	){
				String s = r0[i];
				if(!isNumeric(s) && s!=null){
					if(r0[i].charAt(0)=='R'){
						s = "I"+iteration+"r"+r0[i].substring(1);
						
						cluster = new Cluster(s);
						//cluster.addRule(rules.get(Integer.valueOf(r0[i].substring(1))));
						cluster.addRule(rules.get(Integer.valueOf(r0[i].substring(1))));
						//r0[i]=s;
					    currentClusters.put(s, cluster);	
					}

				}
			}
			RuleDistanceMatrix rdm = new RuleDistanceMatrix(blocks,currentClusters,minBlocks,minLink);
			
		
	}

		private void mapToPreviousR0() {
			int currentIdx = 0;
		   
		         /*
		          * add on mapToPreviousR0;
		          */
		          for(int i=0;i<r0.length;i++){
		        	  
		        	  if(r0[i]==null)
		        		  {
		        		  throw new IndexOutOfBoundsException(i+" : "+r0[i]+":"+currentIdx+" expandRule:  ");//+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);
		        		  	
		        		  }
		        	 
		        	  else if(r0[i].charAt(0)=='R')
		        		  {
		        		  	Integer currentRule = Integer.valueOf(r0[i].substring(1));
		        		  
		        		  
		        		  	mapToPreviousR0.add(currentIdx);
		        	//	  	//1System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
		        		  	int length1 = rules.get(currentRule).getRuleYield();
		        		  	int length2 = countSpaces(rules.get(currentRule).getExpandedRuleString());
		        		   // currentIdx = currentIdx + rules.get(currentRule).getRuleYield();
		        		  	currentIdx = currentIdx + length2;
		        	//	  	//1System.out.println("CurrentIdx = "+currentIdx +" i= "+i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);
		        		  	
		        		  	if(currentIdx>mapToOriginalTS.size()||length1!=length2)
		        		    	
		        		    {
				        		  throw new IndexOutOfBoundsException(i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);

		        		    }
		        		   
		        		    
		        		  }
		        	  else
		        		  {
		        			  
		        		  mapToPreviousR0.add(currentIdx);

		        	//	  //1System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
		        		  
		        		  	currentIdx++;
		        		  	if(currentIdx>mapToOriginalTS.size())
		        		    	
		        		    {
				        //		  //1System.out.println(i+" : "+r0[i]+":"+currentIdx);

		        		    }
		        		  }
		        		  
		          }
		          /* above
			          * add on mapToPreviousR0;
			          */
		          
		          
		          
		          //1System.out.println();
	          		
	}

		private void replaceBack() {
			//1System.out.println("filter");
			printArrayList(filter);
			//1System.out.println("filterMap: ");
			//1System.out.println(filterMap);
			for(int i = realRuleSize; i<rules.size(); i++){
				String[] ruleString = rules.get(i).getRuleString().split(" ");
				//1System.out.println("Rule "+i+" : "+rules.get(i).getRuleString()+"   filterMap: "+filterMap.get(rules.get(i).getRuleNumber()) );
	        	  if(!clusterMap.containsKey(filterMap.get(rules.get(i).getRuleNumber()))){
	        		
	        		//1System.out.print("replace Back Rule String: [ ");
	        		int r0Pos = rules.get(i).getR0Occurrences().get(0);
	        		for(int j = 0; j<ruleString.length;j++){
	        			r0[r0Pos+j] = ruleString[j];
	        			//1System.out.print(ruleString[j]+" ");
	        		}
	        	    //1System.out.println("]");
	        	  }
	        	
	        		  
	          }
			ArrayList<String> r0new = new ArrayList<String>();
			
			for(int i = 0; i<r0.length; i++)
				if(r0[i]!=null)
					{
						r0new.add(r0[i]);
					}
			r0 = new String[r0new.size()];
			
			//1System.out.print("r0new: [");
			for(int i = 0; i<r0.length; i++)
				
					{ 
						r0[i] = r0new.get(i);
						//1System.out.print(" "+r0[i]);
					}
			//1System.out.println("]");
	}

		private void setR0Occ() {
			HashMap<String, Integer> hm = new HashMap<String, Integer>();
	          
	          
	          
	          
	          
			  
			//  allR0.add(r0);
	          for(int i = 0; i<rules.size();i++){
	        	  String key = rules.get(i).getRuleName();
	        	 // String expandedString = rules.get(i).getExpandedRuleString();
	        	//  //1System.out.println(rules.get(i));
	        	  hm.put(key, 0);
	          }
	          
	         // //1System.out.println("R0: "+rule0.getRuleString());
	          //1System.out.print("r0: ");
	          
	          for(int i = 0; i<r0.length;i++){
	        	          
	          //1System.out.print(r0[i]+" ");
	          }
	          //1System.out.println();
	      //    //1System.out.println(r0);
	          int currentIdx = 0;
	       //   int[] indexes = new int[r0.length];
	          
	          
	          
	          
	          
	         /*
	          * add on mapToPreviousR0;
	          */
	          for(int i=0;i<r0.length;i++){
	        	  
	        	  if(r0[i]==null)
	        		  {
	        	//	  mapToPreviousR0.add(currentIdx);
	        			
	  		        //		  //1System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
	  		        		  
	  		        		 // 	currentIdx++;
	  		        		  	if(currentIdx>mapToOriginalTS.size())
	  		        		    	
	  		        		    {
	  				        		  //1System.out.println(i+" : "+r0[i]+":"+currentIdx);

	  		        		    }
	        		  	
	        		  }
	        	 
	        	  else if(r0[i].charAt(0)=='R')
	        		  {
	        		  	Integer currentRule = Integer.valueOf(r0[i].substring(1));
	        		  	hm.put(r0[i], hm.get(r0[i])+1);
	        		  	rules.get(currentRule).addR0Occurrence(currentIdx); // setOccurenceInR0
	        		  	//mapToPreviousR0.add(currentIdx);
	        		  //	//1System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
	        		  	int length1 = rules.get(currentRule).getRuleYield();
	        		  	int length2 = countSpaces(rules.get(currentRule).getExpandedRuleString());
	        		   // currentIdx = currentIdx + rules.get(currentRule).getRuleYield();
	        		  	currentIdx = currentIdx + length2;
	        		  	//1System.out.println("CurrentIdx = "+currentIdx +" i= "+i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);
	        		  	if(currentIdx>mapToOriginalTS.size())
	        		    	
	        		    {
			        		  throw new IndexOutOfBoundsException(i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);

	        		    }
	        		  }
	        	  else
	        		  {
	        			  
	        	//	  mapToPreviousR0.add(currentIdx);

	        //		  //1System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
	        		  
	        		  	currentIdx++;
	        		  	if(currentIdx>mapToOriginalTS.size())
	        		    	
	        		    {
			        		  //1System.out.println(i+" : "+r0[i]+":"+currentIdx);

	        		    }
	        		  }
	        		  
	          }
	          /* above
		          * add on mapToPreviousR0;
		          */
	          
	          
	          
	          //1System.out.println();
	    //      //1System.out.print("mapToPreviousR0: ");
	       //   printArrayList(mapToPreviousR0);
	          
	          for(int i = 1; i<rules.size();i++){
	        	  
	        	  String key = rules.get(i).getRuleName();
	        	  rules.get(i).setFrequencyInR0((hm.get(key)).intValue());
	          }
	          
	          
	          
	          
	          		
	}

		private void clusterRules() {
			
		      /*
	         * Postprocessing merge, connect
	         */
			
			/*
			 * Warning: rules in clusters are real rules, rules in clusterMap are filter rules. 
			 * 
			 */
		    
		    
	        /* print all rule details
	         */
	         
	          for(int i=0;i<rules.size();i++){
	        	  //1System.out.println("Rule number: "+rules.getRuleRecord(i).getRuleNumber()+" Fre in R0: "+rules.get(i).frequencyInR0()+" LEVEL: "+rules.get(i).getRuleLevel()+" "+rules.get(i)+" StringOccurence: "+rules.getRuleRecord(i).occurrencesToString()+"OccurenceInR0: "+rules.get(i).r0OccurrencesToString()+" Rule String: "+rules.getRuleRecord(i).getExpandedRuleString()+" Rule Positions: "+rules.getRuleRecord(i).getR0Intervals());
	          }
	        
	       /*  */
	         
	        filterMap = new HashMap<Integer,Integer>();
	        for (int i = 1; i<rules.size();i++){
	        	//1System.out.println("Before filter: Frequency in R0: "+ rules.get(i).frequencyInR0()+"  Yield: "+rules.get(i).getRuleYield()+" string: "+rules.get(i).getExpandedRuleString());
					if ((rules.get(i).frequencyInR0()>=1&&countSpaces(RuleDistanceMatrix.parseRule(rules.get(i).getExpandedRuleString()))>=1))//||
						//	(originalRules.get(i).frequencyInR0()>1&&originalRules.get(i).getR0Intervals().size()>2&&originalRules.get(i).getRuleYield()>=minBlocks))
						{
						//HashSet<Integer> set = new HashSet<Integer>();
						//1System.out.println("Yield: "+rules.get(i).getRuleYield()+" string: "+rules.get(i).getExpandedRuleString());
						filterMap.put(i, filter.size());
						filter.add(i);
					/*	
						if(rules.get(i).getR0Intervals().size()<2)
							//1System.out.println("Bug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+i);
						*/
						}
					
				}
	        //1System.out.println("filter Size = "+filter.size());
	        allFilters.add(filter);
	        if(filter.size()>1){
	        //HashMap<Integer,ArrayList<Integer>> mergeRecord = new HashMap<Integer, ArrayList<Integer>>();
	        long t1s = System.currentTimeMillis();
	        RuleDistanceMatrix rdm;
	        //1System.out.println("AlphabetSize="+this.alphabetSize);
	        rdm = new RuleDistanceMatrix(blocks,rules, filter,minBlocks, minLink); 
	        long t1e = System.currentTimeMillis();
	        long buildMatrixTime = t1e-t1s;
	        
	  //      clusters = new ArrayList<HashSet<Integer>>(); 
	       /*
	        for(int i = 0; i<rdm.filter.size();i++){
	      	 families.add(new HashSet<Integer>());
	      	 
	      	 
	      	 families.get(i).add(i);
	      	//  mergeRecord.put(i, family.add(i) );
	        }
	       */
	        long t2s =System.currentTimeMillis();
	        NumberFormat formatter = new DecimalFormat("#0.00");
	        //1System.out.println("rdm.pq.size(): "+rdm.pq.size());
	        int mergableCount = 0;
	        while(rdm.pq.size()>0){
	      	  PairDistance pair = rdm.pq.remove();
	      	  int lineSize;
	      	  int colSize;
	      	  int totalSize;
	      	  if(isMergable(rdm.matrix,clusters,pair.getLine(),pair.getCol(),clusterMap, minLink)){
	      		  hasNewCluster=true;
	      		  mergableCount++;
	      	//	  merge(rules,rdm.filter.get(pair.getLine()),rdm.filter.get(pair.getCol()));
	      		  if(clusterMap.containsKey(pair.getLine())||clusterMap.containsKey(pair.getCol()))
	      		  {
	      			  if(!clusterMap.containsKey(pair.getLine())){
	      				  clusters.get(clusterMap.get(pair.getCol())).add(filter.get(pair.getLine()));
	      				  clusterMap.put(pair.getLine(), clusterMap.get(pair.getCol()));
	      			//	  //1System.out.println("Adding Line  to a cluster, Line:"+pair.getLine()+" Colu:"+pair.getCol()+clusters.get(clusterMap.get(pair.getCol())));
	      				//  //1System.out.println("Map:"+clusterMap);
	      				  
	      			  }
	      			  else if(!clusterMap.containsKey(pair.getCol())){
	      				  clusters.get(clusterMap.get(pair.getLine())).add(filter.get(pair.getCol()));
	      				  clusterMap.put(pair.getCol(), clusterMap.get(pair.getLine()));
	      			//	  //1System.out.println("Adding Colum to a cluster,Colum:"+pair.getCol()+" Colu:"+pair.getCol()+clusters.get(clusterMap.get(pair.getLine())));
	      			//	  //1System.out.println("Map:"+clusterMap);
	      			  }
	      			  else{
	      				  if(!clusterMap.get(pair.getLine()).equals(clusterMap.get(pair.getCol())))
	      				  {
	      				//  //1System.out.println("Before Merge, line in cluster:"+clusterMap.get(pair.getLine())+clusters.get(clusterMap.get(pair.getLine()))+" colu in cluster:"+clusterMap.get(pair.getCol())+clusters.get(clusterMap.get(pair.getCol())));
	      				  lineSize = clusters.get(clusterMap.get(pair.getLine())).size();
	      				  colSize = clusters.get(clusterMap.get(pair.getCol())).size();
	      				  clusters.get(clusterMap.get(pair.getLine())).addAll(clusters.get(clusterMap.get(pair.getCol())));
	      				  int colCluster = clusterMap.get(pair.getCol());
	      				  for(int v : clusters.get(clusterMap.get(pair.getCol())))
	      					  {
	      				//	  //1System.out.print("v: "+v+" ");
	      					  clusterMap.put(filterMap.get(v), clusterMap.get(pair.getLine()));
	      					  clusters.get(clusterMap.get(pair.getLine())).add(v);
	      					  }
	      				  ////1System.out.println();
	      				  clusters.get(colCluster).clear();
	      				 // //1System.out.println("After  Merge, Line:"+pair.getLine()+clusters.get(clusterMap.get(pair.getLine()))+" Colu:"+pair.getCol()+clusters.get(colCluster));
	      				 // //1System.out.println("Map:"+clusterMap);
	      				  totalSize = clusters.get(clusterMap.get(pair.getLine())).size();
	      				  //if((lineSize+colSize)!=totalSize){
	      					//  //1System.out.println("Error Candidate here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	      				  //}
	      				  }
	      				  //else
	      					// //1System.out.println("Same Cluster! "+clusterMap.get(pair.getLine())+","+clusterMap.get(pair.getCol()));
	      			  }
	      		  }
	      		  else{
	      			  HashSet<Integer> set = new HashSet<Integer>();
	      			  set.add(filter.get(pair.getLine()));            
	      			  set.add(filter.get(pair.getCol()));
	      			  clusters.add(set);
	      			  clusterMap.put(pair.getLine(), clusters.size()-1);
	      			  clusterMap.put(pair.getCol(), clusters.size()-1);
	      			//  //1System.out.println("Created a cluster: "+clusters.get(clusters.size()-1));
	      			//  //1System.out.println("Map:"+clusterMap);
	      		  }
	      		  
	      		  /*
	      		  clusters.get(pair.getLine()).addAll(clusters.get(pair.getCol()));
	      		  clusters.get(pair.getCol()).addAll(clusters.get(pair.getLine()));
	      		  
	      		  for(int i: families.get(pair.getCol()))
	      			  families.get(pair.getLine()).add(i);
	      		  for(int i: families.get(pair.getLine()))
	      			  families.get(pair.getCol()).add(i);
	      			  */	        		
	      	//	  //1System.out.print("Merged Pair: <"+pair.getLine()+", "+pair.getCol()+"> = "+rdm.matrix[pair.getLine()][pair.getCol()]);
	      	//	  //1System.out.print(" all distances: ");
	      		  /*
	      		  for (int i : clusters.get(clusterMap.get(pair.getLine())))
	    				for(int j : clusters.get(clusterMap.pair.getCol()))
	    				{
	    				
	    				//1System.out.print(formatter.format(rdm.matrix[i][j])+", ");
	    					
	    				}*/
	      	//	  //1System.out.println();
	      	  }
	        }
	        
		  
	        
	        
	        //1System.out.println("MergableCount: "+mergableCount);
	        
	        /*
	        ArrayList<HashSet<Integer>> tempCluster = new ArrayList<HashSet<Integer>>();
	        for(int i=0;i<clusters.size();i++)
	        {
	      	  if(clusters.get(i).size()>0)
	      		  tempCluster.add(clusters.get(i));
	        }
	        clusters = tempCluster; // be aware!!!! hashMap did not update here, but who cares?
	        */
	        allClusters.add(clusters);
	        long t2e = System.currentTimeMillis();
	        long clusterTime = t2e -t2s;
	    	//1System.out.println("build matrix: "+(double)(buildMatrixTime/1000.0));
			  //1System.out.println("Clustering Time: "+(clusterTime/1000.0));
	        /*
	        for(int i = 0; i<clusters.size();i++){
	      	  //1System.out.println("i = "+i+" : "+clusters.get(i));
	        }
	        */
	       
			  //1System.out.println("cluster map size = "+ clusterMap.size());
			    //1System.out.println("clusterMap:   "+clusterMap);
	        }
	        
		}
private void mergeTerminals() {
	int c = 0;
	ArrayList<String> r0New = new ArrayList<String>(); 
	while(c<r0.length){
		String s = r0[c];
	
		String ruleString;
		String expandedRuleString;
		Integer posR0 = 0;
		int length = 0;
    	//	  //1System.out.println(minBlocks+"  = getNextNonTerminal(i) = "+ i +" =  " +getNextNonTerminal(i)+" = "+r0[i]);
			
	    	if(isNumeric(s)&&Integer.valueOf(s)>=0){
	    	
	    		int numStartPos;
	    		int numEndPos;
	    		if(c<1||isNumeric(r0[c-1])||r0[c-1]==null){
	    			numStartPos= mapToOriginalTS.get(c);
	    			length = Math.min(10, getNextNonTerminal(c)-c);
	    		
	    			StringBuffer sb = new StringBuffer();
	    			for(int i = 0; i< length; i++){
	    				sb.append(r0[c+i]+" ");
	    			}
	    			ruleString = sb.toString();
	    			expandedRuleString = sb.toString();
	    			numEndPos = mapToOriginalTS.get(c+length);
	    		    posR0 = c;
	    		    r0New.add("R"+rules.size());
	    		    r0[c] = "R"+rules.size();
	    			
	    			for (int i = c+1; i<c+length ; i++)
	    				r0[i] = null;
	    			
	    		}
	    		else{
	    			String p = r0[c-1];
	    			numStartPos = mapToOriginalTS.get(c-1);
	    			length = Math.min(10, getNextNonTerminal(c)-c);
	    			StringBuffer sb = new StringBuffer();
	    		
	    			for(int i = 0; i< length; i++){
	    				sb.append(r0[c+i]+" ");
	    			
	    			}
	    			ruleString = r0[c-1] + " "+ sb.toString();
	    			//1System.out.println("p: "+p);
	    			if(p.charAt(0)=='R')
	    				expandedRuleString = rules.get(Integer.valueOf(p.substring(1))).getExpandedRuleString()+sb.toString();
	    			else
	    				expandedRuleString = p+" "+sb.toString();
	    			numEndPos = mapToOriginalTS.get(c+length);
	    			r0[c-1] = "R"+rules.size();
	    		
	    			for (int i = c; i<c+length; i++)
	    				r0[i] = null;
	    			r0New.set(r0New.size()-1,"R"+rules.size());
	    			
	    		    posR0 = c-1;
	    		    
	    		}
	//	    	//1System.out.println("Rule String: " + ruleString);
		//    	//1System.out.println("expe String: "+expandedRuleString);
	    		
	    		c = c + length;
	    		GrammarRuleRecord newRule = new GrammarRuleRecord(rules.size(),ruleString, expandedRuleString, posR0);//, numStartPos, numEndPos);
	    		rules.addRule(newRule,rules.size());
	    		
	    	}
	    	else
	    		{
	    			r0New.add(s);
	    			c++;
	    		}
	}
	
	r0Ori = r0;
	//r0 = new String[r0New.size()];
	
	
	
//	r0 = new String[r0New.size()];
	//1System.out.print("r0Ori = [");
	for (int i = 0; i<r0Ori.length; i++)
		{
	
		//1System.out.print(r0Ori[i]+", ");
		}
	
	//1System.out.println("]");
	
	//1System.out.print("r0New = [");
	for (int i = 0; i<r0New.size(); i++)
		{
	//	r0[i] = r0New.get(i);
		//1System.out.print(r0New.get(i)+", ");
		}
	
	//1System.out.println("]");
	for (int i = 0; i<r0.length; i++)
	{
	//1System.out.print(r0[i]+", ");
	}

//1System.out.println("]");
	}

	/*
	private void AnomalyDetection() {
		ArrayList<RuleInterval> anomalyCandidate = new ArrayList<RuleInterval>();
		int start = 0;
		int end = 0;
		for(int i = 0; i<r0.length; i++){
			
			if(!isNumeric(r0[i])){
				end = getPositionInOriginalTrimedString(i)-1;
				if(end>0) 
					{
						RuleInterval ruleInterval = new RuleInterval(start,end);
						anomalyCandidate.add(ruleInterval);
						//1System.out.println("r0[i]: "+r0[i]+ruleInterval);//anomalyCandidate.get(anomalyCandidate.size()-1));
					}
				start = end + 2;	
			}
			else if(Integer.valueOf(r0[i])<0){
				end = getPositionInOriginalTrimedString(i)-1;
				if(end>0) 
					{
						RuleInterval ruleInterval = new RuleInterval(start,end);
						anomalyCandidate.add(ruleInterval);
						//1System.out.println("r0[i]: "+r0[i]+ruleInterval);//anomalyCandidate.get(anomalyCandidate.size()-1));
					}
				start = end + 2;
			}
		}
	}
	*/
	/*
	private int getPositionInOriginalTrimedString(int index) {
		int ans = -1;
		int idx = index;
		int iter = allMapToPreviousR0.size()-1;
		//int ans = -1;
		ArrayList<Integer> map = new ArrayList<Integer>();
		while(iter>= 0){
			map = allMapToPreviousR0.get(iter);
			idx = map.get(idx);
		}
		ans = 
		return ans;
	}

*/
private void drawOnMap(){
	  // Generate All Motifs and record them on files respectively.
  trueAnomalyCount = 0;
  falsePositiveCount = 0;
  trueNegativeCount = 0;
  falseNegativeCount = 0;
	for(int i = 0; i<isCovered.length;i++)
		{
			isCovered[i] = true;
			ruleCovered[i] = false;
		}
	finalIntervals = new HashMap<String, ArrayList<RuleInterval>>();
	ruleIntervals = new ArrayList<ArrayList<RuleInterval>>();
	anomalyIntervals = new ArrayList<RuleInterval>();
	routes = new ArrayList<ArrayList<Route>>();
	int anomalyCount = 0;
	int totalRuleCount = 0;
	immergableRuleCount = 0;
  //for (int i = 0 ; i<r0.length; i++){
	int i = 0;
	int cnt = 0;
	int totalRuleLength = 0;
	int amountR0RuleLength = 0;
	int nonTerminalCounter = 0;
	int trajCursor = 1;
  int[] anomalyPerTraj = new int[251];
	int[] pointsPerTraj = new int[251]; 
	
	anomalyPerTraj[0] = 0;
	pointsPerTraj[0] = 0;
	int trueAnomalyTraj = 0;
	int falseAnomalyTraj = 0;
	int falseNegativeTraj = 0;
	int trueNegativeTraj = 0;
	int startTraj = 0;
//	int totalNonTerminal = 0;
  while (i<r0.length){
//   	//???//1System.out.println("i:"+i);
		String s = r0[i];
	//	  //???//1System.out.println(minBlocks+"  = getNextNonTerminal(i) = "+ i +" =  " +getNextNonTerminal(i)+" = "+r0[i]);

  	if(!isNumeric(s)){
  	  nonTerminalCounter++;	
  	  amountR0RuleLength = amountR0RuleLength + countSpaces(RuleDistanceMatrix.parseRule(s));  	
  	 // if(countSpaces(RuleDistanceMatrix.parseRule(s))>=minBlocks){
	    	//  if(countSpaces(RuleDistanceMatrix.parseRule(s))>=2){

  	  if(true){
  	  //  	//???//1System.out.println("r0: "+i+" : "+r0[i]+" : "+RuleDistanceMatrix.parseRule(s));

        int startPos = mapToOriginalTS.get(i);
        int endPos;
     //   //???//1System.out.println("r0.length = "+r0.length);
     //   //???//1System.out.println("r0[i] = "+r0[i]);
        		
        if(isNumeric(r0[i+1])&&Integer.valueOf(r0[i+1])<0)
  	     endPos = mapToOriginalTS.get((i+1))-1;
        else
      	 endPos = mapToOriginalTS.get((i+1));
  	  /*
        int endPos;
        
        if(i+2<mapToOriginalTS.size()&&isNumeric(r0[i+1])&&Integer.valueOf(r0[i+1])>0)
  	   {
      	  endPos = mapToOriginalTS.get((i+2))-1;
      	  i++;
  	   }
        else
      	  endPos = mapToOriginalTS.get((i+1))-1;
      	*/    
        RuleInterval interval = new RuleInterval(startPos,endPos);
  	  for (int a = startPos; a<=endPos; a++){
  		  ruleCovered[a] = true;
  	  }
  	  	if (!finalIntervals.containsKey(s)){
  		finalIntervals.put(s, new ArrayList<RuleInterval>());
  		finalIntervals.get(s).add(interval);
  	  	}
  	    else{
  		finalIntervals.get(s).add(interval);
  	  	}
  	  }
  	  i++;
  	  /*
  	   * 
  	   *   Don't consider subtrajectories < minBlocks as anomalies.
  	   * 
  	   */
  	  /*
  	  else{
  		  int unsatisfiedStartPos = mapToOriginalTS.get(i);
	    	  int unsatisfiedEndPos;
	    	  if(i==(r0.length-1))
	    		  unsatisfiedEndPos = mapToOriginalTS.get(i);
	    	  else
	    	  {
	    		  unsatisfiedEndPos = mapToOriginalTS.get((i+1))-1;
	    	  }
	    	  for(int pos = unsatisfiedStartPos; pos<=unsatisfiedEndPos; pos++)
	    		{
	    		  
	    		  isCovered[pos] = false;
		    	  anomalyCount++;
	
	    		}
  	  }
  	  */
  	}
  	else{
  		int numStartPos = mapToOriginalTS.get(i);
	    	  int numEndPos;
	    	  if(Integer.valueOf(r0[i])>=0){
	    		//  //???//1System.out.println(minBlocks+"  = getNextNonTerminal(i) = "+ i +" =  " +getNextNonTerminal(i)+" = "+r0[i]);
	    			  if ((getNextNonTerminal(i)-i)>=minBlocks){
  	     
	    	//  if((Integer.valueOf(r0[i])>=0)&&(getNextNonTerminal(i)-i)>=alphabetSize/30){
  		  int nextNonTerminal = getNextNonTerminal(i);
  		  
  		  
  		  if(nextNonTerminal>=r0.length)
  			  nextNonTerminal = r0.length-1;
  		  
  		//  //???//1System.out.println("ii:"+i);
  		//  //???//1System.out.println(nextNonTerminal + "MapToOriginalTS.get(nextNonTerminal) = "+mapToOriginalTS.get(nextNonTerminal));
  		  if(isNumeric(r0[nextNonTerminal])) // negative
  			  numEndPos = mapToOriginalTS.get(nextNonTerminal)-1;
  		  else
  			  numEndPos = mapToOriginalTS.get(nextNonTerminal);
  		/*
  		  //???//1System.out.print(cnt+": [");
  		  cnt++;
  		  for (int a = i; a<=nextNonTerminal; a++)
  			  {
  			  	//???//1System.out.print(" "+parseRule(r0[a]));
  			  
  			  }
  		  //???//1System.out.println("]");
  		  */
  		//  //???//1System.out.println("i_nextNon : "+i+":"+nextNonTerminal+"["+numStartPos+"-"+numEndPos);
  		//  numEndPos = mapToOriginalTS.get((i+minBlocks))-1;
  		  RuleInterval ri = new RuleInterval(numStartPos,numEndPos);
	  	   	  anomalyIntervals.add(ri);
	  	   	  anomalyPerTraj[trajCursor] = anomalyPerTraj[trajCursor] + (numEndPos-numStartPos+1);
  	  for(int pos = numStartPos; pos<=numEndPos; pos++)
  		{
  		  
  		  isCovered[pos] = false;
	    	  anomalyCount++;

  		}
  		
  	  i = nextNonTerminal;
	    			}
	    			  else
	    				  i = getNextNonTerminal(i);
	    }
	    else //Negative Number
	    {
	    	pointsPerTraj[trajCursor] = numStartPos-startTraj;
	    	trajCursor++;
	    	i++;
	    	startTraj = numStartPos + 1;
	    }
	    	  
  	}
  }
  //???//1System.out.println("r0.length="+r0.length);
  
  	    coverCount = 0;
	Iterator it = finalIntervals.entrySet().iterator();
	while (it.hasNext()){
		@SuppressWarnings("unchecked")
		Map.Entry<String,ArrayList<RuleInterval>> pair = (Map.Entry<String,ArrayList<RuleInterval>>)it.next();
  	  totalRuleLength = totalRuleLength + countSpaces(RuleDistanceMatrix.parseRule(pair.getKey()));

		ruleIntervals.add(pair.getValue());
	}
	
	totalSubTrajectory = 0;
	for (int i1 = 0; i1<ruleIntervals.size();i1++){
		totalSubTrajectory = totalSubTrajectory + ruleIntervals.get(i1).size();
		ArrayList<RuleInterval> positions = ruleIntervals.get(i1);//chartData.getRulePositionsByRuleNum(filteredRuleMap.get(i));
		int counter = 0;
		ArrayList<Route> route = new ArrayList<Route>();
			
		for (int k=0;k<positions.size();k++)
				  {
					  Route singleRoute = new Route();
					  int startPos = positions.get(k).getStartPos();
						int endPos = positions.get(k).getEndPos();
						/*
						for(int index=startPos; index<=endPos;index++)
							isCovered[index]=true;
							*/
		//				//???//1System.out.println("startPos: "+startPos);
		//				//???//1System.out.println("endPos: " +endPos);
						
					//	//???//1System.out.print("track#: "+counter+":       ");
						for (int j = startPos; j<=endPos; j++){
							
							Location loca = new Location(lat.get(j),lon.get(j));
				
							singleRoute.addLocation(lat.get(j), lon.get(j));
								
						
							
						}
						route.add(singleRoute);
						
						counter++;
				  }
		//		  //???//1System.out.println("position size: "+positions.size());
			//	  //???//1System.out.println("route size: "+route.size());
				
				  //  if(route.size()>2)
				     routes.add(route);

		  }	
	
			
			int startAnomalyPos = 0;
			int endAnomalyPos = 0;
			int anomalyCount1 = 0;
			
			
			
			
			for (int a = 0; a<isCovered.length;a++){
				if(!isCovered[a]){
					anomalyCount1++;
					////???//1System.out.println("i: "+a+"\t block: "+blocks.findBlockIdForPoint(new Location(lat.get(a),lon.get(a))));
				}
				/*
				if(isCovered[a] && a<breakPoint)
					trueNegativeCount++;
				if(isCovered[a] && a>=breakPoint)
					falseNegativeCount++;
				if(!isCovered[a]&& a<breakPoint)
					falsePositiveCount++;
				if(!isCovered[a]&& a>=breakPoint)
					trueAnomalyCount++;
					*/
			}
			int i1 = 0;
			
			
			
			
			while  (i1<isCovered.length){
			//	//???//1System.out.println(i + " isCovered :"+isCovered[i]);
				  if(isCovered[i1])
					  {
					  	
					  	coverCount++;
					  	i1++;
					  	if(i1<isCovered.length && lat.get(i1)>-999 && !isCovered[i1]){
					  		startAnomalyPos = i1;
					  		endAnomalyPos = i1;
					  		
					  		while(i1<isCovered.length && lat.get(i1)>-999&&!isCovered[i1]){
					  			endAnomalyPos = i1;
					  			
					  			i1++;
					  			////???//1System.out.println("inner loop :"+i);
					  		}
					  		
					  	//	RuleInterval ri = new RuleInterval(startAnomalyPos,endAnomalyPos);
				  		//	anomalyIntervals.add(ri);
				  	//	//???//1System.out.println("new intervals :"+anomalyIntervals.size()+" : " + ri);
					  	
					  	}
					  	
					  }
				  /*
				  if(ruleCovered[i1]){
					  
				  }
				  */
				  else
					  i1++;
				
			  }
			int ruleCoverCount = 0;
			
			
			
			for (int a = 0;a<ruleCovered.length;a++){
				if(ruleCovered[a])
					ruleCoverCount++;
			}
			  drawAnomaly();
		/*	  //???//1System.out.println("Cover Count: "+ coverCount);
			  //???//1System.out.println("Anomaly Count/count1: "+ anomalyCount+","+ anomalyCount1 );
			  
			  //???//1System.out.println("isCover rate: " +(double)coverCount/(isCovered.length-trajCounter));
			  //???//1System.out.println("satisfied rules: "+finalIntervals.size()+ " longRuleRate: ");
			  //???//1System.out.println("RuleCoverCount: "+ruleCoverCount+" RuleCoverRate: "+(double)ruleCoverCount/(ruleCovered.length-trajCounter));
			  //???//1System.out.println("total number of rules in R0: "+finalIntervals.size()+ "avg rule length: "+ (double)totalRuleLength/finalIntervals.size());
			  //???//1System.out.println("total number of nonterminals in R0: "+nonTerminalCounter+ " avg rule length in R0: "+ (double)amountR0RuleLength/nonTerminalCounter);
			  //???//1System.out.println("latSize = "+lat.size()+"  normalcount = "+breakPoint+"   anomalyCount = "+(lat.size()-breakPoint));
			*/
			 /* //???//1System.out.println("Confusion Matrix:");
			  //???//1System.out.println("True Anomaly:\t"+ trueAnomalyCount+"\t"+ falseNegativeCount);
			  //???//1System.out.println("False Anomaly:\t"+ falsePositiveCount+"\t"+ trueNegativeCount);
			  */
			  
			  
			  
			  for(int k = 1; k<251; k++){
				  if (((double)anomalyPerTraj[k])/(double)pointsPerTraj[k]>EVAL_ANOMALY_THRESHOLD)             //classified as anomaly
					  {
					  if(k%50==1)
					  
						  trueAnomalyTraj++;
					  else
						  falseAnomalyTraj++;
					  }
				  else   //classified as normal
					  {
					  if(k%50==1)
						  falseNegativeTraj++;
					  else
						  trueNegativeTraj++;
					  }
			  }
			  
			totalTP = totalTP+trueAnomalyTraj;
			totalFN = totalFN+falseNegativeTraj;
			totalFP = totalFP+falseAnomalyTraj;
			totalTN = totalTN+trueNegativeTraj;
			
		
			System.out.println("Confusion Matrix:");
			System.out.println("True Anomaly:\t"+ trueAnomalyTraj+"\t"+ falseNegativeTraj);
			System.out.println("False Anomaly:\t"+ falseAnomalyTraj+"\t"+ trueNegativeTraj);
			  
			  
		//	evaluateResult();
		
}


	private int getNextNonTerminal(int i) {
		int j = i+1;
	//	//1System.out.println("j: "+j);
		while(j<r0.length&&isNumeric(r0[j])&&Integer.valueOf(r0[j])>=0)
		{
			j++;
			
		}
		
		return j;
	}

	private boolean isNumberAhead( int i) {
		
		for (int j = i; j<r0.length&&j<(i+minBlocks);j++)
			{
				if(!isNumeric(r0[j])||Integer.valueOf(r0[j])<0)
			
				return false;
			}
		
		//1System.out.print("r0_"+i+"_"+(i+minBlocks-1)+": [");
		for (int j = i; j<r0.length&&j<(i+minBlocks);j++)
		{
			//1System.out.print(r0[j]+" ");
		}
		//1System.out.println("]");
		
		return true;
	}

	private Integer getPositionsInTS(ArrayList<Integer> mapToPreviousR0,ArrayList<Integer> previousMapToOriginalTS, int index) {
		
	//	if(previousMapToOriginalTS.get(mapToPreviousR0.get(index))==108)
		/*
				//1System.out.println("index = "+index+"    r0"+r0[index]);
				//1System.out.println(" mapToPreviousR0.get(index) ="+mapToPreviousR0.get(index));
				//1System.out.println("  previousMapToOriginalTS.get(mapToPreviousR0.get(index))  ="+previousMapToOriginalTS.get(mapToPreviousR0.get(index)));
	   */
		return previousMapToOriginalTS.get(mapToPreviousR0.get(index));
	}

	/*
	   * evaluation
	   */
	private void evaluateResult() {
		
	
	   
		  double[] evalResult = evaluateMotifs(routes);
		  double avgIntraDistance = evalResult[0];
		  double avgIntraDistanceStdDev = evalResult[1];
		  double minInterDistance = evalResult[2];
		  double avgSilhouetteCoefficient = evalResult[3];
		  
		  
	//	  eBlocks = new Blocks(EVAL_RESOLUTION,latMin, latMax, lonMin, lonMax);    // establish a grid map to evaluate the similarity
		  
	//	  String evalHead = "DataName,PaaSize,AlphabetSize,MinimalContinuousBlocks,NoiseCancellationThreshold\n";
		  
		  try{
		  File evalFile = new File("./evaluation/"+"evaluate_"+(int)(minLink*1000)+"_"+alphabetSize+"_"+minBlocks+"_"+noiseThreshold+"_"+lat.size()+"_"+fileNameOnly);
		  FileWriter fr = new FileWriter(evalFile);
		  String sb1; // = new StringBuffer();
		  //sb1.append(fileNameOnly+",");
		  sb1 = (fileNameOnly+","+minLink+","+alphabetSize+","+minBlocks+","+noiseThreshold+","+runTime+","+avgIntraDistance+","+avgIntraDistanceStdDev+","+ minInterDistance+","+avgSilhouetteCoefficient+","+routes.size()+","+lat.size()+','+totalSubTrajectory+","+coverCount+","+immergableRuleCount+"\n");
		  fr.append(sb1);
		  //1System.out.println(EVALUATION_HEAD);
		  //1System.out.println(sb1);
		 // .append("running time: "+runTime+"\n");
		  //fr.append(sb1);
		//  fr.append(Average distances amon)
		  this.log(EVALUATION_HEAD);
		  this.log(sb1.toString());
		  
		  fr.flush();
		  fr.close();
		  }
		  catch (IOException e){
			 
			  e.printStackTrace();
		  }
		  
		 		
	}

	

	
	
	
	

	  /*
	   * Generate All Motifs and record them on files respectively.
	   */
/*
	private void drawOnMap() {
		 // Generate All Motifs and record them on files respectively.
		 // String header = "type,latitude,longitude";
//		  //1System.out.println("Total rules:"+chartData.getRulesNumber());
		  
		//  ArrayList<SAXMotif> allMotifs = chartData.getAllMotifs();
		//  for (int i=1; i<chartData.getRulesNumber();i++){
		    // create merged rule interval data structure corresponding to "clusters" 
		    //ruleIntervals = new ArrayList<ArrayList<RuleInterval>>();
	//        mapToOriginRules = new ArrayList<HashSet<Integer>>();
	        
		    int totalRuleCount = 0;
		    immergableRuleCount = 0;
		    
		    for(int i=0;i<filter.size();i++){
		    	// getRulePositions() was modified to show the intervals only occurred in R0
		  //  	if(!clusterMap.containsKey(i)&&chartData.getRulePositionsByRuleNum(filter.get(i)).size()>=minBlocks) 
		    		{
		    		    ArrayList<RuleInterval> ri = chartData.getRulePositionsByRuleNum(filter.get(i));
		    		  
		    			{ruleIntervals.add(chartData.getRulePositionsByRuleNum(filter.get(i)));
		    			HashSet<Integer> set = new HashSet<Integer>();
		    			set.add(filter.get(i));
		    			mapToOriginRules.add(set);
		    			totalRuleCount++;
		    			immergableRuleCount++;
		    			}
		    			  if(ri.size()<=2)
		    		    	//1System.out.println("Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+ri);
		    		}
		    }
		    
		    	
		    for(int i = 0; i< clusters.size();i++){
		    	
		    	ArrayList<RuleInterval> mergedIntervals = new ArrayList<RuleInterval>();
		    	HashSet<Integer> set = new HashSet<Integer>();
		    	if(clusters.get(i).size()>0)
		    	{
		    		////1System.out.println("cluster "+i+" : {" +clusters.get(i)+"}");
		    		totalRuleCount = totalRuleCount+clusters.get(i).size();
		    	
		    	for(int r : clusters.get(i)){
		    	
		    		int rule = r; //filter.get(r);
		    		set.add(rule);
		    		if(mergedIntervals.size()==0)
		    			mergedIntervals.addAll(chartData.getRulePositionsByRuleNum(rule));
		    		else{
		    			ArrayList<RuleInterval> newIntervals = chartData.getRulePositionsByRuleNum(rule);
		    			for(int j = 0; j<newIntervals.size();j++){
		    				RuleInterval newComer = newIntervals.get(j);
		    				boolean hasMerged = false;
		    				for(int k = 0; k<mergedIntervals.size();k++){
		    					if(RuleInterval.isMergable(mergedIntervals.get(k),newComer)){
		    						//1System.out.println("I still need merge here.");
		    						//1System.out.println("1:" +mergedIntervals.get(k) +" 2:"+ newComer );
		    							
		    						RuleInterval newInterval = RuleInterval.merge(mergedIntervals.get(k), newComer);
		    						mergedIntervals.set(k, newInterval);
		    						hasMerged = true;
		    						//mergedIntervals.remove(k);
		    						//mergedIntervals.add(newInterval);
		    						//break;
		    						
		    					}
		    					
		    						
		    				}
		    				if(!hasMerged)
		    					mergedIntervals.add(newComer);
		    				
		    			}
		    		}
		    		}
		    	if(mergedIntervals.size()>2)
		    	{	
		    	ruleIntervals.add(mergedIntervals);
		    	mapToOriginRules.add(set);
		    	}
		    	else
		    		//1System.out.println("mergedIntervals.size = "+mergedIntervals.size());
		    	}
		    }
		    //1System.out.println("Immergable Rule  = "+ immergableRuleCount);
		    //1System.out.println("Total Rule Count = "+totalRuleCount);
		    boolean[] isCovered = new boolean[lat.size()];
		    coverCount = 0;
		    for (int i = 0; i<isCovered.length;i++)
		    	isCovered[i] = false;
		    totalSubTrajectory = 0;
			for (int i = 0; i<ruleIntervals.size();i++){
		  	totalSubTrajectory = totalSubTrajectory + ruleIntervals.get(i).size();
		  	
			  {//(countSpaces(chartData.getRule(i).getExpandedRuleString())>minBlocks){
			  ArrayList<RuleInterval> positions = ruleIntervals.get(i);//chartData.getRulePositionsByRuleNum(filteredRuleMap.get(i));
			  
			//  ArrayList<RuleInterval> positions = chartData.getRulePositionsByRuleNum(i);	  
		//	  //1System.out.println("rule" + i+" :  "+ positions);//.get(0).toString());
			  
			  
			  if(true)//(positions.size()>2)
				  //&&chartData.getRule(i).getMeanLength()>1)
			  {
				  
				  
					//  File fname = new File("./rules/motif_"+i+".csv");
					//  FileWriter motifPos = new FileWriter(fname);
				//Generating evaluation file
					  
					  int counter = 0;
					  ArrayList<Route> route = new ArrayList<Route>();
					 // Integer route0Id =-1;
					 // Integer route1Id =-1;
				  for (int k=0;k<positions.size();k++)
				  {
					  Route singleRoute = new Route();
			//		  motifPos.append(header+"\n");
					  int startPos = positions.get(k).getStartPos();
						int endPos = positions.get(k).getEndPos();
						
						for(int index=startPos; index<=endPos;index++)
							isCovered[index]=true;
		//				//1System.out.println("startPos: "+startPos);
		//				//1System.out.println("endPos: " +endPos);
						boolean firstPoint = true;
						
					//	//1System.out.print("track#: "+counter+":       ");
						for (int j = startPos; j<=endPos; j++){
							
				//			motifPos.append("T,"+lat.get(j)+","+lon.get(j));
						//	if(counter<2){
							Location loca = new Location(lat.get(j),lon.get(j));
						//	  blocks.addPoint2Block(loc);
					//		  Integer idss = new Integer(blocks.findBlockIdForPoint(loca));
					//		  //1System.out.print(idss+", ");
							singleRoute.addLocation(lat.get(j), lon.get(j));
								
						//	}
							if(firstPoint)
							{
					//			motifPos.append(",Track "+counter+",red\n");
								firstPoint = false;
								
							}
					//		else motifPos.append("\n");
							
						}
						route.add(singleRoute);
						
						counter++;
					//	//1System.out.println();
				  }
		//		  //1System.out.println("position size: "+positions.size());
			//	  //1System.out.println("route size: "+route.size());
				
				  //  if(route.size()>2)
				     routes.add(route);
				    
				    
				    
				    
				  //	motifPos.flush();
				  //	motifPos.close();
				  //	//1System.out.println(fname.getName());
			//	  	route.get(0).print();
				  	
				 
			  }
			  
			//  //1System.out.println("motif index: "+motif.getRuleIndex()+"   " +motif.toString());
			  //FileWriter motifPos = new FileWriter(new File("./motif_"+motif.getRuleIndex()+".csv"));
			 
		  	}
		  }	
			for (int i = 0;i<isCovered.length;i++){
				  if(isCovered[i]==true)
					  coverCount++;
			  }
			  //1System.out.println("Cover Count: "+ coverCount);
			  //1System.out.println("cover rate: " +(double)coverCount/isCovered.length);
		 
	}
*/

	public static void printArrayList(ArrayList<Integer> al) {
		if(al == null || al.size()==0)
			System.out.println("Null or empty ArrayList");
		else 
		{	
		//	//1System.out.print("[ ");
			for (int i = 0; i<al.size();i++)
				;
				//1System.out.println(al.get(i)+" ");
			//1System.out.println();
		}
	}

	private boolean isMergable(double[][] distance, ArrayList<HashSet<Integer>> families, int x, int y, HashMap<Integer, Integer> map, double minLink) {
		//boolean mergable = true;
		double dist = 0;
		int counter = 0;
		if(map.containsKey(x)||map.containsKey(y)){
			if(!map.containsKey(x)){
				
				for(int j: families.get(map.get(y)))
				{
					int i = filterMap.get(j);
					dist = dist + distance[x][i]; 
					counter++;
					//if(distance[x][i]>(minLink*2))
						//return false;
			
				}
				dist = dist/counter;
			}
			else if(!map.containsKey(y)){
				for(int j: families.get( map.get(x)))
					{
					int i = filterMap.get(j);
					dist = dist + distance[i][y];
					counter++;
					//if(distance[i][y]>(minLink*2))
					
						//return false;
			
					}
				dist = dist/counter;
			}
			else
			{	
			for (int m : families.get(map.get(x)))
				for(int n : families.get(map.get(y)))

				{
					
				int i =filterMap.get(m);
				int j =filterMap.get(n);
		//		int xSibling = families.get(x).get(i);
		//		int ySibling = families.get(y).get(j);
			//	if(distance[i][j]>(minLink*2))
				//	return false;
				dist = dist+distance[i][j];
				counter++;
				}
			dist = dist/counter;
			}
			if (dist>minLink)
				return false;
		}
		else if(distance[x][y]>(minLink))
			return false;
		
		return true;
	}

	/*
	   * evaluate the distances intr
	   */
	  private double[] evaluateMotifs(ArrayList<ArrayList<Route>> routes) {
		double[] result = new double[4];  
		StringBuffer sb = new StringBuffer();
		eBlocks = new Blocks(EVAL_RESOLUTION,latMin, latMax, lonMin, lonMax);    // establish a grid map to evaluate the similarity
		ArrayList<Double> allDistances = new ArrayList<Double>();     // The ArrayList of Average Distances of each motif
		ArrayList<Double> allStdDev = new ArrayList<Double>();
		ArrayList<Double> allMinimalInterDistances = new ArrayList<Double>();
		ArrayList<ArrayList<ArrayList<Integer>>> allRules = new ArrayList<ArrayList<ArrayList<Integer>>>();
		for (int i=0;i<routes.size();i++){   // iterate all motifs
			ArrayList<ArrayList<Integer>> allTracks = new ArrayList<ArrayList<Integer>>();
			for(int j = 0; j<routes.get(i).size();j++){			// iterate all tracks under each motif
				ArrayList<Integer> trackIds = new ArrayList<Integer>();
				Route tracks = routes.get(i).get(j);
				Integer previousId = (Integer)(-1);
				for(int k = 0;k<tracks.getLats().size();k++){
					Location loc = new Location(tracks.getLats().get(k),tracks.getLons().get(k));
					Integer id = new Integer(eBlocks.findBlockIdForPoint(loc));
					if(!id.equals(previousId)){
						trackIds.add(id);
						previousId = id;
					}
				}
				allTracks.add(trackIds);
			}
			allRules.add(allTracks);   // after the whole loop all trajectory should be represented in seq of Ids
			ArrayList<Double> pairwiseDistances = getSimilarities(allTracks);
			double sums = 0;
			for (int x = 0; x<pairwiseDistances.size(); x++)
				{
					sums = sums+pairwiseDistances.get(x);
		//			//1System.out.print(pairwiseDistances.get(x)+", ");
				}
		//	//1System.out.println();
		//	//1System.out.println("sum of pairwise distance: "+sums);
		//	//1System.out.println("pairSize = "+pairwiseDistances.size());
			double avgDistance = avg(pairwiseDistances);
			allDistances.add(avgDistance);
			
			Double stdDev = (Double)dev(pairwiseDistances);
			
			allStdDev.add(stdDev);
		/*	//1System.out.println("pairwire distances of motif "+i+": mean = "+avgDistance+",  Std.Dev ="+stdDev);
			for (int m = 0; m<pairwiseDistances.size();m++)
				//1System.out.print(" "+pairwiseDistances.get(m));
			//1System.out.println();
			*/
			
		}
		
		// evaluate distances inter-rules
		for(int i = 0; i<allRules.size();i++){
			ArrayList<Double> pairwiseInterDistances = new ArrayList<Double>();
			for(int j=0; j<allRules.size();j++){
				if(i!=j)
				 pairwiseInterDistances.add(avg(getSimilaritiesInterRules(allRules.get(i),allRules.get(j))));
			
			}
			if(pairwiseInterDistances.size()>0)
				allMinimalInterDistances.add(min(pairwiseInterDistances));

		}
		
//		//1System.out.println("average distances among all motifs: "+avg(allDistances));
	//	//1System.out.println("average standard deviation among all motifs: "+avg(allStdDev));
		//sb.append(avg(allDistances)+","+avg(allStdDev)+"\n");
		
		//return sb.toString();
		result[0] = avg(allDistances);  // avg intra distances
		result[1] = avg(allStdDev);  // avg std. dev. intra distances
		result[2] = avg(allMinimalInterDistances);  // avg minimal inter distances
		ArrayList<Double> silhouetteCoefficients = new ArrayList<Double>();
		for (int i = 0; i<routes.size(); i++)
		{
			
			double sc;
			
			if(allMinimalInterDistances.size()>0&&allDistances.get(i)<allMinimalInterDistances.get(i)){
		//		if(allDistances.get(i)<allMinimalInterDistances.get(i)){

				sc = 1 - allDistances.get(i)/allMinimalInterDistances.get(i);
			}
			else{
				if(allDistances.get(i)==0||allMinimalInterDistances.size()==0)
					sc = 1;
				else
				sc = allMinimalInterDistances.get(i)/allDistances.get(i) - 1;
				
			}
	//		//1System.out.println("compare: "+allDistances.get(i)+"/"+allMinimalInterDistances.get(i)+" = "+sc);
			silhouetteCoefficients.add(sc);
		}
		result[3] = avg(silhouetteCoefficients);
		return result;
	  }

	private double min(ArrayList<Double> list) throws NullPointerException {
		double min;// = -1000000000;
		if(list == null||list.size()==0)
			throw new NullPointerException();
		
			else
			{
				min = list.get(0);
				for(int i = 1; i<list.size();i++)
					if(min>list.get(i))
						min = list.get(i);
				
			}
		return min;
	}

	private ArrayList<Double> getSimilaritiesInterRules(
			ArrayList<ArrayList<Integer>> rule1,
			ArrayList<ArrayList<Integer>> rule2) {
		ArrayList<Double> pairwiseInterDistances = new ArrayList<Double>(); 
		for(int i = 0; i<rule1.size(); i++)
			for (int j=0;j<rule2.size();j++){
				pairwiseInterDistances.add(avgDTWDistance(eBlocks, rule1.get(i),rule2.get(j)));
			}
		
		return pairwiseInterDistances;
	}

	private double dev(ArrayList<Double> pairwiseDistances) {
		double avg = avg(pairwiseDistances);
		double sum = 0;
		for (int i=0; i<pairwiseDistances.size(); i++)
			sum = sum + (pairwiseDistances.get(i)-avg)*(pairwiseDistances.get(i)-avg);
		return Math.sqrt(sum);
	}

	private ArrayList<Double> getSimilarities(ArrayList<ArrayList<Integer>> allTracks) {
		ArrayList<Double> pairwiseDistance = new ArrayList<Double>();
		for(int i = 0; i<allTracks.size();i++){
			for(int j=i+1;j<allTracks.size();j++){
			//	//1System.out.println("i="+i+" j="+j);
				double similarity = avgDTWDistance(eBlocks, allTracks.get(i),allTracks.get(j));
				pairwiseDistance.add(similarity);
			}
		}
		return pairwiseDistance;
	}

	private double avgDTWDistance(Blocks blocks, ArrayList<Integer> s,
			ArrayList<Integer> t) {
		
	//	//1System.out.print("s::::::::::size:"+s.size());
		/*
		for(int i=0; i<s.size();i++)
			//1System.out.print(" "+s.get(i));
			*/
	//	//1System.out.println();
	//	//1System.out.print("t::::::::::size:"+t.size()+"   ");
		/*
		for(int i=0; i<t.size();i++)
			//1System.out.print(" "+t.get(i));
			*/
	//	//1System.out.println();
		int n = s.size();
		int m = t.size();
		double[][] DTW = new double[n+1][m+1];
		double cost = 0;
		for(int i=0;i<n;i++)
			DTW[i+1][0]=Double.MAX_VALUE;
		for(int i=0;i<m;i++)
			DTW[0][i+1] = Double.MAX_VALUE;
		DTW[0][0] = 0;
		for (int i=0;i<n;i++){
			for(int j=0;j<m;j++){
				cost = blocks.distance(s.get(i),t.get(j));
			//	//1System.out.println("cost_"+i+","+j+": "+cost);
				DTW[i+1][j+1]=cost+minimum(DTW[i][j+1],		// insertion
										   DTW[i+1][j], 	// deletion
										   DTW[i][j]);	// match
			}
		}
	//	//1System.out.println("DTW:::::"+DTW[n][m]);
		int step = 1;
		int x = n;
		int y = m;
		while(!((x==1)&&(y==1))){
			step = step + 1;
			switch(min(DTW[x-1][y-1],DTW[x-1][y],DTW[x][y-1])){
			case 1: x--; y--; break;
			case 2: x--; break;
			case 3: y--; break;
			default: //1System.out.println("Error!!!!");
			}
			
		}
	//	//1System.out.println("step: "+step);
		double avg = DTW[n][m]/step;
	//	//1System.out.println("avgDTW:::::"+avg);
		return avg;
	}

	private int min(double d, double e, double f) {
		if(d<=e&&d<=f)
			return 1;
		if(e<=d&&e<=f)
			return 2;
		if(f<=e&&f<=d)
			return 3;
		return 0;
	}

	private double minimum(double a, double b, double c) {
		return Math.min(Math.min(a, b), Math.min(b, c));
	}

	

	private boolean isNoise(Integer id, int i,int noiseThreshold) {
		  if(i<1)
			  return false;
		  if(id.intValue()<0)
	  		{
	  		////1System.out.println("id: "+id);
	  		return false;
	  		}
		  if((i+noiseThreshold)>lat.size())
			  return false;
		  for(int j = 1; j<noiseThreshold;j++)
		  {	Location loc = new Location(lat.get(i+j),lon.get(i+j));
		  	//blocks.addPoint2Block(loc);
		  	Integer currentId = new Integer(blocks.findBlockIdForPoint(loc));
		  	
		  	if(!currentId.equals(id))
		  		{
		  	//	//1System.out.println("id   currentId:  "+id+"       "+currentId);	
		  		return true;
		  		}
		  }
		return false;
	}
	
	 /* 

	private String getTrimedIds(
			ArrayList<NumerosityReductionMapEntry> track) {
		
	//	  String ans;
		  StringBuffer sb = new StringBuffer();
			for (int i = 0; i<track.size();i++){
				sb.append(track.get(i).getValue());
				sb.append(" ");
				
			}
			return sb.toString();
	}
	*/

/*
 * Compute the avg. value of the given ArrayList
 */
	private Double avg(ArrayList<Double> list){
		Double sum= new Double(0);
		for (int i = 0; i< list.size(); i++){
			if(Double.isNaN(list.get(i)))
				throw new NullPointerException();
			sum = sum + list.get(i);
	//		//1System.out.print(list.get(i)+" ");
		}
	//	//1System.out.println();
		
	//	//1System.out.println("sum = "+sum+" avg = "+sum/list.size());
		if(list.size()>0)
		return sum/list.size();
		else
			return -88888.0;
	}
	@SuppressWarnings("rawtypes")
	public ArrayList<Integer> getTrimedPositions(
			ArrayList<NumerosityReductionMapEntry> track) {
		  ArrayList<Integer> ans = new ArrayList<Integer>();
		for (int i = 0; i<track.size();i++){
			ans.add((Integer) track.get(i).getKey());
		}
		return ans;
	}
    public static ArrayList<ArrayList<Route>> getMotifs(){
    	return routes;
    }
    public static ArrayList<Route> getRawTrajectory(){
    	return rawRoutes; 
    }
    public static ArrayList<Route> getAnomaly(){
    	return anomalyRoutes; 
    }
    private void drawAnomaly() {
    	this.anomalyRoutes = new ArrayList<Route>();
		
		
			
  		for (int k=0;k<anomalyIntervals.size();k++)
  				  {
  					  Route singleRoute = new Route();
  					  int startPos = anomalyIntervals.get(k).getStartPos();
  						int endPos = anomalyIntervals.get(k).getEndPos();
  						double distance = 0;
  						/*
  						for(int index=startPos; index<=endPos;index++)
  							isCovered[index]=true;
  							*/
  		//				//1System.out.println("startPos: "+startPos);
  		//				//1System.out.println("endPos: " +endPos);
  						
  					//	//1System.out.print("track#: "+counter+":       ");
  						
  						Location loca = new Location(ncLat.get(startPos),ncLon.get(startPos));
  						//Location endLoc = new Location(lat.get(startPos),lon.get(startPos));
  						
  						for (int j = startPos; j<=endPos; j++){
  							Location previousLoc =loca;
  							loca = new Location(ncLat.get(j),ncLon.get(j));
  							distance = distance + blocks.distance(blocks.findBlockIdForPoint(previousLoc), blocks.findBlockIdForPoint(loca)); 
  							singleRoute.addLocation(ncLat.get(j), ncLon.get(j));
  								
  						
  							
  						}
  						if (distance>0.1) // remove the false anomalies in the same block.
  						{
  						anomalyRoutes.add(singleRoute);
  						}
  				  }
  		//		  //1System.out.println("position size: "+positions.size());
  			//	  //1System.out.println("route size: "+route.size());
  				
  				  //  if(route.size()>2)

  		  	
}

	public static String map2String(HashMap map){
		  String string = new String();
		  Iterator it = map.entrySet().iterator();
		  while(it.hasNext()){
			  
			  string.concat(it.next().toString());
		  }
		  return string;
	  }
		  
	  /**
	   * Performs logging messages distribution.
	   * 
	   * @param message the message to log.
	   */
	  private void log(String message) {
	    this.setChanged();
	    notifyObservers(new SequiturMessage(SequiturMessage.STATUS_MESSAGE, "model: " + message));
	  }

	  /**
	   * Counts spaces in the string.
	   * 
	   * @param str The string.
	   * @return The number of spaces.
	   */
	  public static int countSpaces(String str) {
	    int counter = 0;
	    
	    for (int i = 0; i < str.length(); i++) {
	      if (str.charAt(i) == ' ') {
	        counter++;
	      }
	    }
	//    //1System.out.println("string: "+str+"   length = "+counter);
	    return counter;
	  }
	  public static boolean isNumeric(String str)  
	  {  
		  if(str == null)
			  return false;
	    try  
	    {  
	      double d = Double.parseDouble(str);  
	    }  
	    catch(NumberFormatException nfe)  
	    {  
	      return false;  
	    }  
	    return true;  
	  }

	public static String parseRule(String string) {
		StringBuffer sb = new StringBuffer();
		////1System.out.println("string: "+string);
		ArrayList<String> sa = new ArrayList<String>();
		String[] stringArray = string.split(" ");
		for (String s:stringArray){
			if (s.charAt(0)=='I')
			{
				if(s.contains("r")){
					int rIndex = s.indexOf("r");
					Integer iteration = Integer.valueOf(s.substring(1, rIndex));
					Integer rule = Integer.valueOf(s.substring(rIndex+1));
				//	//1System.out.println("s: "+s+" iteration: "+iteration+" rule: "+rule);
					String subRule = parseRule(allRules.get(iteration).get(rule).getExpandedRuleString());
					sa.add(subRule);
			//		//1System.out.println(s+" = "+subRule );
					
				}
				else if(s.contains("C")){
					int cIndex = s.indexOf("C");
					Integer iteration = Integer.valueOf(s.substring(1, cIndex));
					Integer cluster = Integer.valueOf(s.substring(cIndex+1));
				//	//1System.out.println("s: "+s+" iteration: "+iteration+" cluster: "+cluster);
					Integer ruleInCluster = (Integer)allClusters.get(iteration).get(cluster).toArray()[0];
					String subRule = parseRule(allRules.get(iteration).get(ruleInCluster).getExpandedRuleString());
					sa.add(subRule);
				//	//1System.out.println(s+" = "+subRule );
	
				}
			}
			else if (s.charAt(0)=='R'){
				throw new IllegalArgumentException("expect 'I' encounter 'R'");
			//	int idx = Integer.valueOf( s.substring(1));
				//return rules.get(idx).getExpandedRuleString();
			}
			
			else	//Base Case
			{
				Integer test = Integer.valueOf(s);
				sa.add(s);
		//		//1System.out.println("s: "+ s);
			}
		}
		for (int i = 0; i<sa.size()-1;i++){
			sb.append(sa.get(i));
			sb.append(" ");
		}
		if(sa.size()>0)
		   sb.append(sa.get(sa.size()-1));
		////1System.out.println("sb: "+sb.toString());
		String ans = sb.toString();
		return ans;
	}

}
