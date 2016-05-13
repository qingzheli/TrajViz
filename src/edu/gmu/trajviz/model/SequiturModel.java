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
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Observable;
import java.util.PriorityQueue;
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
import edu.gmu.trajviz.util.Tools;
import net.sf.javaml.core.Instance;
public class SequiturModel extends Observable {
	
	/*
	 * 
	 */
	public static int top1EuDistanceCalled;
	public static int allEuDistanceCalled;
	public static ArrayList<ArrayList<Double>> rawtrajX,rawtrajY,oldtrajX, oldtrajY, trajX, trajY,allTimeline;
	
	public static HashMap<Integer, ArrayList<ArrayList<Center>>> allCenterPoints;
	public static ArrayList<ArrayList<Integer>> rawsaxstrings, saxstrings;
	public static ArrayList<ArrayList<Center>> paats,saxseqs;
//	public static HashMap<Integer, ArrayList<ArrayList<Double>>> allCenterX, allCenterY;
//	public static HashMap<Integer, HashMap<Center, ArrayList<Center>>> allNeighbors; 
	public static HashMap<Integer, ArrayList<Cluster>> allTrajClusters;  //<length, arraylist of clusters>
	public static HashMap<Integer, ArrayList<Cluster>> allMotifs;   // this should be a subset of allTrajClusters with size>1;
	public static HashMap<Integer, ArrayList<String>> allSubseq;
	public static HashMap<Integer, PriorityQueue<RoutePair>> mergablePair;
	public static HashMap<Integer, HashMap<String, Cluster>> clusterMaps;
	public static int longest3Traj[];
	public static double distCut;
	public static ArrayList<ArrayList<Boolean>> isAnomaly;
	public static double R;
	public static double rs;
	public static HashMap<String, HashSet<String>> motifMatches;
	public static HashMap<String, Cluster> subtrajClusterMap;
	public static int slidePoint;
	public static int coverPoint;
	public static HashMap<Integer, ArrayList<String>> anomalyMap;
	public static Map<Integer, ArrayList<String>> sortedAnomalyMap;
	/*
	 * 
	 */
	
	
//	public static double MINLINK = 0.0;
//	public final static double (minLink*2) = 0.0;
	
	public final static int EVAL_RESOLUTION = 100;
	
	final static Charset DEFAULT_CHARSET = StandardCharsets.UTF_8;
	public final static String EVALUATION_HEAD = "DataName,MinLink,AlphabetSize,MinBlocks,NCThreshold,RunningTime,AvgDistance,AvgeStdDev,MinInterDistance,SilhouetteCoefficient,TotalRules,TotalDataPoints, TotalSubTrajectories,CoveredPoints, ImmergableRuleCount\n";
//	public static final int ALPHABETSIZE = 50;
//	public static final int CONTINUALBLOCKTHRESHOLD = 10;
	//public static final int paaSize = 10;
	//public static RuleInterval[][][] allSubTrajectories;
	public double diagnalDistance;
	public long totalPoints;
	private ArrayList<Integer> status;
	private static final String SPACE = " ";
	private static final String CR = "\n";
	private static final int STEP = 2;
	private static final int DEFAULT_TIME_GAP = 6;//180;
//	private static final int DEFAULT_TIME_GAP = 180;
//	private static final int DEFAULT_TIME_GAP = 3;
    private boolean[] isCovered;
    private boolean[] groundTruth;
    private int breakPoint; // the positions<breakPoint are normal, otherwise are abnormal.
    private int trueAnomalyCount;
    private int falsePositiveCount;
    private int trueNegativeCount;
    private int falseNegativeCount;
    private boolean[] ruleCovered;
    
//	private static final int NOISYELIMINATIONTHRESHOLD = 5;
	public static int alphabetSize;
	public static double minLink;
	private int noiseThreshold;
	private static GrammarRules rules;
	public static HashMap<String, ArrayList<String>> allPostions;
	public static ArrayList<GrammarRules> allRules;
	public static TreeMap<String, GrammarRuleRecord> sortedRuleMap;
	public static ArrayList< ArrayList<HashSet<Integer>>> allClusters;
	private ArrayList<HashSet<Integer>> clusters;
//	private Cluster cluster;
//	private HashMap<String, Cluster> currentClusters;
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
	private String discordFileName;
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

	public ArrayList<Double> reTime;
	//public static ArrayList<Double> ncLat = new ArrayList<Double>();
	//public static ArrayList<Double> ncLon = new ArrayList<Double>();

	public static ArrayList<Double> ncLat = new ArrayList<Double>();
	public static ArrayList<Double> ncLon = new ArrayList<Double>();

	//public ArrayList<Double> paaLat;
	//public ArrayList<Double> paaLon;
	public ArrayList<Double> lon;
	public ArrayList<Double> latOri;
	public ArrayList<Double> lonOri;
	public ArrayList<Double> timeLine;
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
	public static int minBlocks; 
	private static Logger consoleLogger;
	private Map<Integer, Double> accumulatDistance;
	public static HashMap<String,Double> allDiscordDistances;
	  private static Level LOGGING_LEVEL = Level.DEBUG;
	public static ArrayList<Double> travelDistance;
	
	  static {
	    consoleLogger = (Logger) LoggerFactory.getLogger(SequiturModel.class);
	    consoleLogger.setLevel(LOGGING_LEVEL);
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

		    consoleLogger.info("setting the file " + filename + " as current data source");

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
			  this.log("unable to load data - no data source select yet");
			  return;
		  }
		  Path path = Paths.get(this.dataFileName);
		  if (!Files.exists(path)){
			  this.log("file"+ this.dataFileName + "doesn't exist.");
			  return;
		  }
	//	  System.out.println("this.dataFileName: "+this.dataFileName);
	//	  System.out.println("this.fileNameOnly : "+this.fileNameOnly);
		  // read the input
		  // init the data array
		  ArrayList<Double> data = new ArrayList<Double>();
		  ArrayList<Double> data1 = new ArrayList<Double>();
		  rawtrajX = new ArrayList<ArrayList<Double>>();
		  rawtrajY = new ArrayList<ArrayList<Double>>();
		  
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
				  rawtrajX.add(new ArrayList<Double>());
				  rawtrajY.add(new ArrayList<Double>());
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
						  rawtrajX.get(rawtrajX.size()-1).add(value);
						  rawtrajY.get(rawtrajY.size()-1).add(value1);
					  
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
						  rawtrajX.add(new ArrayList<Double>());
						  rawtrajY.add(new ArrayList<Double>());
						  rawtrajX.get(rawtrajX.size()-1).add(value);
						  rawtrajY.get(rawtrajY.size()-1).add(value1);
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
			  totalPoints = lineCounter;
		  }
		  catch (Exception e){
			  String stackTrace = StackTrace.toString(e);
			  System.err.println(StackTrace.toString(e));
			  this.log("error while trying to read data from " + this.dataFileName + ":\n" + stackTrace);
		  }
		  
		
			  this.latOri = new ArrayList<Double>();
			  this.lonOri = new ArrayList<Double>();
			  this.lat = new ArrayList<Double>();
			  this.lon = new ArrayList<Double>();
			  
		latMax = Double.valueOf(data.get(0));
		lonMax = Double.valueOf(data1.get(0));
		latMin = Double.valueOf(data.get(0));
		lonMin = Double.valueOf(data1.get(0));
		  for(int i = 0; i<data.size(); i++){
			  double temp_latitude = Double.valueOf(data.get(i));
			  double temp_longitude = Double.valueOf(data1.get(i));
			//  System.out.println("i = "+i+": "+temp_latitude+","+temp_longitude);
			  this.latOri.add(temp_latitude);
			  this.lonOri.add(temp_longitude);
			//  this.lat.add(temp_latitude);
			//  this.lon.add(temp_longitude);
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
		//	  System.out.println(this.lat.get(i)+", "+this.lon.get(i)+","+data.get(i)+", "+data1.get(i));
		  }
		  data = new ArrayList<>();
		  data1 = new ArrayList<>();
		  lat_center = (latMax+latMin)/2;
		  lon_center = (lonMax+lonMin)/2;
		  diagnalDistance = Tools.euDist(latMax, lonMax, latMin, lonMin);
		  System.out.println("lonMax:  "+lonMax+"       lonMin: "+lonMin);
		  System.out.println("latMax:  "+latMax+"       latMin: "+latMin);
		  System.out.println("Number of trajectories: "+trajCounter);
		  consoleLogger.debug("loaded " + this.latOri.size() + " points and "+trajCounter+" Trajecoties... and lineCounter = "+totalPoints);
		  this.log("loaded " + this.latOri.size() + " points from " + this.dataFileName);
		  /*
		  for(int i = 0; i<rawtrajX.size(); i++)
		  {
			  System.out.println(i+"X:"+rawtrajX.get(i));
			  System.out.println(i+"Y:"+rawtrajY.get(i));
		  }
		  */
		  
		  setChanged();
		  notifyObservers(new SequiturMessage(SequiturMessage.TIME_SERIES_MESSAGE, this.latOri,this.lonOri));
		  
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
		  this.allTrajClusters = new HashMap<Integer, ArrayList<Cluster>>();
		  this.allMotifs = new HashMap<Integer, ArrayList<Cluster>>();
		  this.reTime = new ArrayList<Double>();
		  this.isLastIteration = false;
		  this.minLink = minLink;
		  this.minBlocks = minBlocks;
		  this.noiseThreshold = noiseThreshold;
		  this.alphabetSize = alphabetSize;
		  this.allRules = new ArrayList<GrammarRules>();
		  this.allFilters = new ArrayList<ArrayList<Integer>>();
		  this.allClusters = new ArrayList<ArrayList<HashSet<Integer>>>();
		  this.allR0 = new ArrayList<String[]>();
		  this.allMapToPreviousR0 = new ArrayList<ArrayList<Integer>>();
		  this.allMapToOriginalTS = new ArrayList<ArrayList<Integer>>();
		  this.rawRoutes = new ArrayList<Route>();
		  this.anomalyRoutes = new ArrayList<Route>();
		  this.allTrajClusters = new HashMap<Integer, ArrayList<Cluster>>();
		  this.rawAllIntervals = new ArrayList<RuleInterval>();
		  this.subtrajClusterMap = new HashMap<String,Cluster>();
		//  this.currentClusters = new HashMap<String,Cluster>();
		  this.lat = new ArrayList<Double>();
		  this.lon = new ArrayList<Double>();
		  StringBuffer sb = new StringBuffer();
		  if (null == this.latOri ||null == this.lonOri|| this.latOri.size()==0 || this.lonOri.size()==0 ){
			  this.log("unable to \"Process data\" - no data were loaded...");
		  }
		  else{
			  consoleLogger.info("setting up GI with params: ");
			  sb.append(" algorithm: Sequitur");
			  sb.append(" MinLink: ").append(minLink);
			  sb.append(" Alphabet size: ").append(alphabetSize);
			  sb.append(" Minimal Continuous Blocks: ").append(minBlocks);
			  sb.append(" Noise Cancellation Threshold: ").append(noiseThreshold);
			  consoleLogger.info(sb.toString());
			 
			 
			  this.log(sb.toString());
		  }
	//	  buildModel();
		//  drawRawTrajectories();
		  long time = System.currentTimeMillis()/1000;
		  leftPanelRaw();
		  System.out.println("leftPanelRaw time: "+(System.currentTimeMillis()/1000-time));
		   time = System.currentTimeMillis();
		  oldtrajX = rawtrajX;
		  oldtrajY = rawtrajY;
		  resampling();
		  System.out.println("resampling time: "+(System.currentTimeMillis()-time));
	
		 // buildModel();
		  int minLength = minBlocks;
		  int maxLength = minLength;
		//  int maxLength= longest3Traj[2];
		  SequiturModel.R = SequiturModel.distCut*this.minBlocks*SequiturModel.minLink;
		//  R = minLength*distCut*minLink;
		//  rs = R*R;
		//  int maxLength = minLength;
		  System.out.println(longest3Traj[2]);
		  System.out.println("minLength:maxLength = "+minLength+":"+maxLength);
		
		  time = System.currentTimeMillis();
		  
		  findAllMotifSax();
		  System.out.println("findAllMotif time: "+(System.currentTimeMillis()-time));
		   
		 
		//  findMotifs(minLength,maxLength);
		//  findBestMotifs(minLength,maxLength);
		//  findHierarchicalMotifs(minLength, maxLength);
		//  findAllMotifs(minLength,maxLength);
		//  drawMotifs(minLength,maxLength);
		  time = System.currentTimeMillis()/1000;
		  
		 // drawMotifsax();
		  drawAnomalyOnly();
		  System.out.println("drawMotifsax time: "+(System.currentTimeMillis()/1000-time));
		  sortedAnomalyMap = new TreeMap<Integer, ArrayList<String>>(Collections.reverseOrder());
		  Iterator it = anomalyMap.entrySet().iterator();
		  while(it.hasNext()){
			  Entry<Integer,ArrayList<String>> entry = (Entry<Integer, ArrayList<String>>) it.next();
			  sortedAnomalyMap.put(entry.getKey(), entry.getValue());
		  }
		  
		  
		  System.out.println("allMotifs.size = "+allMotifs.size());
		  System.out.println("AnomalyMap = "+anomalyMap);
		  System.out.println("SortedAnomalyMap = "+sortedAnomalyMap);
		  it = sortedAnomalyMap.entrySet().iterator();
		  while(it.hasNext()){
			  Entry<Integer,ArrayList<String>> entry = (Entry<Integer, ArrayList<String>>) it.next();
			  System.out.println(entry.getKey()+","+entry.getValue());
		  }
		//  System.out.println("silcoefMap: "+Evaluation.silCoefMap(allMotifs));
	//	  findTop1Motifs(minLength,maxLength);
		  // evaluation silcoef 2016Apr
		  
		//  notifyObservers(new SequiturMessage(SequiturMessage.CHART_MESSAGE, allMotifs));
		 TrajDiscords.getDiscordsEvaluation();
		  setChanged();
		  notifyObservers(new SequiturMessage(SequiturMessage.CHART_MESSAGE, allMotifs));
		  /*
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
						//	System.out.println("r1: "+r1+" iteration: "+iteration1+" rule1: "+rule1);
							
					//		System.out.println(s+" = "+subRule );
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
						//	System.out.println("r2: "+r2+" iteration2: "+iteration2+" rule2: "+rule2);
							
					//		System.out.println(s+" = "+subRule );
							
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
			  consoleLogger.info("setting up GI with params: ");
			  sb.append(" algorithm: Sequitur");
			  sb.append(" MinLink: ").append(minLink);
			  sb.append(" Alphabet size: ").append(alphabetSize);
			  sb.append(" Minimal Continuous Blocks: ").append(minBlocks);
			  sb.append(" Noise Cancellation Threshold: ").append(noiseThreshold);
			  consoleLogger.info(sb.toString());
			 
			 
			  this.log(sb.toString());
		  }
		  rawAllIntervals = new ArrayList<RuleInterval>();
		  long beginTime = System.currentTimeMillis();
			 // beginTime  =  System.nanoTime();
			  System.out.println("begin time: "+beginTime);
		  buildModel();
		 
		 
		  
		  drawRawTrajectories();
		  
		  allDiscordDistances = new HashMap<String,Double>();
		  TrajDiscords.getAllDiscords();
		  
		  
          allMapToPreviousR0.add(mapToPreviousR0);
          
          
		//  runSequitur();
		
		  
          
            
          //     replace rules with rules' ids and clusters' ids
            
           
		  
		  
          iteration = 0;
          
        /*  System.out.println("before:");
          for (int d = 0; d<words.size(); d++)
        	  System.out.print(words.get(d)+ " ");
          System.out.println();
          */
          
          /*
           * run the algorithm
           */
         /*
          while(hasNewCluster){
        	  
	  		  int lastIteration = iteration;
	  		  hasNewCluster = false;
	  		
        	  System.out.println("Iteration: "+iteration);
        	  
        //  if(hasNewCluster)
        	  
        	  runSequitur(iteration);
        	  iteration = iteration + 1;
        	 // this.minLink = this.minLink*2;
        	  //this.minLink = minLink*(iteration+1);
        	  this.isLastIteration = true;
        	   drawOnMap();
           	System.out.println("total anomalies: "+anomalyRoutes.size());

      }
      
          
          this.isLastIteration = true;
         
       //   drawOnMap();
          
        //  drawOnMap();
         //	System.out.println("total anomalies: "+anomalyRoutes.size());
	  
	  //end while
		  
         
		  System.out.println("Sorted Map.size = "+ sortedRuleMap.size()+ "sortedCounter = "+sortedCounter);
	
	//	  AnomalyDetection();
  
		  
		  
		  
		  this.log("processed data, painting on map");
		  consoleLogger.info("process finished");
		  setChanged();

		
		 


		
		 
		  this.accumulatDistance = sortByValue(this.accumulatDistance);
		  int rank =1;
		  
		  for(Map.Entry<Integer, Double> entry: this.accumulatDistance.entrySet()){
			  int trajId = 0-entry.getKey()-1000;
			  int rankInDiscord = TrajDiscords.allDiscords.get(16).findDiscordByTrajId(trajId);
			  System.out.println("Rank in anomaly: " + rank+" Rank in discords: "+ rankInDiscord + " Trajctory: "+trajId +" Accumulated anomaly distance = "+entry.getValue());
			  rank++;
		  }
		  

		  System.out.println("running time: "+runTime);
		  System.out.println("finalInteravals: "+finalIntervals.size());
		  ArrayList<Integer> frequency = new ArrayList<Integer>();
		 */ 	

		  /*
		  for (int i=0;i<ruleIntervals.size();i++)
			  frequency.add(ruleIntervals.get(i).size());
			  
		  notifyObservers(new SequiturMessage(SequiturMessage.CHART_MESSAGE, this.chartData, ruleIntervals));///, mapToOriginRules));//, frequency ));
		  */
	  //evaluateResult();
	  }
private void findAllMotifSax(){
	paats = this.getCenterArrayList(this.minBlocks);
	allSubseq = new HashMap<Integer, ArrayList<String>>();
	allTrajClusters = new HashMap<Integer,ArrayList<Cluster>>();
	/*
	for(int i=0; i<paats.size();i++){
		System.out.println(i+" paats = "+paats.get(i));
	}
	*/
	paa2saxseqs();
	for(int i = 0; i<blocks.blocks.size(); i++){
		//blocks.blocks.get(i).findHierarchicalMotifs();
		blocks.blocks.get(i).findAnomaly();
	}
	/*
	Iterator it = allTrajClusters.keySet().iterator();
	
	while(it.hasNext()){
		Integer id = (Integer) it.next();
		ArrayList<Cluster> current  = allTrajClusters.get(id);
		*/
	/*
		if(current.size()>0)
		System.out.println(id+":"+current);
		
	}
	*/
	
}
	
	 
				/*
	private void findAllMotifs(int min, int max) {
		allEuDistanceCalled = 0;
	//	allCenterX = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
	//	allCenterY = new HashMap<Integer, ArrayList<ArrayList<Double>>>();
		allCenterPoints = new HashMap<Integer, ArrayList<ArrayList<Center>>>();
	//	allNeighbors = new HashMap<Integer, HashMap<Center, ArrayList<Center>>>(); 
		//establish center ts
		for(int len = min; len<=max; len = len+min){
			slidePoint = (int) (len*minLink);
			R = slidePoint * distCut;
			coverPoint = slidePoint+len;
	//		allNeighbors.put(len, new HashMap<Center, ArrayList<Center>>());
		//	allCenterX.put(len, new ArrayList<ArrayList<Double>>());
		//	allCenterY.put(len, new ArrayList<ArrayList<Double>>());
			ArrayList<ArrayList<Center>> centerArraylist = getCenterArrayList(len);
			allCenterPoints.put(len, centerArraylist);
			
		//	coverPoint = slidePoint+len;
			//	ArrayList<ArrayList<Double>> xcenters = allCenterX.get(len);
			
			
				
		// 
			double minRatio = 100;
			double maxRatio = 0;
			for(int i = 0; i<allCenterPoints.get(len).size(); i++){
				Center c1 = allCenterPoints.get(len).get(i);
				for(int j = i+1; j<allCenterPoints.get(len).size();j++){
					Center c2 = allCenterPoints.get(len).get(j);
					double centerSquareEd = c1.squareDistance(c2);
					double centerEd = Tools.euDist(c1.x, c1.y, c2.x, c2.y);
					Route r1 = Tools.getSubroute(c1.traj, c1.s, c1.l);
					Route r2 = Tools.getSubroute(c2.traj, c2.s, c2.l);
					double avgEd = Tools.routeEuDist(r1, r2);
					double ratio = centerEd/avgEd;
					if(ratio<minRatio)
						minRatio = ratio;
					if(ratio>maxRatio)
						maxRatio = ratio;
					
					System.out.println(i+ "= i :"+c1);
					System.out.println(j+ "= j :"+c2);
					System.out.println("centerED(c1,c2) = "+centerEd);
					
					System.out.println("AvgED(t1,t2) = "+avgEd);
					System.out.println("centerEd/AvgEd = "+ratio);
					
					
				}
				
			}

			System.out.println("minRatio = "+minRatio);
			System.out.println("maxRatio = "+maxRatio);
			*/
			/*
			ArrayList<Center> centers = allCenterPoints.get(len);
			HashMap<Center, ArrayList<Center>> neighbors = new HashMap<Center, ArrayList<Center>>();
	//		double center_R = 1.5*(rs+coverPoint*coverPoint*distCut*distCut);
//			double center_R = distCut*minLink*alphabetSize;
			double center_R = 2*R*R;
			// top
			for(int i = 0; i<centers.size(); i++){
				Center c1 = centers.get(i);
				for(int j = i+1; j<centers.size(); j++){
					Center c2 = centers.get(j);
					if(c1.traj!=c2.traj){
						double squareDistance = c1.squareDistance(c2);
						if(squareDistance<=center_R)
						{
							if(c1.isClose(c2, R/this.noiseThreshold, slidePoint)){
						/*
								System.out.println(c1);
								System.out.println(c1.neighborMap);
							//	c2.isClose(c1, R, len);
								System.out.println(c2);
								System.out.println(c2.neighborMap);
								
							}
						}
							
							
							/*
						{
							if(!neighbors.containsKey(c1)){
								neighbors.put(c1, new ArrayList<Center>());
								neighbors.get(c1).add(c2);
								}
							else{
								boolean isTrivial = false;
								for(int p = 0; p<neighbors.get(c1).size(); p++){
									if(neighbors.get(c1).get(p).traj == c2.traj){
										if(squareDistance<c1.squareDistance(neighbors.get(c1).get(p))){
											neighbors.get(c1).set(p, c2);
										}
										isTrivial = true;
										break;
									}
								}
								if(!isTrivial)
									neighbors.get(c1).add(c2);
							}
							
							if(!neighbors.containsKey(c2)){
								neighbors.put(c2, new ArrayList<Center>());
								neighbors.get(c2).add(c1);	
							}
							else{
								boolean isTrivial = false;
								for(int p = 0; p<neighbors.get(c2).size(); p++){
									if(neighbors.get(c2).get(p).traj == c1.traj){
										if(squareDistance<c2.squareDistance(neighbors.get(c2).get(p))){
											neighbors.get(c2).set(p, c1);
											isTrivial = true;
										}
										break;
									}
								}
								if(!isTrivial)
									neighbors.get(c2).add(c1);
							}
													
						}
						
					}
				}
			}
			
			Iterator it = neighbors.keySet().iterator();
			while(it.hasNext())
			{	
				Center current = (Center) it.next();
				System.out.println(current+" size = "+neighbors.get(current).size()+neighbors.get(current));
			}
			System.out.println("neighbors.size() = "+neighbors.size());
			
			for(int i = 0; i<centers.size(); i++){
				Center c1 = centers.get(i);
				System.out.println(c1+" size = "+c1.neighborMap.size());
				System.out.println(c1.neighborMap);
			}
		}
		
		System.out.println("find all motif R = "+R);
		System.out.println("all Motif EuDistance called = "+this.allEuDistanceCalled);
		
		
	}
	*/

	private ArrayList<ArrayList<Center>> getCenterArrayList(int len) {
		ArrayList<ArrayList<Center>> allCenters = new ArrayList<ArrayList<Center>>();
		for(int traj = 0; traj<oldtrajX.size(); traj++)
		{	
			
			ArrayList<Center> centers = new ArrayList<Center>();
		//	ArrayList<Double> centerX = new ArrayList<Double>();
		//	ArrayList<Double> centerY = new ArrayList<Double>();
			ArrayList<Double> currentX = oldtrajX.get(traj);
			ArrayList<Double> currentY = oldtrajY.get(traj);
			
			
			
			double xAmount = 0;
			double yAmount = 0;
			for(int s = 0; s<len; s++){
				xAmount = xAmount+currentX.get(s);
				yAmount = yAmount+currentY.get(s);
			}
			Center center = new Center(traj,0,xAmount/len,yAmount/len );
			centers.add(center);
			for(int s = 1; s<=currentX.size()-len; s++){
				xAmount = xAmount-currentX.get(s-1)+currentX.get(s+len-1);
				yAmount = yAmount-currentY.get(s-1)+currentY.get(s+len-1);
				center = new Center(traj,s, xAmount/len,yAmount/len );
			//	System.out.println(s + ": xAmount = "+ xAmount/len+", "+yAmount/len);
				centers.add(center);
			}
		
		//	System.out.println(traj+" center: "+centers);
			allCenters.add(centers);
			
		}
		return allCenters;
			
		
	}
	private void drawAnomalyOnly() {
		  anomalyRoutes = new ArrayList<Route>();
		 anomalyMap = new HashMap<Integer, ArrayList<String>>();
		
			  for(int traj = 0; traj<isAnomaly.size();traj++){
				  int pos = 0;
				  int endPos = 0;
				  int startPos = pos;
				  while(pos<isAnomaly.get(traj).size()){
					  if(isAnomaly.get(traj).get(pos))    // is the start position of anomaly
					  {
						  startPos = pos;
						  int prePos =pos;
						  Route singleAnomaly = new Route();
						  double x = oldtrajX.get(traj).get(pos);
						  double y = oldtrajY.get(traj).get(pos);
						  singleAnomaly.addLocation(x, y);
						  pos++;
						  while(pos<isAnomaly.get(traj).size()&&isAnomaly.get(traj).get(pos)){
							 /*
							  if(euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)>distCut*10){
								  System.out.println(rawtrajX.get(traj));
								  System.out.println(oldtrajX.get(traj));
								  System.out.println(rawtrajY.get(traj));
								  System.out.println(oldtrajY.get(traj));
								  
								  throw new IllegalArgumentException("traj = "+traj+"    pos = "+ pos+"    dist = "+euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)+"prevous pos/distCut: "+distCut);
							  }
							  */
							  x =  oldtrajX.get(traj).get(pos);
							  y = oldtrajY.get(traj).get(pos);
							  singleAnomaly.addLocation(x, y);
							  
							  prePos = pos;
							  pos++;
						  }
						  if(singleAnomaly.getLats().size()>=minBlocks)
						  {
							  Integer l = singleAnomaly.getLats().size();
							  anomalyRoutes.add(singleAnomaly);
							  String anomalyString = "T"+traj+"S"+startPos+"L"+l;
							  if(anomalyMap.containsKey(l)){
								  anomalyMap.get(l).add(anomalyString);
								  
							  }
							  else{
								  ArrayList<String> anomalies = new ArrayList<String>();
								  anomalies.add(anomalyString);
								  anomalyMap.put(l, anomalies);
							  }
					//		  System.out.println("Anomalous Trajectory: "+traj+"-"+startPos+"-"+pos);
						  }
					  }
					  else
						  pos++;
				  }
			  }
		  
		  
	}
	private void drawMotifsax() {
		  anomalyRoutes = new ArrayList<Route>();
		 
		  //separate anomalies and motifs
		  for(int id = 0; id<blocks.blocks.size();id++){
			  ArrayList<Cluster> motif = new ArrayList<Cluster>();
			  for(int i = 0; i<allTrajClusters.get(id).size(); i++){
				  if(allTrajClusters.get(id).get(i).getSize()<2)//Anomalies
				  {
					  /*
					  Iterator it = allTrajClusters.get(len).get(i).trajIds.iterator();
					  String s = it.next().toString();
					  int[] subTraj = parseTrajId(s);  //subTraj = {Traj#,Start,End};
					  Route singleAnomaly = new Route();
					  for (int pos = subTraj[1]; pos<=subTraj[2]; pos++){
						  double x = oldtrajX.get(subTraj[0]).get(pos);
						  double y = oldtrajY.get(subTraj[0]).get(pos);
						  singleAnomaly.addLocation(x, y);
					  }
					  System.out.println("Anomalous Trajectory: "+subTraj[0]+"-"+subTraj[1]+"-"+subTraj[2]);
					  anomalyRoutes.add(singleAnomaly);
					  */
				  }
				  //if (allTrajClusters.get(len).get(i).getSize()>=2){
				  else{
				  motif.add(allTrajClusters.get(id).get(i));
				  
				  }
				  
			  }
			  
			  
			  if(motif.size()>0)
				  {
				  	for(int x =0; x<motif.size(); x++){
				  		System.out.println(motif.get(x).trajIds);
				  	}
				  	allMotifs.put(id, motif);
				  }
			  for(int traj = 0; traj<isAnomaly.size();traj++){
				  int pos = 0;
				  int startPos = pos;
				  while(pos<isAnomaly.get(traj).size()){
					  if(isAnomaly.get(traj).get(pos))    // is the start position of anomaly
					  {
						  startPos = pos;
						  int prePos =pos;
						  Route singleAnomaly = new Route();
						  double x = oldtrajX.get(traj).get(pos);
						  double y = oldtrajY.get(traj).get(pos);
						  singleAnomaly.addLocation(x, y);
						  pos++;
						  while(pos<isAnomaly.get(traj).size()&&isAnomaly.get(traj).get(pos)){
							 /*
							  if(euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)>distCut*10){
								  System.out.println(rawtrajX.get(traj));
								  System.out.println(oldtrajX.get(traj));
								  System.out.println(rawtrajY.get(traj));
								  System.out.println(oldtrajY.get(traj));
								  
								  throw new IllegalArgumentException("traj = "+traj+"    pos = "+ pos+"    dist = "+euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)+"prevous pos/distCut: "+distCut);
							  }
							  */
							  x =  oldtrajX.get(traj).get(pos);
							  y = oldtrajY.get(traj).get(pos);
							  singleAnomaly.addLocation(x, y);
							  
							  prePos = pos;
							  pos++;
						  }
						  if(singleAnomaly.getLats().size()>=minBlocks)
						  {
							  anomalyRoutes.add(singleAnomaly);
					//		  System.out.println("Anomalous Trajectory: "+traj+"-"+startPos+"-"+pos);
						  }
					  }
					  else
						  pos++;
				  }
			  }
		  
		  }
	//	routes = new ArrayList<ArrayList<Route>>();
		  
		  
		// print clusters
		  
		  /*
		for(int l = min; l<=max; l++){
			System.out.println(l+"Begin to draw motifs: allMotifs.get(min).size()"+allMotifs.get(l).size());
			for(int i = 0; i<allMotifs.get(l).size(); i++)
			{
				System.out.println("Cluster # "+i+" Size: "+allMotifs.get(l).get(i).getSize());
				Iterator it = allMotifs.get(l).get(i).trajIds.iterator();
				while (it.hasNext()){
					String s = it.next().toString();
					System.out.println(s);
					
				}
			}
		}
		*/
		
	}
	private void drawMotifs(int min, int max) {
		  anomalyRoutes = new ArrayList<Route>();
		 
		  //separate anomalies and motifs
		  for(int len = min; len<=max;len = len+min){
			  ArrayList<Cluster> motif = new ArrayList<Cluster>();
			  for(int i = 0; i<allTrajClusters.get(len).size(); i++){
				  if(allTrajClusters.get(len).get(i).getSize()<2)//Anomalies
				  {
					  /*
					  Iterator it = allTrajClusters.get(len).get(i).trajIds.iterator();
					  String s = it.next().toString();
					  int[] subTraj = parseTrajId(s);  //subTraj = {Traj#,Start,End};
					  Route singleAnomaly = new Route();
					  for (int pos = subTraj[1]; pos<=subTraj[2]; pos++){
						  double x = oldtrajX.get(subTraj[0]).get(pos);
						  double y = oldtrajY.get(subTraj[0]).get(pos);
						  singleAnomaly.addLocation(x, y);
					  }
					  System.out.println("Anomalous Trajectory: "+subTraj[0]+"-"+subTraj[1]+"-"+subTraj[2]);
					  anomalyRoutes.add(singleAnomaly);
					  */
				  }
				  //if (allTrajClusters.get(len).get(i).getSize()>=2){
				  else{
				  motif.add(allTrajClusters.get(len).get(i));
				  
				  }
				  
			  }
			  
			  
			  if(motif.size()>0)
				  {
				  	for(int x =0; x<motif.size(); x++){
				  		System.out.println(motif.get(x).trajIds);
				  	}
				  	allMotifs.put(len, motif);
				  }
			  for(int traj = 0; traj<isAnomaly.size();traj++){
				  int pos = 0;
				  int startPos = pos;
				  while(pos<isAnomaly.get(traj).size()){
					  if(isAnomaly.get(traj).get(pos))    // is the start position of anomaly
					  {
						  startPos = pos;
						  int prePos =pos;
						  Route singleAnomaly = new Route();
						  double x = oldtrajX.get(traj).get(pos);
						  double y = oldtrajY.get(traj).get(pos);
						  singleAnomaly.addLocation(x, y);
						  pos++;
						  while(pos<isAnomaly.get(traj).size()&&isAnomaly.get(traj).get(pos)){
							 /*
							  if(euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)>distCut*10){
								  System.out.println(rawtrajX.get(traj));
								  System.out.println(oldtrajX.get(traj));
								  System.out.println(rawtrajY.get(traj));
								  System.out.println(oldtrajY.get(traj));
								  
								  throw new IllegalArgumentException("traj = "+traj+"    pos = "+ pos+"    dist = "+euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)+"prevous pos/distCut: "+distCut);
							  }
							  */
							  x =  oldtrajX.get(traj).get(pos);
							  y = oldtrajY.get(traj).get(pos);
							  singleAnomaly.addLocation(x, y);
							  
							  prePos = pos;
							  pos++;
						  }
						  if(singleAnomaly.getLats().size()>=minBlocks)
						  {
							  anomalyRoutes.add(singleAnomaly);
							  System.out.println("Anomalous Trajectory: "+traj+"-"+startPos+"-"+pos);
						  }
					  }
					  else
						  pos++;
				  }
			  }
		  
		  }
	//	routes = new ArrayList<ArrayList<Route>>();
		  
		  
		// print clusters
		  
		  /*
		for(int l = min; l<=max; l++){
			System.out.println(l+"Begin to draw motifs: allMotifs.get(min).size()"+allMotifs.get(l).size());
			for(int i = 0; i<allMotifs.get(l).size(); i++)
			{
				System.out.println("Cluster # "+i+" Size: "+allMotifs.get(l).get(i).getSize());
				Iterator it = allMotifs.get(l).get(i).trajIds.iterator();
				while (it.hasNext()){
					String s = it.next().toString();
					System.out.println(s);
					
				}
			}
		}
		*/
		
	}
	
	public static HashMap<Integer, ArrayList<Cluster>> getAllMotifs(){
		return allMotifs;
	}

	

	
	
	private void filterMotifs(int min, int max) {
		  for(int len = min; len<=max;len++){
			  ArrayList<Cluster> motif = new ArrayList<Cluster>();
			  for(int i = 0; i<allTrajClusters.get(len).size(); i++){
				  if(allTrajClusters.get(len).get(i).getSize()<2){
					  
				  }
				  if (allTrajClusters.get(len).get(i).getSize()>3){
					  motif.add(allTrajClusters.get(len).get(i));
				  }
				  
			  }
			  allMotifs.put(len, motif);
		  }
		  
		
	}
	
	private void findHierarchicalMotifs(int minLength,int maxLength) {
		mergablePair = new HashMap<Integer,PriorityQueue<RoutePair>>();
		clusterMaps = new HashMap<Integer, HashMap<String, Cluster>>();
		
		allSubseq = new HashMap<Integer,ArrayList<String>>();
		  for(int len = minLength; len<=maxLength; len = minLength +len){
			  allTrajClusters.put(len, new  ArrayList<Cluster>());
			  allSubseq.put(len, new ArrayList<String>());
			  mergablePair.put(len, new PriorityQueue<RoutePair>());
			  clusterMaps.put(len, new HashMap<String, Cluster>());
		  }
		  for(int i = 0; i<oldtrajX.size(); i++){
			  for (int length = minLength; length<=maxLength; length = length+minLength){
				  for(int s = 0; s<=oldtrajX.get(i).size()-length; s = (int)(s+length*this.minLink) ){
					  String subseqId = "T"+i+"S"+s+"L"+length;
					  allSubseq.get(length).add(subseqId);
				  }
			  }
		  }
		  for(int len = minLength; len<=maxLength; len= len+minLength){
			  
			  for(int i = 0; i<allSubseq.get(len).size(); i++){
				  R = distCut*(len)*minLink;
				//  System.out.println("Threshold = " + R);
				  for(int j = i+1; j<allSubseq.get(len).size();j++){
					//  System.out.println("R" + R);
					  RoutePair pair = new RoutePair(allSubseq.get(len).get(i),allSubseq.get(len).get(j),R);
					  if(pair.dist<=R){
					 // if(!pair.isTrivial&&pair.dist<=R){
					//	  System.out.println(R +" Pair: "+pair);
						  mergablePair.get(len).add(pair);
					  }
				  }
				  
			  }
			 System.out.println("PQ = "+mergablePair.get(len));
			  System.out.println("len = "+ len+ "size = "+ mergablePair.get(len).size()	 + " peak = "+mergablePair.get(len).peek() );
		  }
		  
		  
		  /*
		   *  start hierarchical clustering
		   */
		  for(int len = minLength; len<=maxLength; len= len+minLength){
			  PriorityQueue<RoutePair> pq = mergablePair.get(len);
			  HashMap<String, Cluster> findCluster = clusterMaps.get(len);
			  while(!pq.isEmpty()){
				  RoutePair pair = pq.poll();
				  if(!findCluster.containsKey(pair.r1)&&!findCluster.containsKey(pair.r2)) // if neither belongs to a cluster, create new cluster
				  {
					Cluster cluster = new Cluster(len);
					int[] t1 = Tools.parseTrajId(pair.r1);
					int[] t2 = Tools.parseTrajId(pair.r2);
					cluster.add(t1[0], t1[1]);
					cluster.add(t2[0], t2[1]);
				  findCluster.put(pair.r1, cluster);
				  findCluster.put(pair.r2, cluster);
				  allTrajClusters.get(len).add(cluster);
				  System.out.println("Clusters "+(allTrajClusters.get(len).size()-1)+"\n"+cluster);
				  }
				  else if(!findCluster.containsKey(pair.r1)&&findCluster.containsKey(pair.r2)){    // if r1 is not in cluster r2 is in a cluster
					  int[] t1 = Tools.parseTrajId(pair.r1);
					  Cluster cluster = findCluster.get(pair.r2);
					  cluster.add(t1[0], t1[1]);
					  findCluster.put(pair.r1, cluster);
					  
				  }
				  else if(findCluster.containsKey(pair.r1)&&!findCluster.containsKey(pair.r2)){   // if r2 is not in cluster r1 is in a cluster
					  int[] t2 = Tools.parseTrajId(pair.r2);
					  Cluster cluster = findCluster.get(pair.r1);
					  cluster.add(t2[0], t2[1]);
					  findCluster.put(pair.r2, cluster);
					  
				  }
				  else  
				  {
					 Cluster c1 = findCluster.get(pair.r1);
					 Cluster c2 = findCluster.get(pair.r2);
					 if(c1!=c2){
						// System.out.println("Before Merege C1: "+c1.trajIds);
						// System.out.println("Before Merege C2: "+c2.trajIds); 
						 c1.merge(c2,R,findCluster);
						 if(c1.getRoutes().size()<1){
							 
							 allTrajClusters.get(len).remove(c1);
							 
						 }
						 if(c2.getRoutes().size()<1){
							 allTrajClusters.get(len).remove(c2);
						 }
						 /*
						 if(c1!=null)
							 System.out.println("After  Merege C1: "+c1.trajIds);
						 if(c2!=null)
						 System.out.println("After  Merege C2: "+c2.trajIds);
						 */
					 }
				  }
				  
			  }
			  
			  
		  }
		  
		  
	}
	

	private void findBestMotifs(int minLength,int maxLength) {
		
		  for(int len = minLength; len<=maxLength; len = minLength +len){
			  allTrajClusters.put(len, new  ArrayList<Cluster>());
		  }
		  for (int i = 0; i<oldtrajX.size(); i++){
			  for(int length = minLength; length<=maxLength; length =length+minLength){
				  
				 // for(int s = 0; s<oldtrajX.get(i).size()-length; s = s+this.noiseThreshold){
				//  for(int s = 0; s<oldtrajX.get(i).size()-length; s = s+this.minBlocks/2){
				  for(int s = 0; s<=oldtrajX.get(i).size()-length; s = (int) (s+length*this.minLink)){
				//	  for(int s = 0; s<=oldtrajX.get(i).size()-length; s = (int) (s+length)){
				//  for(int s = 0; s<=oldtrajX.get(i).size()-length; s++){
					      addToBestTrajClusters(i,s,length);
				     
				  }
			  }
		  }
		
		
	}
	private void findTop1Motifs(int minLength,int maxLength) {
		top1EuDistanceCalled = 0;
		int len = minLength;
		int bestMotifCount = 0;
	//	motifMatches = new HashMap<String, HashSet<String>>();
		String bestMotifLocation = null;
		
		
			//  allTrajClusters.put(len, new  ArrayList<Cluster>());
		  for (int i = 0; i<oldtrajX.size(); i++){
	//		  for(int length = minLength; length<=maxLength; length =length+minLength){
				  
				 // for(int s = 0; s<oldtrajX.get(i).size()-length; s = s+this.noiseThreshold){
				//  for(int s = 0; s<oldtrajX.get(i).size()-length; s = s+this.minBlocks/2){
			  	  
			  for(int s = 0; s<=oldtrajX.get(i).size()-len; s++){
				//	  for(int s = 0; s<=oldtrajX.get(i).size()-length; s = (int) (s+length)){
				//  for(int s = 0; s<=oldtrajX.get(i).size()-length; s++){
				  HashMap<String, HashSet<String>> pointer = new HashMap<String, HashSet<String>>();
				  		
					     String current = "T"+i+"S"+s+"L"+len;
					     int count = 0;
					     pointer.put(current, new HashSet<String>());
					     Route currentRoute = Tools.getSubroute(i, s, len);
					     for(int j = 0; j<oldtrajX.size();j++){
					    	 if(i!=j){
					    		 for(int sj = 0; sj<=oldtrajX.get(j).size()-len; sj++){
					    			 String name = "T"+j+"S"+sj+"L"+len;
					    			 Route route = Tools.getSubroute(j, sj, len);
					    			 top1EuDistanceCalled++;
					    			 if(Tools.routeEuDist(currentRoute, route)<=R/this.noiseThreshold){
					    				 {
					    					 count++;
					    					 pointer.get(current).add(name);
					    					 break;
					    				 }
					    				 
					    			 }
					    			 
					    		 }
					    	 }
					     }
					     if(count>bestMotifCount){
					    	 bestMotifCount = count;
					    	 bestMotifLocation = current;
					    	 motifMatches = pointer;
					     }
				     
		//		  }
			  }
			  
		  }
		  System.out.println("top 1 motif Eudistance called = "+this.top1EuDistanceCalled);
		  System.out.println("find top 1 motif R = "+R);
		  System.out.println("bestMotifCount = "+bestMotifCount);
		  System.out.println("bestMotifLocation = "+bestMotifLocation);
		  System.out.println("motifMatches = "+motifMatches);
		  /*
		   Collections.sort(allTrajClusters.get(len));
		 for(int i=0; i<allTrajClusters.get(len).size(); i++){
		//	 System.out.println(allTrajClusters.get(len).get(i));
			 Cluster resultCluster = allTrajClusters.get(len).get(i);
			 int[] subtraj = Tools.parseTrajId(bestMotifLocation);
			 if(resultCluster.findTraj(subtraj)){
				 System.out.println("result Cluster: Rank = "+(allTrajClusters.get(len).size()-i));
				 System.out.println("Size = "+resultCluster.getSize());
				 break;
			 }
		 }
		 */
		
	}
	private void findMotifs(int minLength,int maxLength) {
		
		  for(int len = minLength; len<=maxLength; len = minLength +len){
			  allTrajClusters.put(len, new  ArrayList<Cluster>());
		  }
		  for (int i = 0; i<oldtrajX.size(); i++){
			  for(int length = minLength; length<=maxLength; length =length+minLength){
				  
				 // for(int s = 0; s<oldtrajX.get(i).size()-length; s = s+this.noiseThreshold){
				//  for(int s = 0; s<oldtrajX.get(i).size()-length; s = s+this.minBlocks/2){
				  for(int s = 0; s<=oldtrajX.get(i).size()-length; s = (int) (s+length*this.minLink)){
				//	  for(int s = 0; s<=oldtrajX.get(i).size()-length; s = (int) (s+length)){
				//  for(int s = 0; s<=oldtrajX.get(i).size()-length; s++){

				/*
					  ArrayList<Double> currentSubX = new ArrayList<Double>();
				      ArrayList<Double> currentSubY = new ArrayList<Double>();
				      
						for (int index = 0; index<length; index++){
							currentSubX.add(oldtrajX.get(i).get(s+length));
							currentSubY.add(oldtrajY.get(i).get(s+length));
						}
					  addToAllTrajClusters(currentSubX, currentSubY);
					  */
				      addToAllTrajClusters(i,s,length);
				      /*
				      if(length==48)
				      {
				    	  System.out.println("Length="+length);
				    	  throw new IllegalArgumentException();
				      }
				      */
				  }
			  }
		  }
		  /*
		for( int i = 0; i<oldtrajX.size();i++){
			
				
			for (int j = i+1; j<oldtrajX.size();j++){
				int bestS1 = 0;
				double bestDist = Math.sqrt((oldtrajX.get(i).get(0)-oldtrajX.get(j).get(0))*(oldtrajX.get(i).get(0)-oldtrajX.get(j).get(0))+(oldtrajX.get(j).get(0)-oldtrajX.get(j).get(0))*(oldtrajY.get(j).get(0)-oldtrajY.get(j).get(0)));
				int bestS2 = 0;
				for (int s1 = 0; s1<=oldtrajX.get(i).size()-length; s1++){    // s1 is the start position of existing trajectory
					for (int s2 = 0; s2<oldtrajX.get(j).size()-length; s2++){
						double dist = Math.sqrt((oldtrajX.get(i).get(s1)-oldtrajX.get(j).get(s1))*(oldtrajX.get(i).get(s1)-oldtrajX.get(j).get(s1))+(oldtrajX.get(j).get(s2)-oldtrajX.get(j).get(s2))*(oldtrajY.get(j).get(s2)-oldtrajY.get(j).get(s2)));
	                	if(bestDist>dist){
	                		bestDist = dist;
	                		bestS1 = s1;
	                		bestS2 = s2;
	                		
	                	}
					}
				}
				if(isClose(i,j,bestS1,bestS2,length)){
					clusterTogether(i,j,bestS1,bestS2,length);
				}
				
			}
		}*/
		
	}

	private void addToAllTrajClusters(int trajId, int s, int length) {
		boolean isAdded = false;
		//ArrayList<Cluster> clusterArrayList = allTrajClusters.get(length);
		for(int i = 0; i<allTrajClusters.get(length).size(); i++){
			if(allTrajClusters.get(length).get(i).add(trajId, s))
			{	isAdded = true;
			break;  //once added to one cluster, avoid added to another cluster.
		
			}
		}
		if(!isAdded){
			Cluster cluster = new Cluster(length);
			cluster.add(trajId, s);
			allTrajClusters.get(length).add(cluster);
		}
		
		
		
	}
	private void addToBestTrajClusters(int traj, int s, int length){
		boolean isAdded = false;
		double minDist = Double.MAX_VALUE;
		int minCluster = 0;
		for(int i = 0; i<allTrajClusters.get(length).size(); i++){
			double dist = 0;
			for(int index = 0; index<allTrajClusters.get(length).get(i).repLineX.size(); index++){
				double pairDist = Tools.euDist(allTrajClusters.get(length).get(i).repLineX.get(index),allTrajClusters.get(length).get(i).repLineY.get(index),oldtrajX.get(traj).get(s+index),oldtrajY.get(traj).get(s+index));
				if(pairDist>distCut*(length)*minLink)
					{
					dist = Double.MAX_VALUE;
					break;
					}
				else
				dist = dist+pairDist;
			}
			if (dist<minDist)
				{
					minDist = dist;
					minCluster = i;
				}
			if(length!=allTrajClusters.get(length).get(i).repLineX.size())
				throw new IllegalArgumentException(length+"         repLineX.size = "+allTrajClusters.get(length).get(0).repLineX.size());
			
			
		}
		
		if((minDist/length)<=(distCut*(length)*minLink)){
			//isAdded = true;
		System.out.println("trajId = "+traj);
		System.out.println("minDist/length = "+minDist/length);
		System.out.println("R      = "+distCut*(length)*minLink);
		System.out.println("minCluster = "+ minCluster);
			allTrajClusters.get(length).get(minCluster).addBest(traj, s);
		}
		else{
			Cluster cluster = new Cluster(length);
			cluster.addBest(traj, s);
			allTrajClusters.get(length).add(cluster);
		}
		
	}

	
	private boolean isClose(int i, int j, int bestS1, int bestS2, int length) {
		
		for(int index = 0; index<length; index++){
			if (Tools.euDist(trajX.get(i).get(bestS1+index),trajY.get(i).get(bestS1+index),trajX.get(j).get(bestS2+index),trajY.get(j).get(bestS2+index))>this.minLink)
				return false;
			
		}	
		return true;
	}
	
private void paa2saxseqs() {
	      saxseqs =  new ArrayList<ArrayList<Center>>();
	//	  rawsaxstrings = new ArrayList<ArrayList<Integer>>();
	//	  saxstrings = new ArrayList<ArrayList<Integer>>();
		  blocks = new Blocks(alphabetSize,latMin,latMax,lonMin,lonMax);
		  double latCut = blocks.latCut;
		  double lonCut = blocks.lonCut;
		
		  words = new ArrayList<String>();
		  
	
		  
		  for (int traj = 0; traj<paats.size();traj++){
			  
			  ArrayList<Integer> rawsaxstring = new ArrayList<Integer>();
			  ArrayList<Center> paa = paats.get(traj);
			  ArrayList<Center> saxseq = new ArrayList<Center>();
			  Integer previousid = -1;
			//  int previousPos = 0;
			  for(int i = 0; i<paa.size();i++){
			   
			  Location loc = new Location(paa.get(i).x,paa.get(i).y);
			//  blocks.addPoint2Block(loc); this should not work here because the point will change if it is a noisy point.
			  Integer id = new Integer(blocks.findBlockIdForPoint(loc));
			  Center center = paats.get(traj).get(i);
			  center.blockid = id;
			  rawsaxstring.add(id);
			  if(!id.equals(previousid)){
				  
				  saxseq.add(center);
				  blocks.addCenter2Block(center);
				  if(saxseq.size()>1){
					  int start = saxseq.get(saxseq.size()-1).s;
					  int end = paats.get(traj).get(i).s;
					  
					  saxseq.get(saxseq.size()-2).e = end;
					/*  
					  for(int pos = saxseq.get(saxseq.size()-1).s; pos<=end; pos++){
						  
					  }
					  blocks.blocks.get(previousid).addSubtrajectory(traj,);
					  */
				  }
				  previousid = id;  

			  }
			 
		
			  }
			  saxseq.get(saxseq.size()-1).e = paa.size()-1;
			  saxseqs.add(saxseq);
		  }
		  //test
		  /*
		  for (int traj = 0; traj<paats.size();traj++){
			
			  System.out.println(traj+" paats = "+paats.get(traj));
			  System.out.println(traj+" saxse= "+saxseqs.get(traj));
		  }	
		  
		  for(int i = 0; i<blocks.blocks.size(); i++){
			  System.out.println(i+" Block Center Count = "+blocks.blocks.get(i).centers.size());
		  }
*/
		/*
			 
			words.add(id.toString());
			
		//	  System.out.println("previousId, id:  "+previousId+",   "+id+"        i:   "+i+"   Lat,Lon: "+paaLat.get(i)+","+paaLon.get(i));
			  Integer trimedIndex = 0;
			  if (!id.equals(previousId))
			  {
				  
				  NumerosityReductionMapEntry<Integer, String> entry = new NumerosityReductionMapEntry<Integer, String>(new Integer(i),id.toString());
		//		  System.out.println("entry: "+i+","+id);
				  trimedTrack.add(entry);
				  mapTrimed2Original.add(i);
				  mapToOriginalTS.add(i);
				  //put the new <index,id> pair into a map 
				  //NumerosityReductionMapEntry entry = new NumerosityReductionMapEntry(i,id);
			//	  trackMap.put(i, id);
				  previousId = id;
			  }
			  
			  				  
		  }
		  */
		  System.out.print("mapTrimed2Original: ");
	//	  printArrayList(mapTrimed2Original);		  
		  /*
		   * Following is put the cleaned location data into block again
		   */
		  /*
		  for(int i=0; i<20;i++)
		  System.out.println("Orignal String: " + words.get(i));
		  */
		 // System.out.println("StringTrimedTrack:  "+trimedTrack);
	/*
		  for(int i=0; i<trimedTrack.size();i++){
			  System.out.println(i+" : "+trimedTrack.get(i).getValue()+" ");
		  }
		  */
		  System.out.println();
		 
		 		
	}
	private void buildModel() {
		
		routes = new ArrayList< ArrayList<Route>>();
	
		//  paaLat = new ArrayList<Double>();
		//  paaLon = new ArrayList<Double>();
		//  ArrayList<Double> latBuffer=new ArrayList<Double>();
		//  ArrayList<Double> lonBuffer=new ArrayList<Double>();
		  double avgLat;
		  double avgLon;

		   // System.out.println("Should not see this msg.");
		  
		  
		  
		  
		  

	//	  System.out.println("oriLon" + lon);
		 
		  
		  blocks = new Blocks(alphabetSize,latMin,latMax,lonMin,lonMax);
		  double latCut = blocks.latCut;
		  double lonCut = blocks.lonCut;
		 // ncLat = new ArrayList<Double>();
		 // ncLon = new ArrayList<Double>();
		 // resample();
		//  lat = latOri;
		 // lon = lonOri;
		  
		  isCovered= new boolean[lat.size()];
		  ruleCovered = new boolean[lat.size()];
		  for(int i=0;i<lat.size();i++){
			  isCovered[i] = true;
			 // System.out.println(lat.get(i)+" , "+lon.get(i));
			
			  
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
		  }
		/*
			 
			words.add(id.toString());
			
		//	  System.out.println("previousId, id:  "+previousId+",   "+id+"        i:   "+i+"   Lat,Lon: "+paaLat.get(i)+","+paaLon.get(i));
			  Integer trimedIndex = 0;
			  if (!id.equals(previousId))
			  {
				  
				  NumerosityReductionMapEntry<Integer, String> entry = new NumerosityReductionMapEntry<Integer, String>(new Integer(i),id.toString());
		//		  System.out.println("entry: "+i+","+id);
				  trimedTrack.add(entry);
				  mapTrimed2Original.add(i);
				  mapToOriginalTS.add(i);
				  //put the new <index,id> pair into a map 
				  //NumerosityReductionMapEntry entry = new NumerosityReductionMapEntry(i,id);
			//	  trackMap.put(i, id);
				  previousId = id;
			  }
			  
			  				  
		  }
		  */
		  System.out.print("mapTrimed2Original: ");
	//	  printArrayList(mapTrimed2Original);		  
		  /*
		   * Following is put the cleaned location data into block again
		   */
		  /*
		  for(int i=0; i<20;i++)
		  System.out.println("Orignal String: " + words.get(i));
		  */
		 // System.out.println("StringTrimedTrack:  "+trimedTrack);
	/*
		  for(int i=0; i<trimedTrack.size();i++){
			  System.out.println(i+" : "+trimedTrack.get(i).getValue()+" ");
		  }
		  */
		  System.out.println();
		 
		 		
	}
      private void resampling(){
    	  oldtrajX = new ArrayList<ArrayList<Double>>();
		  oldtrajY = new ArrayList<ArrayList<Double>>();
		  travelDistance = new ArrayList<Double>();
		  isAnomaly = new ArrayList<ArrayList<Boolean>>();
		  allTimeline = new ArrayList<ArrayList<Double>>();
		  ArrayList<ArrayList<Double>> rawTimeResample = new ArrayList<ArrayList<Double>>(); //  present reTime in ArrayList
		  double amountDist = 0;  
		  longest3Traj = new int[3];
		  longest3Traj[0] = 0;
		  longest3Traj[1] = 0;
		  longest3Traj[2] = 0;
		  double longestDist = 0;
			for (int i = 0; i<rawtrajX.size(); i++){
				rawTimeResample.add(new ArrayList<Double>());
				rawTimeResample.get(i).add(0.0);
				double trajDist = 0;
				for(int j = 1; j<rawtrajX.get(i).size(); j++){
					double dist = Math.sqrt(squareDist(rawtrajX.get(i).get(j-1),rawtrajY.get(i).get(j-1),rawtrajX.get(i).get(j),rawtrajY.get(i).get(j)));
					double eudist = Tools.euDist(rawtrajX.get(i).get(j-1),rawtrajY.get(i).get(j-1),rawtrajX.get(i).get(j),rawtrajY.get(i).get(j));
					if(eudist!=dist){
						throw new IllegalArgumentException(dist+" : "+
					eudist);
					}
					
					
					trajDist = trajDist+dist;
					rawTimeResample.get(i).add(trajDist);
				}
				if(trajDist>longestDist)
					longestDist = trajDist;
				travelDistance.add(trajDist);
				amountDist = amountDist+trajDist;
			//	System.out.println("amountDist = "+amountDist);
			//	System.out.println(i+"  :rawTimeResample.get(i) = "+rawTimeResample.get(i));
			}
				
	//		double distCut=amountDist/(alphabetSize*1000);
		//	double distCut=amountDist/(alphabetSize*rawtrajX.size());
		//	double distCut=amountDist/(totalPoints*alphabetSize);
			distCut = diagnalDistance/this.noiseThreshold;
			int maxPoints = (int) (longestDist/distCut);
			System.out.println("MaxPoint = "+maxPoints);
			if(maxPoints>1000){
				throw new IllegalArgumentException("maxPoints is too large: "+maxPoints);
			}
			if(maxPoints<10){
				throw new IllegalArgumentException("maxPoints is too small: "+maxPoints);
			}
			System.out.println("distCut = "+distCut);
			double minDist = distCut*(minBlocks);
			ArrayList<ArrayList<Double>> timeResample = new ArrayList<ArrayList<Double>>();
			for(int i = 0; i<rawtrajX.size();i++){
				if(travelDistance.get(i)>minDist){
					
				//allTimeline.add(new ArrayList<Double>)
				timeResample.add(new ArrayList<Double>());
				oldtrajX.add(new ArrayList<Double>());
				oldtrajY.add(new ArrayList<Double>());
				isAnomaly.add(new ArrayList<Boolean>());
				int j = oldtrajX.size()-1;
				oldtrajX.get(j).add(rawtrajX.get(i).get(0));
				oldtrajY.get(j).add(rawtrajY.get(i).get(0));
				double time = 0;
				int index = 1;
				timeResample.get(j).add(time);
	//			System.out.println("i===================================================================================="+i);
				while(time<rawTimeResample.get(i).get(rawTimeResample.get(i).size()-1)){
					
					while(time<=rawTimeResample.get(i).get(index)){
						
						if(time+distCut<=rawTimeResample.get(i).get(index)){
						time = time+distCut;
					//	System.out.println(rawTimeResample.get(i).get(index)+"=rawTimeResample.get(i)  |  time = "+time);
						timeResample.get(j).add(time);
						//double ratio = (time-rawTimeResample.get(i).get(index-1))/distCut;
						double ratio = (time-rawTimeResample.get(i).get(index-1))/(rawTimeResample.get(i).get(index)-rawTimeResample.get(i).get(index-1));
						double latitude = rawtrajX.get(i).get(index-1)+(rawtrajX.get(i).get(index)-rawtrajX.get(i).get(index-1))*ratio;//latOri.get(index-1)+(latOri.get(index)-latOri.get(index-1))*ratio;
						double longitude = rawtrajY.get(i).get(index-1)+(rawtrajY.get(i).get(index)-rawtrajY.get(i).get(index-1))*ratio; 
						/*
						
								System.out.println("ratio = "+ratio);
								//System.out.println("resample dist = "+euDist(x,y,rawtrajX.get(i).get(index-1),rawtrajY.get(i).get(index-1)));
								System.out.println("resample dist = "+euDist(x,y,oldtrajX.get(i).get(oldtrajX.get(i).size()-1),oldtrajY.get(i).get(oldtrajY.get(i).size()-1)));
								System.out.println("original dist = "+euDist(rawtrajX.get(i).get(index),rawtrajY.get(i).get(index),rawtrajX.get(i).get(index-1),rawtrajY.get(i).get(index-1)));
								System.out.println("i="+i+"  index = "+index);
								System.out.println(x+","+y);
								System.out.println(rawtrajX.get(i).get(index-1)+","+rawtrajY.get(i).get(index-1));
								System.out.println(rawtrajX.get(i).get(index)+","+rawtrajY.get(i).get(index));
								*/
							if(Tools.euDist(latitude,longitude,oldtrajX.get(j).get(oldtrajX.get(j).size()-1),oldtrajY.get(j).get(oldtrajY.get(j).size()-1))>10*distCut)
						  {	
								throw new IllegalArgumentException();
							  }

						oldtrajX.get(j).add(latitude);
						oldtrajY.get(j).add(longitude);
						}
						else 
							{
							time = time+distCut;
							//	System.out.println(rawTimeResample.get(i).get(index)+"=rawTimeResample.get(i)  |  time = "+time);
								timeResample.get(j).add(time);
								oldtrajX.get(j).add(rawtrajX.get(i).get(index));
								oldtrajY.get(j).add(rawtrajY.get(i).get(index));
							//	break;
							}
					//	timeLine.add(time);
						
						//time = time +distCut;
					}
					index++;
				}
				if(oldtrajX.get(j).size()>longest3Traj[0]){
					longest3Traj[2] = longest3Traj[1];
					longest3Traj[1] = longest3Traj[0];
					longest3Traj[0] = oldtrajX.get(j).size();						
				}
				else if(oldtrajX.get(j).size()>longest3Traj[1]){
					longest3Traj[2] = longest3Traj[1];
					longest3Traj[1] = oldtrajX.get(j).size();
				}
				else if(oldtrajX.get(j).size()>longest3Traj[2]){
					longest3Traj[2] = oldtrajX.get(j).size();
				}
		//		System.out.println(oldtrajX.get(i).size() +"   longest3Traj[2] = "+longest3Traj[2]);
				for(int k = 0; k<oldtrajX.get(j).size();k++){
					isAnomaly.get(j).add(true);
				}
				}	
			}
			
			
			System.out.println("distCut = "+distCut);
			/*
			for(int j = 0; j<oldtrajX.size(); j++)
			  {
				
				  System.out.println(oldtrajX.get(j).size()+"points. "+j+"X:"+oldtrajX.get(j));
				  System.out.println(j+"rawX:"+rawtrajX.get(j));
				  System.out.println(j+"Y:"+oldtrajY.get(j));
				  System.out.println(j+"rawY:"+rawtrajY.get(j));
				  
				  System.out.println(oldtrajX.get(j).size()+"points. "+j+"X:");//+oldtrajX.get(j));
				
			  }		  
			  */
			  
      }
	/*
	 * resample original ts with the same speed.
	 */
	  private  void resample() {
		  oldtrajX = new ArrayList<ArrayList<Double>>();
		  oldtrajY = new ArrayList<ArrayList<Double>>();
		  oldtrajX.add(new ArrayList<Double>());
		  oldtrajY.add(new ArrayList<Double>());
		  timeLine = new ArrayList<Double>();
			double amountDist = 0;  
			reTime.add(0.0);
			int i = 1;
			while(i<latOri.size())
		//	for(int i = 1; i< latOri.size(); i++)
			{	if((latOri.get(i)<=-1000)||(lonOri.get(i)<=-1000))
				{
					reTime.add(latOri.get(i));
					if(i<latOri.size()-1)
						reTime.add(reTime.get(reTime.size()-2)+0.000001);
					
					i=i+2;
				}
				else
				{
					amountDist = amountDist+Math.sqrt(squareDist(latOri.get(i-1),lonOri.get(i-1),latOri.get(i),latOri.get(i)));
				 	reTime.add(amountDist);
				 	i++;
				}
			}
			System.out.println(reTime);
			double distCut=amountDist/(alphabetSize*1000);
			System.out.println("distCut = "+distCut);
			double currentPos = reTime.get(0);
			double time = 0;
			int index = 1;
			lat.add(latOri.get(0));
			lon.add(lonOri.get(0));
			timeLine.add(0.0);
			double reTimeS = 0;
			double reTimeE = 0;
			while(time<=reTime.get(reTime.size()-2)){
				
				if(reTime.get(index)>-1000){
					time = time+distCut;
				
				while(time<=reTime.get(index))
				{
					
					double ratio = (time-reTime.get(index-1))/(reTime.get(index)-reTime.get(index-1));
					double latitude = latOri.get(index-1)+(latOri.get(index)-latOri.get(index-1))*ratio;
					double longitude = lonOri.get(index-1)+(lonOri.get(index)-lonOri.get(index-1))*ratio;
					oldtrajX.get(oldtrajX.size()-1).add(latitude);
					oldtrajY.get(oldtrajY.size()-1).add(longitude);
					lat.add(latitude);
					lon.add(longitude);
					timeLine.add(time);
					
					time = time +distCut;
					
				}
				index++;
				}
				else{
					oldtrajX.add(new ArrayList<Double>());
					oldtrajY.add(new ArrayList<Double>());
					lat.add(reTime.get(index));
					lon.add(reTime.get(index));
					timeLine.add(reTime.get(index));
					index++;
				}
				
			}
			lat.add(reTime.get(reTime.size()-1));
			lon.add(reTime.get(reTime.size()-1));
			timeLine.add(reTime.get(reTime.size()-1));
			ArrayList<Double> actLat = new ArrayList<Double>();
			ArrayList<Double> actLon = new ArrayList<Double>();
			for(int j =0; j<latOri.size(); j++){
			//	if(latOri.get(j)>-1000)
				{
					actLat.add(latOri.get(j));
					actLon.add(lonOri.get(j));
					System.out.println(reTime.get(j)+","+latOri.get(j)+","+lonOri.get(j));
				}
			}
			ArrayList<Double> lat1 = new ArrayList<Double>();
			ArrayList<Double> lon1 = new ArrayList<Double>();
			for(int j =0; j<lat.size(); j++){
			//	if(lat.get(j)>-1000)
				{
					lat1.add(lat.get(j));
					lon1.add(lon.get(j));
					System.out.println(lat1.size()-1+","+timeLine.get(j)+","+lat.get(j)+","+lon.get(j));
					
				}
			}
			/*
			for(int j = 0; j<oldtrajX.size(); j++)
			  {
				  System.out.println(oldtrajX.get(j).size()+"points. "+j+"X:"+oldtrajX.get(j));
				  System.out.println(j+"rawX:"+rawtrajX.get(j));
				  System.out.println(j+"Y:"+oldtrajY.get(j));
				  System.out.println(j+"rawY:"+rawtrajY.get(j));
			  }		
			  */
	}
	  
	  
	  
	/*  
	private void resample(double latCut, double lonCut) {
		double amountDist = 0;  
		reTime.add(0.0);
		for(int i = 1; i< latOri.size(); i++)
		{	if((latOri.get(i)<=-1000)||(lonOri.get(i)<=-1000))
				reTime.add(latOri.get(i));
			else
			{
				amountDist = amountDist+Math.sqrt(squareDist(latOri.get(i-1),lonOri.get(i-1),latOri.get(i),latOri.get(i)));
			 	reTime.add(amountDist);
			}
		}
		double distCut=amountDist/(this.alphabetSize*1000);
		double currentPos = reTime.get(0);
		double time = 0;
		int index = 1;
		this.lat.add(latOri.get(0));
		this.lon.add(lonOri.get(0));
		double reTimeS = 0;
		double reTimeE = 0;
		while(time<=reTime.get(reTime.size()-2)){
			time = time+distCut;
			while(time<=reTime.get(index))
			{
				double ratio = (time-reTimeS)/(reTime.get(index)-reTime.get(index-1));
				double x = latOri.get(index-1)+(latOri.get(index)-latOri.get(index-1))*ratio;
				double y = lonOri.get(index-1)+(lonOri.get(index)-lonOri.get(index-1))*ratio;
				
			}
			index++;
		}
		
		
		
	}
	*/

	public static double squareDist(double lat1, double lon1, double lat2, double lon2) {
		return ((lat1-lat2)*(lat1 -lat2)+(lon1-lon2)*(lon1-lon2));
	}
	
	private void leftPanelRaw(){
		for (int i =0; i<rawtrajX.size();i++){
			Route singleRoute = new Route();
			for(int j = 0; j<rawtrajX.get(i).size(); j++){
			//	Location loca = new Location(oldtrajX.get(i).get(j),oldtrajY.get(i).get(j));
				singleRoute.addLocation(rawtrajX.get(i).get(j),rawtrajY.get(i).get(j));
			}
			rawRoutes.add(singleRoute);
			
		}
	}
	private void drawRawTrajectories() {
		System.out.println("drawRawTraj: "+rawAllIntervals.size());
		
		
	  			
	  		for (int k=0;k<rawAllIntervals.size();k++)
	  				  {
	  					  Route singleRoute = new Route();
	  					  int startPos = rawAllIntervals.get(k).getStartPos();
	  						int endPos = rawAllIntervals.get(k).getEndPos();
	  						/*
	  						for(int index=startPos; index<=endPos;index++)
	  							isCovered[index]=true;
	  							*/
	  		//				System.out.println("startPos: "+startPos);
	  		//				System.out.println("endPos: " +endPos);
	  						
	  					//	System.out.print("track#: "+counter+":       ");
	  						for (int j = startPos; j<=endPos; j++){
	  							
	  							Location loca = new Location(lat.get(j),lon.get(j));
	  				
	  							singleRoute.addLocation(lat.get(j), lon.get(j));
	  								
	  						
	  							
	  						}
	  						rawRoutes.add(singleRoute);
	  						
	  				  }
	  		//		  System.out.println("position size: "+positions.size());
	  			//	  System.out.println("route size: "+route.size());
	  				
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
				  System.out.println("Input String Length: " + countSpaces(saxFrequencyData.getSAXString(SPACE)));
				  consoleLogger.trace("String: " + saxFrequencyData.getSAXString(SPACE));
				//  System.out.println("String: "+ saxFrequencyData.getSAXString(SPACE));
				  consoleLogger.debug("running sequitur...");
				  
				  SAXRule sequiturGrammar = SequiturFactory.runSequitur(saxFrequencyData.getSAXString(SPACE));
				//  System.out.println("sequiturGrammar: "+sequiturGrammar.toGrammarRulesData().getRuleRecord(1));
				  consoleLogger.debug("collecting grammar rules data ...");
				 // GrammarRules rules1 = sequiturGrammar.toGRD();
				 // System.out.println("rules size: "+ rules1.size());			 
		          rules = sequiturGrammar.toGrammarRulesData();
		          rules.setParsedString();
		          realRuleSize = rules.size();
		     //     allRules.add(rules);
		          System.out.println("real rules size: "+ realRuleSize);
		          //debug
		          
		          
		          consoleLogger.debug("mapping rule intervals on timeseries ...");
		          GrammarRuleRecord rule0 = rules.get(0);
		          
		          //String rule0 = rules.get(0).getRuleString();
		          int length3 = countSpaces(rule0.getRuleString());
		        		  
		          r0 = rule0.getRuleString().split(" ");
		          r0Recover = rule0.getRuleString().split(" ");
		          System.out.println("R0 = "+r0);
		          int length4 = r0.length;
		          if (length3!=length4)
	       		  throw new IndexOutOfBoundsException(length3+":"+length4);
		        /* print all rule details
		         */
		          setR0Occ();
		        //  SequiturFactory.updateRuleIntervals(rules, saxFrequencyData, lat.size());   //Both update intervals and intervals in R0
		          for(int i=0;i<rules.size();i++){
		        	  System.out.println("Rule number: "+rules.getRuleRecord(i).getRuleNumber()+" Fre in R0: "+rules.get(i).frequencyInR0()+" LEVEL: "+rules.get(i).getRuleLevel()+" "+rules.get(i)+" StringOccurence: "+rules.getRuleRecord(i).occurrencesToString()+"OccurenceInR0: "+rules.get(i).r0OccurrencesToString()+" Rule String: "+rules.getRuleRecord(i).getExpandedRuleString()+" Rule Positions: "+rules.getRuleRecord(i).getRuleIntervals());
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
		          
		        	  System.out.print("mapToOriginalTS: ");
		        	  ArrayList<Integer> previousMapToOriginalTS = mapToOriginalTS;
		        	 // printArrayList(previousMapToOriginalTS);
		        	  /*new ArrayList<Integer>();
		        	  */
		        	  for (int i = 0; i<mapToOriginalTS.size();i++)
		        	  	{
		        		 // previousMapToOriginalTS.add(mapToOriginalTS.get(i));
		        	  	  System.out.print( mapToOriginalTS.get(i) + " ");
		        	  	}
		             System.out.println();
		             
		              mapToOriginalTS = new ArrayList<Integer>();
		              for(int i = 0; i<r0.length;i++){
		            	          
		              System.out.print(r0[i]+" ");
		              }
		              System.out.println();
					  trimedTrack = new ArrayList<NumerosityReductionMapEntry>();
				  /*
				   * Replace Rules' Ids with Clusters' Ids
				   */
					  System.out.println("r0.length: "+r0.length);
		          for (int i = 0; i<r0.length;i++){
		        	  NumerosityReductionMapEntry<Integer, String> entry;
		        	  if(r0[i]==null)
		        		  {
		        		  	
		        		 // 	Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
        		  		//	mapToOriginalTS.add(pos);
        	  		//		System.out.println("BlockID: " +r0[i]+" : "+pos);//mapTrimed2Original.get(mapToPreviousR0.get(i)));

        		  		//	entry = new NumerosityReductionMapEntry<Integer, String>(pos, null);
        		  			//trimedTrack.add(entry);
		        		//  	continue;
		        		  }
		        	  else{
		        	//  System.out.println("r0_"+i+"="+r0[i] );
		        	  if (r0[i].charAt(0)=='R')
		        		  {
		        		  //	if(i==0)
		        		  	//	System.out.println("r0[i] = "+r0[i]);
		        		  	Integer ruleNumber = Integer.parseInt(r0[i].substring(1));
		        		  	String currentRule = "I"+iteration+"r"+ruleNumber;
		        		  	sortedRuleMap.put(currentRule, rules.get(ruleNumber));
		        		  //	System.out.println("sortedRuleMap.size() = " + sortedRuleMap.size()+" "+currentRule+" : "+rules.get(ruleNumber)+" "+sortedRuleMap);
		        		  	sortedCounter++;
		        		//  	int cursor = rules.get(ruleNumber).getCursor(); 
	
		        		  	if (clusterMap.containsKey(filterMap.get(ruleNumber))){
		        		  			hasNewCluster = true;
		        		  			String s = "I" + (iteration) + "C" + clusterMap.get(filterMap.get(ruleNumber));
		        		  			r0[i] = s;
		        		  			Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
		        		  			mapToOriginalTS.add(pos);
		        	  		//		System.out.println("BlockID: " +r0[i]+" : "+pos);//mapTrimed2Original.get(mapToPreviousR0.get(i)));
	
		        		  			entry = new NumerosityReductionMapEntry<Integer, String>(pos, s);
		        		  			trimedTrack.add(entry);
		        		  		
		        		  	}
		        		  	else{
		        		  			String s = "I" + (iteration) + "r"+ruleNumber;
		        		  			r0[i] = s;
		        		  			Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
		        		  			mapToOriginalTS.add(pos);
	
		        		  		//	System.out.println("RuleID: " +r0[i]+" : "+pos);
	
		        		  			entry = new NumerosityReductionMapEntry<Integer, String>(pos, s);
		        	  				trimedTrack.add(entry);
			
		        		  	
		        		  	}
		        		  		
		        			
		        		  	
		        		  }
		        	  else
		        		  
		        	  {
				  		
		        		   
		        		  	Integer pos = getPositionsInTS(mapToPreviousR0,previousMapToOriginalTS,i);
				  			mapToOriginalTS.add(pos);
	
			  				entry = new NumerosityReductionMapEntry<Integer, String>(pos, r0[i]);
			  			//	System.out.println("BlockID: " +r0[i]+" : "+pos);//mapTrimed2Original.get(mapToPreviousR0.get(i)));
			  				trimedTrack.add(entry);
	
		        	  }
		          }
		          }
		          
		          /*
		          System.out.println("after:");
		          for (int d = 0; d<words.size(); d++)
		        	  System.out.print(words.get(d)+ " ");
		          System.out.println();
		          */
		                
		          
		          System.out.println();
		          allMapToPreviousR0.add(mapToPreviousR0);
		          
		         
		          
		          consoleLogger.debug("done ...");
		          
		          
		          
		          
		          
		          chartData.setGrammarRules(rules);
		          System.out.println("chartData size: "+ chartData.getRulesNumber());
				
		
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
			 // String header = "type,x,y";
	//		  System.out.println("Total rules:"+chartData.getRulesNumber());
			  
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
			    		    	System.out.println("Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+ri);
			    		}
			    }
			    
			    	
			    for(int i = 0; i< clusters.size();i++){
			    	
			    	ArrayList<RuleInterval> mergedIntervals = new ArrayList<RuleInterval>();
			    	HashSet<Integer> set = new HashSet<Integer>();
			    	if(clusters.get(i).size()>0)
			    	{
			    		//System.out.println("cluster "+i+" : {" +clusters.get(i)+"}");
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
			    						System.out.println("I still need merge here.");
			    						System.out.println("1:" +mergedIntervals.get(k) +" 2:"+ newComer );
			    							
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
			    		System.out.println("mergedIntervals.size = "+mergedIntervals.size());
			    	}
			    }
			    System.out.println("Immergable Rule  = "+ immergableRuleCount);
			    System.out.println("Total Rule Count = "+totalRuleCount);
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
			//	  System.out.println("rule" + i+" :  "+ positions);//.get(0).toString());
				  
				  
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
			//				System.out.println("startPos: "+startPos);
			//				System.out.println("endPos: " +endPos);
							boolean firstPoint = true;
							
						//	System.out.print("track#: "+counter+":       ");
							for (int j = startPos; j<=endPos; j++){
								
					//			motifPos.append("T,"+lat.get(j)+","+lon.get(j));
							//	if(counter<2){
								Location loca = new Location(lat.get(j),lon.get(j));
							//	  blocks.addPoint2Block(loc);
						//		  Integer idss = new Integer(blocks.findBlockIdForPoint(loca));
						//		  System.out.print(idss+", ");
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
						//	System.out.println();
					  }
			//		  System.out.println("position size: "+positions.size());
				//	  System.out.println("route size: "+route.size());
					
					  //  if(route.size()>2)
					     routes.add(route);
					    
					    
					    
					    
					  //	motifPos.flush();
					  //	motifPos.close();
					  //	System.out.println(fname.getName());
				//	  	route.get(0).print();
					  	
					 
				  }
				  
				//  System.out.println("motif index: "+motif.getRuleIndex()+"   " +motif.toString());
				  //FileWriter motifPos = new FileWriter(new File("./motif_"+motif.getRuleIndex()+".csv"));
				 
			  	}
			  }	
				for (int i = 0;i<isCovered.length;i++){
					  if(isCovered[i]==true)
						  coverCount++;
				  }
				  System.out.println("Cover Count: "+ coverCount);
				  System.out.println("cover rate: " +(double)coverCount/isCovered.length);
			 
		}
	*/
	/*
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
	*/

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
		        	//	  	System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
		        		  	int length1 = rules.get(currentRule).getRuleYield();
		        		  	int length2 = countSpaces(rules.get(currentRule).getExpandedRuleString());
		        		   // currentIdx = currentIdx + rules.get(currentRule).getRuleYield();
		        		  	currentIdx = currentIdx + length2;
		        	//	  	System.out.println("CurrentIdx = "+currentIdx +" i= "+i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);
		        		  	
		        		  	if(currentIdx>mapToOriginalTS.size()||length1!=length2)
		        		    	
		        		    {
				        		  throw new IndexOutOfBoundsException(i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);

		        		    }
		        		   
		        		    
		        		  }
		        	  else
		        		  {
		        			  
		        		  mapToPreviousR0.add(currentIdx);

		        	//	  System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
		        		  
		        		  	currentIdx++;
		        		  	if(currentIdx>mapToOriginalTS.size())
		        		    	
		        		    {
				        //		  System.out.println(i+" : "+r0[i]+":"+currentIdx);

		        		    }
		        		  }
		        		  
		          }
		          /* above
			          * add on mapToPreviousR0;
			          */
		          
		          
		          
		          System.out.println();
	          		
	}

		private void replaceBack() {
			System.out.println("filter");
			printArrayList(filter);
			System.out.println("filterMap: ");
			System.out.println(filterMap);
			for(int i = realRuleSize; i<rules.size(); i++){
				String[] ruleString = rules.get(i).getRuleString().split(" ");
				System.out.println("Rule "+i+" : "+rules.get(i).getRuleString()+"   filterMap: "+filterMap.get(rules.get(i).getRuleNumber()) );
	        	  if(!clusterMap.containsKey(filterMap.get(rules.get(i).getRuleNumber()))){
	        		
	        		System.out.print("replace Back Rule String: [ ");
	        		int r0Pos = rules.get(i).getR0Occurrences().get(0);
	        		for(int j = 0; j<ruleString.length;j++){
	        			r0[r0Pos+j] = ruleString[j];
	        			System.out.print(ruleString[j]+" ");
	        		}
	        	    System.out.println("]");
	        	  }
	        	
	        		  
	          }
			ArrayList<String> r0new = new ArrayList<String>();
			
			for(int i = 0; i<r0.length; i++)
				if(r0[i]!=null)
					{
						r0new.add(r0[i]);
					}
			r0 = new String[r0new.size()];
			
			System.out.print("r0new: [");
			for(int i = 0; i<r0.length; i++)
				
					{ 
						r0[i] = r0new.get(i);
						System.out.print(" "+r0[i]);
					}
			System.out.println("]");
	}

		private void setR0Occ() {
			HashMap<String, Integer> hm = new HashMap<String, Integer>();
	          
	          
	          
	          
	          
			  
			//  allR0.add(r0);
	          for(int i = 0; i<rules.size();i++){
	        	  String key = rules.get(i).getRuleName();
	        	 // String expandedString = rules.get(i).getExpandedRuleString();
	        	//  System.out.println(rules.get(i));
	        	  hm.put(key, 0);
	          }
	          
	         // System.out.println("R0: "+rule0.getRuleString());
	          System.out.print("r0: ");
	          
	          for(int i = 0; i<r0.length;i++){
	        	          
	          System.out.print(r0[i]+" ");
	          }
	          System.out.println();
	      //    System.out.println(r0);
	          int currentIdx = 0;
	       //   int[] indexes = new int[r0.length];
	          
	          
	          
	          
	          
	         /*
	          * add on mapToPreviousR0;
	          */
	          for(int i=0;i<r0.length;i++){
	        	  
	        	  if(r0[i]==null)
	        		  {
	        	//	  mapToPreviousR0.add(currentIdx);
	        			
	  		        //		  System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
	  		        		  
	  		        		 // 	currentIdx++;
	  		        		  	if(currentIdx>mapToOriginalTS.size())
	  		        		    	
	  		        		    {
	  				        		  System.out.println(i+" : "+r0[i]+":"+currentIdx);

	  		        		    }
	        		  	
	        		  }
	        	 
	        	  else if(r0[i].charAt(0)=='R')
	        		  {
	        		  	Integer currentRule = Integer.valueOf(r0[i].substring(1));
	        		  	hm.put(r0[i], hm.get(r0[i])+1);
	        		  	rules.get(currentRule).addR0Occurrence(currentIdx); // setOccurenceInR0
	        		  	//mapToPreviousR0.add(currentIdx);
	        		  //	System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
	        		  	int length1 = rules.get(currentRule).getRuleYield();
	        		  	int length2 = countSpaces(rules.get(currentRule).getExpandedRuleString());
	        		   // currentIdx = currentIdx + rules.get(currentRule).getRuleYield();
	        		  	currentIdx = currentIdx + length2;
	        		  	System.out.println("CurrentIdx = "+currentIdx +" i= "+i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);
	        		  	if(currentIdx>mapToOriginalTS.size())
	        		    	
	        		    {
			        		  throw new IndexOutOfBoundsException(i+" : "+r0[i]+":"+currentIdx+" expandRule:  "+rules.get(currentRule).getExpandedRuleString()+" length1:length2 = "+length1+":"+length2);

	        		    }
	        		  }
	        	  else
	        		  {
	        			  
	        	//	  mapToPreviousR0.add(currentIdx);

	        //		  System.out.println(i+" : "+r0[i]+":"+currentIdx+" ");
	        		  
	        		  	currentIdx++;
	        		  	if(currentIdx>mapToOriginalTS.size())
	        		    	
	        		    {
			        		  System.out.println(i+" : "+r0[i]+":"+currentIdx);

	        		    }
	        		  }
	        		  
	          }
	          /* above
		          * add on mapToPreviousR0;
		          */
	          
	          
	          
	          System.out.println();
	    //      System.out.print("mapToPreviousR0: ");
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
	        	  System.out.println("Rule number: "+rules.getRuleRecord(i).getRuleNumber()+" Fre in R0: "+rules.get(i).frequencyInR0()+" LEVEL: "+rules.get(i).getRuleLevel()+" "+rules.get(i)+" StringOccurence: "+rules.getRuleRecord(i).occurrencesToString()+"OccurenceInR0: "+rules.get(i).r0OccurrencesToString()+" Rule String: "+rules.getRuleRecord(i).getExpandedRuleString()+" Rule Positions: "+rules.getRuleRecord(i).getR0Intervals());
	          }
	        
	       /*  */
	         
	        filterMap = new HashMap<Integer,Integer>();
	        for (int i = 1; i<rules.size();i++){
	        	System.out.println("Before filter: Frequency in R0: "+ rules.get(i).frequencyInR0()+"  Yield: "+rules.get(i).getRuleYield()+" string: "+rules.get(i).getExpandedRuleString());
					if ((rules.get(i).frequencyInR0()>=1&&countSpaces(RuleDistanceMatrix.parseRule(rules.get(i).getExpandedRuleString()))>=1))//||
						//	(originalRules.get(i).frequencyInR0()>1&&originalRules.get(i).getR0Intervals().size()>2&&originalRules.get(i).getRuleYield()>=minBlocks))
						{
						//HashSet<Integer> set = new HashSet<Integer>();
						System.out.println("Yield: "+rules.get(i).getRuleYield()+" string: "+rules.get(i).getExpandedRuleString());
						filterMap.put(i, filter.size());
						filter.add(i);
					/*	
						if(rules.get(i).getR0Intervals().size()<2)
							System.out.println("Bug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+i);
						*/
						}
					
				}
	        System.out.println("filter Size = "+filter.size());
	        allFilters.add(filter);
	        if(filter.size()>1){
	        //HashMap<Integer,ArrayList<Integer>> mergeRecord = new HashMap<Integer, ArrayList<Integer>>();
	        long t1s = System.currentTimeMillis();
	        RuleDistanceMatrix rdm;
	        System.out.println("AlphabetSize="+this.alphabetSize);
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
	        System.out.println("rdm.pq.size(): "+rdm.pq.size());
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
	      			//	  System.out.println("Adding Line  to a cluster, Line:"+pair.getLine()+" Colu:"+pair.getCol()+clusters.get(clusterMap.get(pair.getCol())));
	      				//  System.out.println("Map:"+clusterMap);
	      				  
	      			  }
	      			  else if(!clusterMap.containsKey(pair.getCol())){
	      				  clusters.get(clusterMap.get(pair.getLine())).add(filter.get(pair.getCol()));
	      				  clusterMap.put(pair.getCol(), clusterMap.get(pair.getLine()));
	      			//	  System.out.println("Adding Colum to a cluster,Colum:"+pair.getCol()+" Colu:"+pair.getCol()+clusters.get(clusterMap.get(pair.getLine())));
	      			//	  System.out.println("Map:"+clusterMap);
	      			  }
	      			  else{
	      				  if(!clusterMap.get(pair.getLine()).equals(clusterMap.get(pair.getCol())))
	      				  {
	      				//  System.out.println("Before Merge, line in cluster:"+clusterMap.get(pair.getLine())+clusters.get(clusterMap.get(pair.getLine()))+" colu in cluster:"+clusterMap.get(pair.getCol())+clusters.get(clusterMap.get(pair.getCol())));
	      				  lineSize = clusters.get(clusterMap.get(pair.getLine())).size();
	      				  colSize = clusters.get(clusterMap.get(pair.getCol())).size();
	      				  clusters.get(clusterMap.get(pair.getLine())).addAll(clusters.get(clusterMap.get(pair.getCol())));
	      				  int colCluster = clusterMap.get(pair.getCol());
	      				  for(int v : clusters.get(clusterMap.get(pair.getCol())))
	      					  {
	      				//	  System.out.print("v: "+v+" ");
	      					  clusterMap.put(filterMap.get(v), clusterMap.get(pair.getLine()));
	      					  clusters.get(clusterMap.get(pair.getLine())).add(v);
	      					  }
	      				  //System.out.println();
	      				  clusters.get(colCluster).clear();
	      				 // System.out.println("After  Merge, Line:"+pair.getLine()+clusters.get(clusterMap.get(pair.getLine()))+" Colu:"+pair.getCol()+clusters.get(colCluster));
	      				 // System.out.println("Map:"+clusterMap);
	      				  totalSize = clusters.get(clusterMap.get(pair.getLine())).size();
	      				  //if((lineSize+colSize)!=totalSize){
	      					//  System.out.println("Error Candidate here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	      				  //}
	      				  }
	      				  //else
	      					// System.out.println("Same Cluster! "+clusterMap.get(pair.getLine())+","+clusterMap.get(pair.getCol()));
	      			  }
	      		  }
	      		  else{
	      			  HashSet<Integer> set = new HashSet<Integer>();
	      			  set.add(filter.get(pair.getLine()));            
	      			  set.add(filter.get(pair.getCol()));
	      			  clusters.add(set);
	      			  clusterMap.put(pair.getLine(), clusters.size()-1);
	      			  clusterMap.put(pair.getCol(), clusters.size()-1);
	      			//  System.out.println("Created a cluster: "+clusters.get(clusters.size()-1));
	      			//  System.out.println("Map:"+clusterMap);
	      		  }
	      		  
	      		  /*
	      		  clusters.get(pair.getLine()).addAll(clusters.get(pair.getCol()));
	      		  clusters.get(pair.getCol()).addAll(clusters.get(pair.getLine()));
	      		  
	      		  for(int i: families.get(pair.getCol()))
	      			  families.get(pair.getLine()).add(i);
	      		  for(int i: families.get(pair.getLine()))
	      			  families.get(pair.getCol()).add(i);
	      			  */	        		
	      	//	  System.out.print("Merged Pair: <"+pair.getLine()+", "+pair.getCol()+"> = "+rdm.matrix[pair.getLine()][pair.getCol()]);
	      	//	  System.out.print(" all distances: ");
	      		  /*
	      		  for (int i : clusters.get(clusterMap.get(pair.getLine())))
	    				for(int j : clusters.get(clusterMap.pair.getCol()))
	    				{
	    				
	    				System.out.print(formatter.format(rdm.matrix[i][j])+", ");
	    					
	    				}*/
	      	//	  System.out.println();
	      	  }
	        }
	        
		  
	        
	        
	        System.out.println("MergableCount: "+mergableCount);
	        
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
	    	System.out.println("build matrix: "+(double)(buildMatrixTime/1000.0));
			  System.out.println("Clustering Time: "+(clusterTime/1000.0));
	        /*
	        for(int i = 0; i<clusters.size();i++){
	      	  System.out.println("i = "+i+" : "+clusters.get(i));
	        }
	        */
	       
			  System.out.println("cluster map size = "+ clusterMap.size());
			    System.out.println("clusterMap:   "+clusterMap);
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
    	//	  System.out.println(minBlocks+"  = getNextNonTerminal(i) = "+ i +" =  " +getNextNonTerminal(i)+" = "+r0[i]);
			
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
	    			System.out.println("p: "+p);
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
	//	    	System.out.println("Rule String: " + ruleString);
		//    	System.out.println("expe String: "+expandedRuleString);
	    		
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
	System.out.print("r0Ori = [");
	for (int i = 0; i<r0Ori.length; i++)
		{
	
		System.out.print(r0Ori[i]+", ");
		}
	
	System.out.println("]");
	
	System.out.print("r0New = [");
	for (int i = 0; i<r0New.size(); i++)
		{
	//	r0[i] = r0New.get(i);
		System.out.print(r0New.get(i)+", ");
		}
	
	System.out.println("]");
	for (int i = 0; i<r0.length; i++)
	{
	System.out.print(r0[i]+", ");
	}

System.out.println("]");
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
						System.out.println("r0[i]: "+r0[i]+ruleInterval);//anomalyCandidate.get(anomalyCandidate.size()-1));
					}
				start = end + 2;	
			}
			else if(Integer.valueOf(r0[i])<0){
				end = getPositionInOriginalTrimedString(i)-1;
				if(end>0) 
					{
						RuleInterval ruleInterval = new RuleInterval(start,end);
						anomalyCandidate.add(ruleInterval);
						System.out.println("r0[i]: "+r0[i]+ruleInterval);//anomalyCandidate.get(anomalyCandidate.size()-1));
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

		Comparator<Double> doubleComparator = new Comparator<Double>() {
	        @Override public int compare(Double s1, Double s2) {
	            return s1.compareTo(s2);
	        }           
	    };
		    this.accumulatDistance = new TreeMap<Integer,Double>();
		    double currentAmountDistance = 0.0;
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
		  //	int totalNonTerminal = 0;
		  	int nullCounter =0;
		    while (i<r0.length){
		 //   	System.out.println("i:"+i);
		    	
		    	if(r0[i]==null)
		    		{
		    			i++;
		    			nullCounter++;
		    			continue;
		    		}
		    		
		  		String s = r0[i];
	    	//	  System.out.println(minBlocks+"  = getNextNonTerminal(i) = "+ i +" =  " +getNextNonTerminal(i)+" = "+r0[i]);

		    	if(!isNumeric(s)){
		    	  nonTerminalCounter++;	
		    	  amountR0RuleLength = amountR0RuleLength + countSpaces(RuleDistanceMatrix.parseRule(s));  	
		    	 // if(countSpaces(RuleDistanceMatrix.parseRule(s))>=minBlocks){
			    	//  if(countSpaces(RuleDistanceMatrix.parseRule(s))>=2){

		    	  if(true){
		    	  //  	System.out.println("r0: "+i+" : "+r0[i]+" : "+RuleDistanceMatrix.parseRule(s));
		
		          int startPos = mapToOriginalTS.get(i-nullCounter);
		          int endPos;
		          if(isNumeric(r0[i+1])&&Integer.valueOf(r0[i+1])<0)
		    	     endPos = mapToOriginalTS.get((i-nullCounter+1))-1;
		          else
		        	 endPos = mapToOriginalTS.get((i-nullCounter+1));
		          
		          
		        //  currentAmountDistance = currentAmountDistance+allDiscordDistances.get(startPos+","+endPos);
		          
		          
		          
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
>>>>>>> refs/heads/discordsEvaluation
			
<<<<<<< HEAD
=======
			    		}
		    	  }
		    	  */
		    	}
		    	else{
		    		int numStartPos = mapToOriginalTS.get(i-nullCounter);
			    	  int numEndPos;
			    	  if(Integer.valueOf(r0[i])>=0){
			    		//  System.out.println(minBlocks+"  = getNextNonTerminal(i) = "+ i +" =  " +getNextNonTerminal(i)+" = "+r0[i]);
			    	  if ((getNextNonTerminal(i)-i)>=minBlocks){
		    	   
			    	//  if((Integer.valueOf(r0[i])>=0)&&(getNextNonTerminal(i)-i)>=alphabetSize/30){
		    		  int nextNonTerminal = getNextNonTerminal(i);
		    		  
		    		  
		    		  if(nextNonTerminal>=r0.length)
		    			  nextNonTerminal = r0.length-1;
		    		  
		    		//  System.out.println("ii:"+i);
		    		//  System.out.println(nextNonTerminal + "MapToOriginalTS.get(nextNonTerminal) = "+mapToOriginalTS.get(nextNonTerminal));
		    		  if(isNumeric(r0[nextNonTerminal])) // negative
		    			  numEndPos = mapToOriginalTS.get(nextNonTerminal-nullCounter)-1;
		    		  else
		    			  numEndPos = mapToOriginalTS.get(nextNonTerminal-nullCounter);
		    		  
		    		  
		    		  currentAmountDistance = currentAmountDistance+allDiscordDistances.get(numStartPos+","+numEndPos);
		    		/*
		    		  System.out.print(cnt+": [");
		    		  cnt++;
		    		  for (int a = i; a<=nextNonTerminal; a++)
		    			  {
		    			  	System.out.print(" "+parseRule(r0[a]));
		    			  
		    			  }
		    		  System.out.println("]");
		    		  */
		    		//  System.out.println("i_nextNon : "+i+":"+nextNonTerminal+"["+numStartPos+"-"+numEndPos);
		    		//  numEndPos = mapToOriginalTS.get((i+minBlocks))-1;
		    		  RuleInterval ri = new RuleInterval(numStartPos,numEndPos);
			  	   	  anomalyIntervals.add(ri);
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
			    else
			    {
			    	this.accumulatDistance.put(Integer.valueOf(r0[i]), Double.valueOf(currentAmountDistance));
			    	currentAmountDistance = 0;
			    	i++; //Negative Number
			    }
		    	}
		    }
		    System.out.println("r0.length="+r0.length);
		    
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
		  		//				System.out.println("startPos: "+startPos);
		  		//				System.out.println("endPos: " +endPos);
		  						
		  					//	System.out.print("track#: "+counter+":       ");
		  						for (int j = startPos; j<=endPos; j++){
		  							
		  							Location loca = new Location(lat.get(j),lon.get(j));
		  				
		  							singleRoute.addLocation(lat.get(j), lon.get(j));
		  								
		  						
		  							
		  						}
		  						route.add(singleRoute);
		  						
		  						counter++;
		  				  }
		  		//		  System.out.println("position size: "+positions.size());
		  			//	  System.out.println("route size: "+route.size());
		  				
		  				  //  if(route.size()>2)
		  				     routes.add(route);
		
		  		  }	
		  	
		  			
		  			int startAnomalyPos = 0;
		  			int endAnomalyPos = 0;
		  			int anomalyCount1 = 0;
		  			
		  			
		  			
		  			
		  			for (int a = 0; a<isCovered.length;a++){
		  				if(!isCovered[a]){
		  					anomalyCount1++;
		  					//System.out.println("i: "+a+"\t block: "+blocks.findBlockIdForPoint(new Location(lat.get(a),lon.get(a))));
		  				}
		  				if(isCovered[a] && a<breakPoint)
		  					trueNegativeCount++;
		  				if(isCovered[a] && a>=breakPoint)
		  					falseNegativeCount++;
		  				if(!isCovered[a]&& a<breakPoint)
		  					falsePositiveCount++;
		  				if(!isCovered[a]&& a>=breakPoint)
		  					trueAnomalyCount++;
		  			}
		  			int i1 = 0;
		  			
		  			
		  			
		  			
		  			while  (i1<isCovered.length){
		  			//	System.out.println(i + " isCovered :"+isCovered[i]);
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
		  					  			//System.out.println("inner loop :"+i);
		  					  		}
		  					  		
		  					  	//	RuleInterval ri = new RuleInterval(startAnomalyPos,endAnomalyPos);
	  					  		//	anomalyIntervals.add(ri);
	  					  	//	System.out.println("new intervals :"+anomalyIntervals.size()+" : " + ri);
		  					  	
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
		  			  
		  			  
				
		  }

	private int getNextNonTerminal(int i) {
		int j = i+1;
	//	System.out.println("j: "+j);
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
		
		System.out.print("r0_"+i+"_"+(i+minBlocks-1)+": [");
		for (int j = i; j<r0.length&&j<(i+minBlocks);j++)
		{
			System.out.print(r0[j]+" ");
		}
		System.out.println("]");
		
		return true;
	}

	private Integer getPositionsInTS(ArrayList<Integer> mapToPreviousR0,ArrayList<Integer> previousMapToOriginalTS, int index) {
		
	//	if(previousMapToOriginalTS.get(mapToPreviousR0.get(index))==108)
		/*
				System.out.println("index = "+index+"    r0"+r0[index]);
				System.out.println(" mapToPreviousR0.get(index) ="+mapToPreviousR0.get(index));
				System.out.println("  previousMapToOriginalTS.get(mapToPreviousR0.get(index))  ="+previousMapToOriginalTS.get(mapToPreviousR0.get(index)));
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
		  System.out.println(EVALUATION_HEAD);
		  System.out.println(sb1);
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
		 // String header = "type,x,y";
//		  System.out.println("Total rules:"+chartData.getRulesNumber());
		  
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
		    		    	System.out.println("Error!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+ri);
		    		}
		    }
		    
		    	
		    for(int i = 0; i< clusters.size();i++){
		    	
		    	ArrayList<RuleInterval> mergedIntervals = new ArrayList<RuleInterval>();
		    	HashSet<Integer> set = new HashSet<Integer>();
		    	if(clusters.get(i).size()>0)
		    	{
		    		//System.out.println("cluster "+i+" : {" +clusters.get(i)+"}");
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
		    						System.out.println("I still need merge here.");
		    						System.out.println("1:" +mergedIntervals.get(k) +" 2:"+ newComer );
		    							
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
		    		System.out.println("mergedIntervals.size = "+mergedIntervals.size());
		    	}
		    }
		    System.out.println("Immergable Rule  = "+ immergableRuleCount);
		    System.out.println("Total Rule Count = "+totalRuleCount);
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
		//	  System.out.println("rule" + i+" :  "+ positions);//.get(0).toString());
			  
			  
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
		//				System.out.println("startPos: "+startPos);
		//				System.out.println("endPos: " +endPos);
						boolean firstPoint = true;
						
					//	System.out.print("track#: "+counter+":       ");
						for (int j = startPos; j<=endPos; j++){
							
				//			motifPos.append("T,"+lat.get(j)+","+lon.get(j));
						//	if(counter<2){
							Location loca = new Location(lat.get(j),lon.get(j));
						//	  blocks.addPoint2Block(loc);
					//		  Integer idss = new Integer(blocks.findBlockIdForPoint(loca));
					//		  System.out.print(idss+", ");
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
					//	System.out.println();
				  }
		//		  System.out.println("position size: "+positions.size());
			//	  System.out.println("route size: "+route.size());
				
				  //  if(route.size()>2)
				     routes.add(route);
				    
				    
				    
				    
				  //	motifPos.flush();
				  //	motifPos.close();
				  //	System.out.println(fname.getName());
			//	  	route.get(0).print();
				  	
				 
			  }
			  
			//  System.out.println("motif index: "+motif.getRuleIndex()+"   " +motif.toString());
			  //FileWriter motifPos = new FileWriter(new File("./motif_"+motif.getRuleIndex()+".csv"));
			 
		  	}
		  }	
			for (int i = 0;i<isCovered.length;i++){
				  if(isCovered[i]==true)
					  coverCount++;
			  }
			  System.out.println("Cover Count: "+ coverCount);
			  System.out.println("cover rate: " +(double)coverCount/isCovered.length);
		 
	}
*/

	public static void printArrayList(ArrayList<Integer> al) {
		if(al == null || al.size()==0)
			System.out.println("Null or empty ArrayList");
		else 
		{	
		//	System.out.print("[ ");
			for (int i = 0; i<al.size();i++)
				System.out.println(al.get(i)+" ");
			System.out.println();
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
		//			System.out.print(pairwiseDistances.get(x)+", ");
				}
		//	System.out.println();
		//	System.out.println("sum of pairwise distance: "+sums);
		//	System.out.println("pairSize = "+pairwiseDistances.size());
			double avgDistance = avg(pairwiseDistances);
			allDistances.add(avgDistance);
			
			Double stdDev = (Double)dev(pairwiseDistances);
			
			allStdDev.add(stdDev);
		/*	System.out.println("pairwire distances of motif "+i+": mean = "+avgDistance+",  Std.Dev ="+stdDev);
			for (int m = 0; m<pairwiseDistances.size();m++)
				System.out.print(" "+pairwiseDistances.get(m));
			System.out.println();
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
		
//		System.out.println("average distances among all motifs: "+avg(allDistances));
	//	System.out.println("average standard deviation among all motifs: "+avg(allStdDev));
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
	//		System.out.println("compare: "+allDistances.get(i)+"/"+allMinimalInterDistances.get(i)+" = "+sc);
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
			//	System.out.println("i="+i+" j="+j);
				double similarity = avgDTWDistance(eBlocks, allTracks.get(i),allTracks.get(j));
				pairwiseDistance.add(similarity);
			}
		}
		return pairwiseDistance;
	}

	private double avgDTWDistance(Blocks blocks, ArrayList<Integer> s,
			ArrayList<Integer> t) {
		
	//	System.out.print("s::::::::::size:"+s.size());
		/*
		for(int i=0; i<s.size();i++)
			System.out.print(" "+s.get(i));
			*/
	//	System.out.println();
	//	System.out.print("t::::::::::size:"+t.size()+"   ");
		/*
		for(int i=0; i<t.size();i++)
			System.out.print(" "+t.get(i));
			*/
	//	System.out.println();
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
			//	System.out.println("cost_"+i+","+j+": "+cost);
				DTW[i+1][j+1]=cost+minimum(DTW[i][j+1],		// insertion
										   DTW[i+1][j], 	// deletion
										   DTW[i][j]);	// match
			}
		}
	//	System.out.println("DTW:::::"+DTW[n][m]);
		int step = 1;
		int x = n;
		int y = m;
		while(!((x==1)&&(y==1))){
			step = step + 1;
			switch(min(DTW[x-1][y-1],DTW[x-1][y],DTW[x][y-1])){
			case 1: x--; y--; break;
			case 2: x--; break;
			case 3: y--; break;
			default: System.out.println("Error!!!!");
			}
			
		}
	//	System.out.println("step: "+step);
		double avg = DTW[n][m]/step;
	//	System.out.println("avgDTW:::::"+avg);
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
	  		//System.out.println("id: "+id);
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
		  	//	System.out.println("id   currentId:  "+id+"       "+currentId);	
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
	//		System.out.print(list.get(i)+" ");
		}
	//	System.out.println();
		
	//	System.out.println("sum = "+sum+" avg = "+sum/list.size());
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
    	System.out.println("getrawTrajectory: "+rawRoutes.size());
    	return rawRoutes; 
    }
    public static ArrayList<Route> getAnomaly(){
    	return anomalyRoutes; 
    }
    private void drawAnomaly() {
    	anomalyRoutes = new ArrayList<Route>();
		
		
			
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
  						System.out.println("["+startPos+"," +endPos+"]"+"Trajectory: "+getTrajectory(endPos)+"length = "+(endPos-startPos+1));
  						System.out.println();
  						
  					//	System.out.print("track#: "+counter+":       ");
  						
  						Location loca = new Location(lat.get(startPos),lon.get(startPos));
  						//Location endLoc = new Location(lat.get(startPos),lon.get(startPos));
  						
  						for (int j = startPos; j<=endPos; j++){
  							Location previousLoc =loca;
  							loca = new Location(lat.get(j),lon.get(j));
  							distance = distance + blocks.distance(blocks.findBlockIdForPoint(previousLoc), blocks.findBlockIdForPoint(loca)); 
  							singleRoute.addLocation(lat.get(j), lon.get(j));
  								
  						
  							
  						}
  					//	if (distance>0.1) // remove the false anomalies in the same block.
  						{
  						anomalyRoutes.add(singleRoute);
  						}
  				  }
  		//		  System.out.println("position size: "+positions.size());
  			//	  System.out.println("route size: "+route.size());
  				
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
	//    System.out.println("string: "+str+"   length = "+counter);
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
		//System.out.println("string: "+string);
		ArrayList<String> sa = new ArrayList<String>();
		String[] stringArray = string.split(" ");
		for (String s:stringArray){
			if (s.charAt(0)=='I')
			{
				if(s.contains("r")){
					int rIndex = s.indexOf("r");
					Integer iteration = Integer.valueOf(s.substring(1, rIndex));
					Integer rule = Integer.valueOf(s.substring(rIndex+1));
				//	System.out.println("s: "+s+" iteration: "+iteration+" rule: "+rule);
					String subRule = parseRule(allRules.get(iteration).get(rule).getExpandedRuleString());
					sa.add(subRule);
			//		System.out.println(s+" = "+subRule );
					
				}
				else if(s.contains("C")){
					int cIndex = s.indexOf("C");
					Integer iteration = Integer.valueOf(s.substring(1, cIndex));
					Integer cluster = Integer.valueOf(s.substring(cIndex+1));
				//	System.out.println("s: "+s+" iteration: "+iteration+" cluster: "+cluster);
					Integer ruleInCluster = (Integer)allClusters.get(iteration).get(cluster).toArray()[0];
					String subRule = parseRule(allRules.get(iteration).get(ruleInCluster).getExpandedRuleString());
					sa.add(subRule);
				//	System.out.println(s+" = "+subRule );
	
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
		//		System.out.println("s: "+ s);
			}
		}
		for (int i = 0; i<sa.size()-1;i++){
			sb.append(sa.get(i));
			sb.append(" ");
		}
		if(sa.size()>0)
		   sb.append(sa.get(sa.size()-1));
		//System.out.println("sb: "+sb.toString());
		String ans = sb.toString();
		return ans;
	}
	public static int getTrajectory(int endPos) {
		int i = endPos;
		Double traj;
		while(ncLat.get(i)>-999)
			i++;
		traj = -1000-ncLat.get(i);
		return traj.intValue();
	}
	public static <K, V extends Comparable<? super V>> Map<K, V> 
    sortByValue( Map<K, V> map )
{
    List<Map.Entry<K, V>> list =
        new LinkedList<>( map.entrySet() );
    Collections.sort( list, new Comparator<Map.Entry<K, V>>()
    {
        @Override
        public int compare( Map.Entry<K, V> o1, Map.Entry<K, V> o2 )
        {
            return (o2.getValue()).compareTo( o1.getValue() );
        }
    } );

    Map<K, V> result = new LinkedHashMap<>();
    for (Map.Entry<K, V> entry : list)
    {
        result.put( entry.getKey(), entry.getValue() );
    }
    return result;
}

	public static Instance trajToInstance(Integer traj) {
		// TODO Auto-generated method stub
		return null;
	}

}
