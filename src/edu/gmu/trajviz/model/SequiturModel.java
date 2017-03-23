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
import java.util.SortedSet;
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
import edu.gmu.trajviz.sax.datastructures.Cluster;
import edu.gmu.trajviz.sax.datastructures.Interval;
import edu.gmu.trajviz.sax.datastructures.Motif;
import edu.gmu.trajviz.sax.datastructures.SAXRecords;
import edu.gmu.trajviz.timeseries.TSException;
import edu.gmu.trajviz.util.StackTrace;
import edu.gmu.trajviz.util.Tools;
import net.sf.javaml.core.Instance;
public class SequiturModel extends Observable {
	
	/*
	 * 
	 */
//	public static ArrayList<String> allTrueAnomalyString;
	public static int count_works,not_works;  
	public static FileWriter frcsv;
	public static FileWriter frhead;
	public String dataset;
	public double heuristicF1Point;
	public double heuristicF1Overlap;
	public long bruteForceTime;
	public long heuristicTime;
	public static final Double OVERLAP_DEGREE = 0.3;
	private static final int STEP = 2;
	private static final int DEFAULT_TIME_GAP = 5;//180;
//	private static final int DEFAULT_TIME_GAP = 180;
//	private static final int DEFAULT_TIME_GAP = 3;
	
//=======================20170209=========================
	public static double stepDist;
	public ArrayList<Double> rescaleX;
	public ArrayList<Double> rescaleY;
	public static ArrayList<Integer[]> whole2separateTrajMap; // i.e. whole2separateTrajMap.get(indexOfrescaleX) =={trajId,position in rescaledRoutes}
	public int queryResultCounter;
	public int queryRuleCounter;
	public static ArrayList<Motif> motifList;
	public static ArrayList<Integer> motifListRuleMap;
	public static ArrayList<ArrayList<Motif>> trajMotifMap;
	public static HashMap<SortedSet<Integer>, ArrayList<Motif>> startPosMap2Motif;
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	public static int top1EuDistanceCalled;
	public static int allEuDistanceCalled;
	public static ArrayList<ArrayList<Double>> rawtrajX,rawtrajY,oldtrajX, oldtrajY, trajX, trajY,allTimeline;
	
	public static HashMap<Integer, ArrayList<ArrayList<Center>>> allCenterPoints;
	public static ArrayList<ArrayList<Integer>> rawsaxstrings, saxstrings, anomalyGroundTruth;
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
	public static ArrayList<ArrayList<Boolean>> isAnomaly,isLongAnomaly,isDiscord, isTrueAnomaly;
	public static double R;
	public static double rs;
	public static HashMap<String, HashSet<String>> motifMatches;
	public static HashMap<String, Cluster> subtrajClusterMap;
	public static int slidePoint;
	public static int coverPoint;
	public static HashMap<Integer, ArrayList<String>> anomalyMap;
	public static ArrayList<String> allAnomalies, allDiscords,allTrueAnomalyString, allResampledTrueAnomalyString;
	public static Map<Integer, ArrayList<String>> sortedAnomalyMap;
	/*
	 * 
	 */
	
	
//	public static double MINLINK = 0.0;
//	public final static double (maxPointErrorDistance*2) = 0.0;
	
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
	public static double maxPointErrorDistance;
	private static int resampleRate;
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
	private static ArrayList<Route> rawRoutes,rescaleRoutes;  
	private static ArrayList<Route> anomalyRoutes,discordRoutes,trueAnomalyRoutes;
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
	private ArrayList<Integer> words;
	public static Blocks blocks, eBlocks;
	public static int minBlocks; 
	private static Logger consoleLogger;
	private Map<Integer, Double> accumulatDistance;
	private static int identicalMotifStartPosSetCount;
	public static HashMap<String, Integer> subSeqBlockMap;
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
		  anomalyGroundTruth = new ArrayList<ArrayList<Integer>>();
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
				  anomalyGroundTruth.add(new ArrayList<Integer>());
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
						  anomalyGroundTruth.get(anomalyGroundTruth.size()-1).add(value2);
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
						  anomalyGroundTruth.add(new ArrayList<Integer>());
						  rawtrajX.get(rawtrajX.size()-1).add(value);
						  rawtrajY.get(rawtrajY.size()-1).add(value1);
						  anomalyGroundTruth.get(anomalyGroundTruth.size()-1).add(value2);
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
		  
		
			  
			
			  
		latMax = Double.valueOf(data.get(0));
		lonMax = Double.valueOf(data1.get(0));
		latMin = Double.valueOf(data.get(0));
		lonMin = Double.valueOf(data1.get(0));
		  for(int i = 0; i<data.size(); i++){
			  double temp_latitude = Double.valueOf(data.get(i));
			  double temp_longitude = Double.valueOf(data1.get(i));
			//  System.out.println("i = "+i+": "+temp_latitude+","+temp_longitude);
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
		  latOri = data;
		  lonOri = data1;
		  
		  data = new ArrayList<>();
		  data1 = new ArrayList<>();
		  lat_center = (latMax+latMin)/2;
		  lon_center = (lonMax+lonMin)/2;
		  diagnalDistance = Tools.pointEuDist(latMax, lonMax, latMin, lonMin);
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
	
	  
	  
	  
	  public synchronized void processData(double maxErrorSteps, int blockSize, int minContinualBlocks, int resamplingRate)throws IOException{
		  resampleRate = resamplingRate;
		 
		  minBlocks = minContinualBlocks;
		  
		  alphabetSize = blockSize;
		  stepDist = Tools.pointEuDist(latMin, lonMin, latMax, lonMax)/resampleRate;
		  maxPointErrorDistance = maxErrorSteps*stepDist;
		  
		  
		  R = Tools.pointEuDist(latMin, lonMin, latMax, lonMax)/minContinualBlocks; 
		  
		  initializeVariables();
		 
		  StringBuffer sb = new StringBuffer();
		  if (null == this.latOri ||null == this.lonOri|| this.latOri.size()==0 || this.lonOri.size()==0 ){
			  this.log("unable to \"Process data\" - no data were loaded...");
		  }
		  else{
			
			  consoleLogger.info("setting up GI with params: ");
			  sb.append(" algorithm: Sequitur");
			  sb.append(" maxErrorSteps: ").append(maxErrorSteps);
			  sb.append(" blockSize: ").append(blockSize);
			  sb.append(" Minimal Continuous Blocks: ").append(minContinualBlocks);
			  sb.append(" Resampling Rate: ").append(resamplingRate);
			  consoleLogger.info(sb.toString());
			 
			 
			  this.log(sb.toString());
		  }
		  long ct = System.currentTimeMillis();
		  buildModel();
		  System.out.println("buildModel time: "+ (System.currentTimeMillis()-ct));
		  ct = System.currentTimeMillis();
		  leftPanelRaw();
		  System.out.println("leftPanelRaw time: "+ (System.currentTimeMillis()-ct));
		  writeForMatlab();
		  ct = System.currentTimeMillis();
		  runSequitur();
		  System.out.println("runSequitur time: "+ (System.currentTimeMillis()-ct));
		  ct = System.currentTimeMillis();
		  postProcessing();
		  System.out.println("postProcessing time: "+ (System.currentTimeMillis()-ct));
		  ct = System.currentTimeMillis();
		  System.out.println("queryResultCounter = "+queryResultCounter);
		  System.out.println("queryRuleCounter = "+queryRuleCounter);
		/*  
		  for (int i = 0; i<motifList.size(); i++){
			  System.out.println("Motif "+i+" : "+motifList.get(i).trajIds);
			  System.out.println("Motif "+i+" : "+motifList.get(i).startPositions);
		  }
		  */
		  System.out.println("identicalMotifStartPosSetCount = "+this.identicalMotifStartPosSetCount);
		  
		  /*
		  long time = System.currentTimeMillis()/1000;
		  System.out.println("distCut = "+distCut);
			double stepDist = SequiturModel.distCut*minBlocks*maxPointErrorDistance;
			System.out.println("stepDist================="+stepDist);
			  SequiturModel.R = stepDist*stepDist;  	
		  int minLength = minBlocks;
		  int maxLength = minLength;
		 
		  System.out.println("Dataset = "+this.fileNameOnly);
		  System.out.println("MinLink = "+maxPointErrorDistance+", AlphabetSize = "+this.alphabetSize+", slidingWindow = "+this.minBlocks+ ", resampling rate = "+this.resolution);
		  System.out.println("findAllMotif time: "+this.heuristicTime);  
		  System.out.println("drawMotifsax time: "+(System.currentTimeMillis()-time));
		  

		//  System.out.println("R = "+SequiturModel.R);
		  time = System.currentTimeMillis();
		//  TrajDiscords.getThresholdDiscordsEvaluation((int)(minBlocks));
		  this.bruteForceTime = (System.currentTimeMillis()-time);
		  System.out.println("getThresholdDiscordsEvaluation time: "+this.bruteForceTime);
		  System.out.println("distCut = "+distCut);
		//	computeDummyConfusionMatrix();
		  computeResampledConfusionMatrix();
		  computeHitsGroundTruthConfusionMatrix();
		  computePointConfusionMatrix();
		  computeHitsConfusionMatrix();
		  System.out.println("allAnomalies = "+allAnomalies);
		  System.out.println("allDiscords  = "+allDiscords);
		  frcsv.append(this.minLink+","+this.alphabetSize+","+this.minBlocks+","+this.resolution+","+this.heuristicF1Overlap+","+this.heuristicF1Point+","+this.heuristicTime+","+this.bruteForceTime);
		  frcsv.close();
		  frhead.close();
		  System.out.println("Count works = "+count_works);
		  System.out.println("Not works = "+not_works);
		  System.out.println("Done!");
		  */
		  setChanged();
		  notifyObservers(new SequiturMessage(SequiturMessage.CHART_MESSAGE, motifList));
		
	  }
	  
	  
	/*
	 * 2017 new method  
	 */
	  
	private void writeForMatlab() {
		System.out.println("CSV style to matlab~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~``");
		for(int i = 0; i<rescaleRoutes.size(); i++){
		//	System.out.println("traj "+i+" : "+rescaleRoutes.get(i).size());
		//	System.out.println("r"+(i+1)+"x="+rescaleRoutes.get(i).getLats());
		//	System.out.println("r"+(i+1)+"y="+rescaleRoutes.get(i).getLons());
		}
	}

	private void postProcessing() {
		GrammarRuleRecord r0 = rules.get(0);
			
			System.out.println(r0);
			System.out.println("actualRuleYield :"+r0.getActualRuleYield());
			
			
			rules.sortByLength();
		//SortedSet<GrammarRuleRecord> set = rules.sortByValues();
		
		
		for (GrammarRuleRecord rule : rules.getSortedRulesByLength()){
		//	GrammarRuleRecord rule = rules.get(i);
			/*
			System.out.println(rule);
			System.out.println("ruleYield = " +rule.getRuleYield());
			System.out.println("expandedRuleString :"+rule.getExpandedRuleString());
			System.out.println("Rule String list:"+ rule.getRuleStringList());
		    */
			if(rule.getActualRuleYield()>=minBlocks){
			
			queryRuleCounter++;
			query(rule) ;
			}
			
		}
	}

	private void query(GrammarRuleRecord rule) {
		/* Todo 03/14/2017: 
		1. merge intervals among candidate grids
		2. pruning when compute euclidean distance using property 4. 
		*/
		
		
		ArrayList<String> stringList = rule.getRuleStringList();
		
		RuleInterval ruleInterval = rule.getRuleIntervals().get(0);   // get the 1st rule interval first
		int length = ruleInterval.getLength()+1;
		int queryStartPoint = ruleInterval.getStartPos();
		int queryEndPoint = ruleInterval.getEndPos();
		double maxSubtrajSquareEuDist = length*maxPointErrorDistance*maxPointErrorDistance;   //
		//Location queryStartLocation = new Location(rescale)
		
		/*=================================================================================================================
		 * Below is pruning 1: only select points in blocks which overlap with the query circle of query trajectory's start point  
		 */
		
		Block startBlock = blocks.findBlockById(words.get(queryStartPoint));
		Route queryRoute = new Route(rescaleX.subList(queryStartPoint,queryEndPoint+1),rescaleY.subList(queryStartPoint, queryEndPoint+1));
		Motif motif = new Motif(rule.getRuleNumber(), queryRoute,queryStartPoint);
		double[] lowerBoundDistance = startBlock.getLowerBoundDistance2Neighbor(rescaleX.get(queryStartPoint), rescaleY.get(queryStartPoint));
	//	blocks.printBlockMap();
	//	System.out.println("startBlock: "+startBlock);
		for(int i = -1; i<startBlock.nearbyBlocks.length; i++){
			Block nearbyBlock;
			if(i<0){
				nearbyBlock = startBlock;
			}
			else{
			 nearbyBlock = startBlock.nearbyBlocks[i];
			}
			if(i<0||(nearbyBlock!=null && lowerBoundDistance[i]<maxPointErrorDistance)){   // if this neighbor block might contains nearby points
              
       /*==================================================================================================================================
		* Below is pruning 2: only select sub-trajectories which are close to the start and end points located in query circle of
		*  query trajectory's start and end point respectively 
		* 
        */
				for(Interval interval : nearbyBlock.getIntervals()) {
	              int start = interval.getStartIdx();
	              int end = start+length;
	              /*
	              System.out.println(start +"   startTraj" +whole2separateTrajMap.get(start)[0]);	 
	              System.out.println(end + "     endTrajd" +whole2separateTrajMap.get(end)[0]);	 
	              System.out.println(start +" whole2separateTrajMp start==end : "+(whole2separateTrajMap.get(start)[0].equals(whole2separateTrajMap.get(end)[0]) ));
	             */
				  if(start<=interval.getEndIdx()&&
            			  end<rescaleX.size()&&whole2separateTrajMap.get(start)[0].equals(whole2separateTrajMap.get(end)[0])){	  
	            	  //int e = s+length;
	             
            	  Location startPoint = new Location(rescaleX.get(start),rescaleY.get(start));
            	  Location endPoint = new Location(rescaleX.get(end),rescaleY.get(end));
            	  double startDist = Tools.locationDist(startPoint,queryRoute.getStartLocation());
            	  double endDist = Tools.locationDist(endPoint,queryRoute.getEndLocation());

            	  while(start<=interval.getEndIdx()&&
            			  end<rescaleX.size()&&(whole2separateTrajMap.get(start)[0].equals(whole2separateTrajMap.get(end)[0]))
            			  &&(startDist>maxPointErrorDistance||endDist>maxPointErrorDistance)){
            		  int steps = Math.max( Double.valueOf((startDist-maxPointErrorDistance)/this.stepDist).intValue(),Double.valueOf((endDist-maxPointErrorDistance)/this.stepDist).intValue());
            		  steps = steps==0?1:steps;
            		  start = start+steps;
            		  end = end+steps;
            		  startDist = Tools.locationDist(startPoint,queryRoute.getStartLocation());
            		  endDist = Tools.locationDist(endPoint,queryRoute.getEndLocation());
            	  }
            	  /*======================================================================================================================================
            	  
         * Below is pruning 3: the closest sub-trajectory among all trivial sub-trajectory must satisfy the start and end positions must be also the closest pair.  
         */
            	  if(start<=interval.getEndIdx()&& end<rescaleX.size()&&(whole2separateTrajMap.get(start)[0].equals(whole2separateTrajMap.get(start+length)[0]))){
            	
                	//  System.out.println("end = "+end+"        sDist = "+startDist+ "eDist = "+endDist+"    maxPointErrorDistance = "+maxPointErrorDistance);
                	  double minDist = startDist+endDist;
                	  int minStart = start;
                	  
                	  while(start<=interval.getEndIdx()&& end<rescaleX.size()&&(whole2separateTrajMap.get(start)[0].equals(whole2separateTrajMap.get(start+length)[0]))&&
                			  startDist<maxPointErrorDistance&&endDist<maxPointErrorDistance){
                		  if(minDist<startDist+endDist){
                			  minDist = startDist+endDist;
                			  minStart = start;
                		  }
                		  start++;
                		  end++;
                		  startPoint = new Location(rescaleX.get(start),rescaleY.get(start));
                    	  endPoint = new Location(rescaleX.get(end),rescaleY.get(end));
                    	  startDist = Tools.locationDist(startPoint,queryRoute.getStartLocation());
                		  endDist = Tools.locationDist(endPoint,queryRoute.getEndLocation());
                		  
                	  }
                //	  System.out.println("start : end = "+start+" : "+end+ "  minStart  = "+minStart+"   length="+length);
            	  
            	  Route candidateRoute = new Route(rescaleX.subList(minStart, minStart+length),rescaleY.subList(minStart, minStart+length));
            	  double subtrajSquareEuDist = Tools.routeSqrEuDist(queryRoute, candidateRoute, maxSubtrajSquareEuDist);
            	  if(subtrajSquareEuDist<maxSubtrajSquareEuDist){
            		  motif.add(minStart, candidateRoute);
            		  queryResultCounter++;
            		 // System.out.println("found query Result:\n queryStartPoint = "+queryStartPoint+"\n resultStartPoint = "+minStart+"\n subtrajEuDist = "+subtrajSquareEuDist+" maxSubtrajSquareEudist = "+maxSubtrajSquareEuDist);
            	  }
            	  
            	  } // else advance to next interval.
            	  
            	  
                  }
			  }  // end for interval
			}
	//	  System.out.println(startBlock.nearbyBlocks[i]);
		}
		if(startPosMap2Motif.containsKey(motif.startPositions)){
			startPosMap2Motif.get(motif.startPositions).add(motif);
			identicalMotifStartPosSetCount++;
		}
		else{
		if(motif.size()>1){
			motifList.add(motif);
			motifListRuleMap.add(motif.id);
			ArrayList<Motif> tempList = new ArrayList<Motif>();
			tempList.add(motif);
			startPosMap2Motif.put(motif.startPositions, tempList);
		}
		}// end else
		
	//	System.out.println(queryRoute);
	}

	private void initializeVariables() throws IOException {
		  startPosMap2Motif = new HashMap<SortedSet<Integer>,ArrayList<Motif>>();
		  whole2separateTrajMap = new ArrayList<Integer[]>();  // i.e. whole2separateTrajMap.get(indexOfrescaleX) =={trajId,position in rescaledRoutes}
		  sortedCounter = 0;
		  count_works = 0;
		  not_works = 0;
		  identicalMotifStartPosSetCount = 0;
		  subSeqBlockMap = new HashMap<String, Integer>();
		  this.allTrajClusters = new HashMap<Integer, ArrayList<Cluster>>();
		  this.allMotifs = new HashMap<Integer, ArrayList<Cluster>>();
		  this.reTime = new ArrayList<Double>();
		  this.isLastIteration = false;
		 
		  this.allRules = new ArrayList<GrammarRules>();
		  this.allFilters = new ArrayList<ArrayList<Integer>>();
		  this.allClusters = new ArrayList<ArrayList<HashSet<Integer>>>();
		  this.allR0 = new ArrayList<String[]>();
		  this.allMapToPreviousR0 = new ArrayList<ArrayList<Integer>>();
		  this.allMapToOriginalTS = new ArrayList<ArrayList<Integer>>();
		  this.rawRoutes = new ArrayList<Route>();
		  this.rescaleRoutes = new ArrayList<Route>();
		  this.anomalyRoutes = new ArrayList<Route>();
		  this.discordRoutes = new ArrayList<Route>();
		  this.trueAnomalyRoutes = new ArrayList<Route>();
		  this.allTrajClusters = new HashMap<Integer, ArrayList<Cluster>>();
		  this.rawAllIntervals = new ArrayList<RuleInterval>();
		  this.subtrajClusterMap = new HashMap<String,Cluster>();
		  this.allAnomalies = new ArrayList<String>();
		  this.allDiscords = new ArrayList<String>();
		  this.allTrueAnomalyString = new ArrayList<String>();
		  String outputfn = this.fileNameOnly+"_"+this.alphabetSize+".csv";
		  frhead = new FileWriter(fileNameOnly+"_Result.csv");
		  frcsv = new FileWriter(outputfn);
		  this.dataset = this.fileNameOnly;
		  this.heuristicF1Point = -1.0;
		  this.heuristicF1Overlap = -1.0;
		  this.bruteForceTime = 0;
		  this.heuristicTime = 0;
		  
		  
	}

	private void computeHitsGroundTruthConfusionMatrix() {
		  	setTrueAnomalyString();
			int tt1 = 0; //  # of discords hit by anomaly
			int tt2 = 0;//   # of anomalies hit by discord
			int tf = 0;
			int ft = 0;
			int ff = 0;
			for(int i = 0; i<allResampledTrueAnomalyString.size(); i++ ){
				String trueAnomalyString = allResampledTrueAnomalyString.get(i);
				boolean hited = false;
				for(int j = 0; j<allAnomalies.size();j++){
					ArrayList<String> anomaly = new ArrayList<String>();
					String  a = allAnomalies.get(j);
					anomaly.add(a);
					
					while(j<allAnomalies.size()-1&&Tools.parseTrajId(allAnomalies.get(j+1))[0]==Tools.parseTrajId(a)[0]){
						j++;
					anomaly.add(allAnomalies.get(j));
					}
					if(isHit(anomaly,trueAnomalyString))
						{	
			//				System.out.println(anomaly+" hits discord "+discord);
							tt1++;
							hited = true;
							break;
						}
				}
				if(!hited)
					{
			//		System.out.println(" failed discord: "+discord);	
					tf++;
					}
			}
			for(int i = 0; i<allAnomalies.size(); i++ ){
				String anomaly = allAnomalies.get(i);
				boolean hited = false;
				for(int j = 0; j<allResampledTrueAnomalyString.size();j++){
					ArrayList<String> discord = new ArrayList<String>();
					String  d = allResampledTrueAnomalyString.get(j);
					discord.add(d);
					while(j<allResampledTrueAnomalyString.size()-1&&Tools.parseTrajId(allResampledTrueAnomalyString.get(j+1))[0]==Tools.parseTrajId(d)[0]){
						j++;
					discord.add(allResampledTrueAnomalyString.get(j));
					}
					if(isHit(discord,anomaly))
						{
				//			System.out.println(discord+" hits "+anomaly);
							tt1++;
							hited = true;
							break;
						}
				}
				if(!hited&&Tools.parseTrajId(anomaly)[2]>minBlocks)
					{
			//		System.out.println(" failed anomaly: "+anomaly);
						ft++;
					}
			}
			double tt = (tt1+tt2)/2.0;
			double precision = tt/(tt+ft);
			double recall = tt/(tt+tf);
			double f1 = 2*tt/(2*tt+tf+ft);
			double fmeasure = 2* (precision*recall)/(precision+recall);
			this.heuristicF1Overlap = f1;
			
			System.out.println("Resampled Hits Confusion Matris:");
			System.out.println("D/A\t True\t False");
			System.out.println("T\t"+tt+"\t"+tf);
			System.out.println("F\t"+ft+"\t"+ff);
			System.out.println("precision = "+precision );
			System.out.println("recall = "+recall );
			System.out.println("f1 = "+f1 );
			System.out.println("fmeasure = "+fmeasure );
			
				
			}

private void computeHitsConfusionMatrix() {
	int tt1 = 0; //  # of discords hit by anomaly
	int tt2 = 0;//   # of anomalies hit by discord
	int tf = 0;
	int ft = 0;
	int ff = 0;
	for(int i = 0; i<allDiscords.size(); i++ ){
		String discord = allDiscords.get(i);
		boolean hited = false;
		for(int j = 0; j<allAnomalies.size();j++){
			ArrayList<String> anomaly = new ArrayList<String>();
			String  a = allAnomalies.get(j);
			anomaly.add(a);
			
			while(j<allAnomalies.size()-1&&Tools.parseTrajId(allAnomalies.get(j+1))[0]==Tools.parseTrajId(a)[0]){
				j++;
			anomaly.add(allAnomalies.get(j));
			}
			if(isHit(anomaly,discord))
				{	
	//				System.out.println(anomaly+" hits discord "+discord);
					tt1++;
					hited = true;
					break;
				}
		}
		if(!hited)
			{
	//		System.out.println(" failed discord: "+discord);	
			tf++;
			}
	}
	for(int i = 0; i<allAnomalies.size(); i++ ){
		String anomaly = allAnomalies.get(i);
		boolean hited = false;
		for(int j = 0; j<allDiscords.size();j++){
			ArrayList<String> discord = new ArrayList<String>();
			String  d = allDiscords.get(j);
			discord.add(d);
			while(j<allDiscords.size()-1&&Tools.parseTrajId(allDiscords.get(j+1))[0]==Tools.parseTrajId(d)[0]){
				j++;
			discord.add(allDiscords.get(j));
			}
			if(isHit(discord,anomaly))
				{
		//			System.out.println(discord+" hits "+anomaly);
					tt1++;
					hited = true;
					break;
				}
		}
		if(!hited&&Tools.parseTrajId(anomaly)[2]>minBlocks)
			{
	//		System.out.println(" failed anomaly: "+anomaly);
				ft++;
			}
	}
	double tt = (tt1+tt2)/2.0;
	double precision = tt/(tt+ft);
	double recall = tt/(tt+tf);
	double f1 = 2*tt/(2*tt+tf+ft);
	double fmeasure = 2* (precision*recall)/(precision+recall);
	this.heuristicF1Overlap = f1;
	System.out.println("D/A\t True\t False");
	System.out.println("T\t"+tt+"\t"+tf);
	System.out.println("F\t"+ft+"\t"+ff);
	System.out.println("precision = "+precision );
	System.out.println("recall = "+recall );
	System.out.println("f1 = "+f1 );
	System.out.println("fmeasure = "+fmeasure );
	
		
	}

private boolean isHit(ArrayList<String> x, String y) {
	
//		
	int alength = 0;
	int overlap = 0;
	int[] d = Tools.parseTrajId(y);
	int ds = d[1];
	int de = ds+d[2]-1;
	for(int i = 0; i<x.size();i++){
		String anomaly = x.get(i);
		
	int[] a = Tools.parseTrajId(anomaly);
	if(a[0]!=d[0])
		return false;
	int as = a[1];
	alength = alength+a[2];
	int ae = as+a[2]-1;
	
	for(int c = as; c<=ae; c++){
		if(c>=ds&&c<=de)
			overlap++;
	}
	
		
	}
	int length = Math.min(alength, d[2]);
	double degree = ((double)overlap)/length;
	if(degree>1||degree<0)
		throw new IllegalArgumentException("overlap = "+overlap+"  length = "+length);
	if(degree>=OVERLAP_DEGREE)
		return true;
	else
		return false;
	/*
	if(a[0]!=d[0])
		return false;
		 if(as>=ds && as<de)
			overlap = Math.min(ae,de)-as+1; 
	else if(ae>ds && ae<=de)
			overlap = ae - Math.max(as,ds)+1; 
	else if(ds>=as && ds<ae)
			overlap = Math.min(de,de)-ds+1; 
	else if(de>as && de<=ae)
			overlap = de - Math.max(as,ds)+1; 
	else
		return false;
	if((double)overlap/length>=OVERLAP_DEGREE)
		return true;
	else
		return false;
	*/
	
}
private void computeDummyConfusionMatrix(){
	int tt = 0;
	int tf = 0;
	int ft = 0;
	int ff = 0;
	for(int traj = 0; traj<isLongAnomaly.size(); traj++){
		System.out.println(traj +" : "+isLongAnomaly.get(traj));
		for(int s = 0; s<isLongAnomaly.get(traj).size();s++){
			boolean isdiscord = (anomalyGroundTruth.get(traj).get(s)==0);
			boolean isanomaly = isLongAnomaly.get(traj).get(s);
			
			if(isdiscord==true&&isanomaly==true)
				tt++;
			if(isdiscord==true&&isanomaly==false)
				tf++;
			if(isdiscord==false&&isanomaly==true)
				ft++;
			if(isdiscord==false&&isanomaly==false)
				ff++;
		}
	}
	System.out.println("Dummy Confusion Matrix:");
	System.out.println("D/A\t True\t False");
	System.out.println("T\t"+tt+"\t"+tf);
	System.out.println("F\t"+ft+"\t"+ff);
	double accuracy = ((double)tt+ff)/(tt+ff+tf+ft);
	double precision = (double)tt/(tt+ft);
	double recall = (double)tt/(tt+tf);
	double f1 = 2*(double)tt/(2*tt+tf+ft);
	double fmeasure = 2* (precision*recall)/(precision+recall);
	this.heuristicF1Point = f1;
	System.out.println("accuracy = "+accuracy );
	System.out.println("precision = "+precision );
	System.out.println("recall = "+recall );
	System.out.println("f1 = "+f1 );
	System.out.println("fmeasure = "+fmeasure );
}
private void computeResampledConfusionMatrix(){
	int tt = 0;
	int tf = 0;
	int ft = 0;
	int ff = 0;
	for(int traj = 0; traj<isLongAnomaly.size(); traj++){
		System.out.println(traj +" : "+isLongAnomaly.get(traj));
		for(int s = 0; s<isLongAnomaly.get(traj).size();s++){
			boolean istrueanomaly = isTrueAnomaly.get(traj).get(s);
			boolean isanomaly = isLongAnomaly.get(traj).get(s);
			
			if(istrueanomaly==true&&isanomaly==true)
				tt++;
			if(istrueanomaly==true&&isanomaly==false)
				tf++;
			if(istrueanomaly==false&&isanomaly==true)
				ft++;
			if(istrueanomaly==false&&isanomaly==false)
				ff++;
		}
	}
	System.out.println("Resampled Anomaly Confusion Matrix:");
	System.out.println("D/A\t True\t False");
	System.out.println("T\t"+tt+"\t"+tf);
	System.out.println("F\t"+ft+"\t"+ff);
	double accuracy = ((double)tt+ff)/(tt+ff+tf+ft);
	double precision = (double)tt/(tt+ft);
	double recall = (double)tt/(tt+tf);
	double f1 = 2*(double)tt/(2*tt+tf+ft);
	double fmeasure = 2* (precision*recall)/(precision+recall);
	this.heuristicF1Point = f1;
	System.out.println("accuracy = "+accuracy );
	System.out.println("precision = "+precision );
	System.out.println("recall = "+recall );
	System.out.println("f1 = "+f1 );
	System.out.println("fmeasure = "+fmeasure );
}
private void computePointConfusionMatrix() {
	int tt = 0;
	int tf = 0;
	int ft = 0;
	int ff = 0;
	for(int traj = 0; traj<isLongAnomaly.size(); traj++)
		for(int s = 0; s<isLongAnomaly.get(traj).size();s++){
			boolean isdiscord = isDiscord.get(traj).get(s);
			boolean isanomaly = isLongAnomaly.get(traj).get(s);
			
			if(isdiscord==true&&isanomaly==true)
				tt++;
			if(isdiscord==true&&isanomaly==false)
				tf++;
			if(isdiscord==false&&isanomaly==true)
				ft++;
			if(isdiscord==false&&isanomaly==false)
				ff++;
		}
	System.out.println("D/A\t True\t False");
	System.out.println("T\t"+tt+"\t"+tf);
	System.out.println("F\t"+ft+"\t"+ff);
	double accuracy = ((double)tt+ff)/(tt+ff+tf+ft);
	double precision = (double)tt/(tt+ft);
	double recall = (double)tt/(tt+tf);
	double f1 = 2*(double)tt/(2*tt+tf+ft);
	double fmeasure = 2* (precision*recall)/(precision+recall);
	this.heuristicF1Point = f1;
	System.out.println("accuracy = "+accuracy );
	System.out.println("precision = "+precision );
	System.out.println("recall = "+recall );
	System.out.println("f1 = "+f1 );
	System.out.println("fmeasure = "+fmeasure );

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
	for(int i = 0; i<blocks.blocks.size();i++){
		blocks.blocks.get(i).setResidualSet();
	}
	for(int i = 0; i<blocks.blocks.size(); i++){
		//blocks.blocks.get(i).findHierarchicalMotifs();
		blocks.blocks.get(i).findAnomaly();
	}
	for(int i = 0; i<blocks.blocks.size(); i++){
		
		blocks.blocks.get(i).checkNearby();;
	}

	
}
	
	 


	private ArrayList<ArrayList<Center>> getCenterArrayList(int len) {
		ArrayList<ArrayList<Center>> allCenters = new ArrayList<ArrayList<Center>>();
		for(int traj = 0; traj<oldtrajX.size(); traj++)
		{	
			
			ArrayList<Center> centers = new ArrayList<Center>();
		//	ArrayList<Double> centerX = new ArrayList<Double>();
		//	ArrayList<Double> centerY = new ArrayList<Double>();
			ArrayList<Double> currentX = oldtrajX.get(traj);
			ArrayList<Double> currentY = oldtrajY.get(traj);
			
			
		//	System.out.println("debug traj = "+traj);
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
	public static void drawGroundTruth()	
	{
	for(int traj = 0; traj<anomalyGroundTruth.size();traj++){
		System.out.println("anomalyGroundTruth: "+anomalyGroundTruth.get(traj));
		  int pos = 0;
		  int endPos = 0;
		  int startPos = pos;
		  while(pos<anomalyGroundTruth.get(traj).size()){
			  if(anomalyGroundTruth.get(traj).get(pos)==0||(pos+1<anomalyGroundTruth.get(traj).size()&&anomalyGroundTruth.get(traj).get(pos+1)==0))    // is the start position of anomaly
			  {
				  startPos = pos;
				  int prePos =pos;
				  Route singleTrueAnomaly = new Route();
				  double x = rawtrajX.get(traj).get(pos);
				  double y = rawtrajY.get(traj).get(pos);
				  singleTrueAnomaly.addLocation(x, y);
				  pos++;
				  while(pos<anomalyGroundTruth.get(traj).size()&&(anomalyGroundTruth.get(traj).get(pos)==0||(anomalyGroundTruth.get(traj).get(pos-1)==0))){
					 /*
					  if(euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)>distCut*10){
						  System.out.println(rawtrajX.get(traj));
						  System.out.println(oldtrajX.get(traj));
						  System.out.println(rawtrajY.get(traj));
						  System.out.println(oldtrajY.get(traj));
						  
						  throw new IllegalArgumentException("traj = "+traj+"    pos = "+ pos+"    dist = "+euDist(oldtrajX.get(traj).get(pos),oldtrajY.get(traj).get(pos),x,y)+"prevous pos/distCut: "+distCut);
					  }
					  */
					  x =  rawtrajX.get(traj).get(pos);
					  y = rawtrajY.get(traj).get(pos);
					  singleTrueAnomaly.addLocation(x, y);
					  
					  prePos = pos;
					  pos++;
				  }
				 
				  if(singleTrueAnomaly.getLats().size()>=3)//(int)(minBlocks*SequiturModel.maxPointErrorDistance))
				  {
					  
					  Integer l = singleTrueAnomaly.getLats().size();
				  String DiscordString = "T"+traj+"S"+startPos+"L"+l;
					  allTrueAnomalyString.add(DiscordString);
					  
					  trueAnomalyRoutes.add(singleTrueAnomaly);
					  
					  
					  System.out.println("Anomalous Trajectory: "+traj+"-"+startPos+"-"+pos);
				  }
			  }
			  else
				  pos++;
		  }
	  }
}	
	
	//helper method to set allTrueAnomalyString according to is TrueAnomaly;
	public static void setTrueAnomalyString()	
	{
		allResampledTrueAnomalyString = new ArrayList<String>();
	for(int traj = 0; traj<isTrueAnomaly.size();traj++){
		
		  int pos = 0;
		  int endPos = 0;
		  int startPos = pos;
		  while(pos<isTrueAnomaly.get(traj).size()){
			  if(isTrueAnomaly.get(traj).get(pos))    // is the start position of anomaly
			  {
				  startPos = pos;
				  int prePos =pos;
				  Route singleDiscord = new Route();
				  double x = oldtrajX.get(traj).get(pos);
				  double y = oldtrajY.get(traj).get(pos);
				  singleDiscord.addLocation(x, y);
				  pos++;
				  while(pos<isTrueAnomaly.get(traj).size()&&isTrueAnomaly.get(traj).get(pos)){
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
					  singleDiscord.addLocation(x, y);
					  
					  prePos = pos;
					  pos++;
				  }
				 
				  if(singleDiscord.getLats().size()>=3)//(int)(minBlocks*SequiturModel.maxPointErrorDistance))
				  {
					  
					  Integer l = singleDiscord.getLats().size();
				  String resampledTrueAnomalyString = "T"+traj+"S"+startPos+"L"+l;
				  	  allResampledTrueAnomalyString.add(resampledTrueAnomalyString);
//					  allDiscords.add(DiscordString);
//					  discordRoutes.add(singleDiscord);
					  
					  
			//		  System.out.println("Anomalous Trajectory: "+traj+"-"+startPos+"-"+pos);
				  }
			  }
			  else
				  pos++;
		  }
	  }
}
public static void DrawDiscordOnly()	
	{
	for(int traj = 0; traj<isDiscord.size();traj++){
		  int pos = 0;
		  int endPos = 0;
		  int startPos = pos;
		  while(pos<isDiscord.get(traj).size()){
			  if(isDiscord.get(traj).get(pos))    // is the start position of anomaly
			  {
				  startPos = pos;
				  int prePos =pos;
				  Route singleDiscord = new Route();
				  double x = oldtrajX.get(traj).get(pos);
				  double y = oldtrajY.get(traj).get(pos);
				  singleDiscord.addLocation(x, y);
				  pos++;
				  while(pos<isDiscord.get(traj).size()&&isDiscord.get(traj).get(pos)){
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
					  singleDiscord.addLocation(x, y);
					  
					  prePos = pos;
					  pos++;
				  }
				 
				  if(singleDiscord.getLats().size()>=3)//(int)(minBlocks*SequiturModel.maxPointErrorDistance))
				  {
					  
					  Integer l = singleDiscord.getLats().size();
				  String DiscordString = "T"+traj+"S"+startPos+"L"+l;
					  allDiscords.add(DiscordString);
					  discordRoutes.add(singleDiscord);
					  
					  
			//		  System.out.println("Anomalous Trajectory: "+traj+"-"+startPos+"-"+pos);
				  }
			  }
			  else
				  pos++;
		  }
	  }
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
						  if(singleAnomaly.getLats().size()>=3)//(int)(minBlocks*SequiturModel.maxPointErrorDistance))
						  {
							  Integer l = singleAnomaly.getLats().size();
					//		  System.out.println("DEBUG::::::::::::::: l = "+l);
							  for(int i = 0; i<l;i++){
								  isLongAnomaly.get(traj).set(prePos-i, true);
							  }
						//	  System.out.println("DEBUG: isLongAnomalySetTrue!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
							//  if(l>=minBlocks*0.5)
								  anomalyRoutes.add(singleAnomaly);
							  String anomalyString = "T"+traj+"S"+startPos+"L"+l;
							  allAnomalies.add(anomalyString);
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
	/*
	private void addToBestTrajClusters(int traj, int s, int length){
		boolean isAdded = false;
		double minDist = Double.MAX_VALUE;
		int minCluster = 0;
		for(int i = 0; i<allTrajClusters.get(length).size(); i++){
			double dist = 0;
			for(int index = 0; index<allTrajClusters.get(length).get(i).repLineX.size(); index++){
				double pairDist = Tools.euDist(allTrajClusters.get(length).get(i).repLineX.get(index),allTrajClusters.get(length).get(i).repLineY.get(index),oldtrajX.get(traj).get(s+index),oldtrajY.get(traj).get(s+index));
				if(pairDist>distCut*(length)*maxPointErrorDistance)
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
		
		if((minDist/length)<=(distCut*(length)*maxPointErrorDistance)){
			//isAdded = true;
		System.out.println("trajId = "+traj);
		System.out.println("minDist/length = "+minDist/length);
		System.out.println("R      = "+distCut*(length)*maxPointErrorDistance);
		System.out.println("minCluster = "+ minCluster);
			allTrajClusters.get(length).get(minCluster).addBest(traj, s);
		}
		else{
			Cluster cluster = new Cluster(length);
			cluster.addBest(traj, s);
			allTrajClusters.get(length).add(cluster);
		}
		
	}
	*/

	/*
	private boolean isClose(int i, int j, int bestS1, int bestS2, int length) {
		
		for(int index = 0; index<length; index++){
			if (Tools.euDist(trajX.get(i).get(bestS1+index),trajY.get(i).get(bestS1+index),trajX.get(j).get(bestS2+index),trajY.get(j).get(bestS2+index))>this.minLink)
				return false;
			
		}	
		return true;
	}
	*/
	
private void paa2saxseqs() {
	      saxseqs =  new ArrayList<ArrayList<Center>>();
	//	  rawsaxstrings = new ArrayList<ArrayList<Integer>>();
	//	  saxstrings = new ArrayList<ArrayList<Integer>>();
		  blocks = new Blocks(alphabetSize,latMin,latMax,lonMin,lonMax);
		  double latCut = blocks.latCut;
		  double lonCut = blocks.lonCut;
		
		  words = new ArrayList<Integer>();
		  
	
		  
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
			//	  blocks.addCenter2Block(center);
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
		System.out.println("BBBBBBBBBBBBBBBBBBBBBBBBBbb  latMin = "+latMin+"   lonMin = "+lonMin +"   latMax = "+latMax+"  lonMax = "+lonMax );
		blocks = new Blocks(alphabetSize,latMin,latMax,lonMin,lonMax);
		queryResultCounter = 0;
		queryRuleCounter = 0;
		System.out.println("stepDist = "+stepDist);
		
		resample();
		descritization();
		
}
		 		
	private void descritization() {
		  words = new ArrayList<Integer>();
		  // add all points into blocks.
		  Integer previousId=(Integer)(-1);
		  Integer startIndex = 0;
		
		//  HashMap<Integer,Integer> trackMap = new HashMap<Integer, Integer>();
		  trimedTrack = new ArrayList<NumerosityReductionMapEntry>();
		  mapTrimed2Original = new ArrayList<Integer>();  // The index is the position in trimmed array, and the content is position in original time series.
		  mapToOriginalTS = new ArrayList<Integer>();
		  int startPoint = 0;
		  int endPoint = 0;
		  
		  /*
		   * Numerosity Reduction
		   */
		  for (int i = 0; i<rescaleX.size();i++){
			  Integer id;
			  if(rescaleX.get(i)<=-1000){
				  id = rescaleX.get(i).intValue();
			  }
			  else{
				  Location loc = new Location(rescaleX.get(i),rescaleY.get(i));
				  id = new Integer(blocks.findBlockIdForPoint(loc));
				  //blocks.addPoint2Block(loc, id);
			  }
			 
			words.add(id);
		 
		  
			
		//	  System.out.println("previousId, id:  "+previousId+",   "+id+"        i:   "+i+"   Lat,Lon: "+paaLat.get(i)+","+paaLon.get(i));
			  Integer trimedIndex = 0;
			  if (!id.equals(previousId))
			  {
				  
				  NumerosityReductionMapEntry<Integer, String> 
				  entry = new NumerosityReductionMapEntry<Integer, String>(i,id.toString());
		//		  System.out.println("entry: "+i+","+id);
				  trimedTrack.add(entry);
				  mapTrimed2Original.add(i);
				  mapToOriginalTS.add(i);
				  //put the new <index,id> pair into a map 
				  //NumerosityReductionMapEntry entry = new NumerosityReductionMapEntry(i,id);
			//	  trackMap.put(i, id);
				  if(previousId>0){
				  blocks.addInterval2Block(startIndex,i-1,previousId);
				//  System.out.println(previousId+" "+blocks.findBlockById(previousId).getIntervals());
				  }
				  previousId = id;
				  startIndex = i>0?i-1:0;
			  }
			  
			  				  
		  }
		  
		  /*
		  System.out.println("words = "+words);
		  
		  System.out.print("mapTrimed2Original: ");
		  for(int i=0; i<trimedTrack.size();i++){
			  System.out.print(i+" : "+trimedTrack.get(i).getKey()+":"+trimedTrack.get(i).getValue()+"\t");
		  }
		  
		  System.out.println();
		  */
	

	}
    private void dummyResampling(){
    	R = maxPointErrorDistance;
  	  oldtrajX = new ArrayList<ArrayList<Double>>();
		  oldtrajY = new ArrayList<ArrayList<Double>>();
		  travelDistance = new ArrayList<Double>();
		  isAnomaly = new ArrayList<ArrayList<Boolean>>();
		  isLongAnomaly = new ArrayList<ArrayList<Boolean>>();
		  isDiscord = new ArrayList<ArrayList<Boolean>>();
		  allTimeline = new ArrayList<ArrayList<Double>>();
		  ArrayList<ArrayList<Double>> rawTimeResample = new ArrayList<ArrayList<Double>>(); //  present reTime in ArrayList
		  double amountDist = 0;  
		 			ArrayList<ArrayList<Double>> timeResample = new ArrayList<ArrayList<Double>>();
			for(int i = 0; i<rawtrajX.size();i++){
				
				timeResample.add(new ArrayList<Double>());
				oldtrajX.add(new ArrayList<Double>());
				oldtrajY.add(new ArrayList<Double>());
				isAnomaly.add(new ArrayList<Boolean>());
				isLongAnomaly.add(new ArrayList<Boolean>());
				isDiscord.add(new ArrayList<Boolean>());
				
				double time = 0;
				int index = 0;
				
	//			System.out.println("i===================================================================================="+i);
				//while(time<rawTimeResample.get(i).get(rawTimeResample.get(i).size()-1)){
					while(index<rawtrajX.get(i).size()){
					timeResample.get(i).add(time);
					oldtrajX.get(i).add(rawtrajX.get(i).get(index));
					oldtrajY.get(i).add(rawtrajY.get(i).get(index));
					isAnomaly.get(i).add(true);
					isLongAnomaly.get(i).add(false);
					isDiscord.get(i).add(true);
					index++;
				
				}	
			}
			
			
			
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
      /* 20170213
	 * resample original ts with the same speed.
	 */
	  private  void resample() {
		  rescaleX = new ArrayList<Double>();
		  rescaleY = new ArrayList<Double>();
	//	  rescaleX.add(latOri.get(0));
	//	  rescaleY.add(lonOri.get(0));
		//  System.out.println(latOri.subList(0, 100));
		//  System.out.println(lonOri.subList(0, 100));
		  double x=Double.MAX_VALUE;
		  double y=Double.MAX_VALUE;
		  int trajId = -1001;
		  int posInTraj = 0;
		  for(int j = 1; j<latOri.size(); j++){
			 if(latOri.get(j)<=-1000){    // if the current point is the trajectory id
				 
				 rescaleX.add(latOri.get(j));
				 rescaleY.add(lonOri.get(j));
				
				 
				 posInTraj = 0;
				 Integer[] map = {trajId,-1};
				 whole2separateTrajMap.add(map);
				 trajId--;// = latOri.get(j).intValue();
				/*
				 if(j+1<latOri.size()){  // add the start point to rescaled trajectory
					rescaleX.add(latOri.get(j+1));
					rescaleY.add(lonOri.get(j+1));
				}
				*/
				x = Double.MAX_VALUE;
				y = Double.MAX_VALUE;
				j++;
				 continue;  
			 }
			 
			 
			 double x0 = latOri.get(j-1);
			 double y0 = lonOri.get(j-1);
			 double x1 = latOri.get(j);
			 double y1 = lonOri.get(j);
			 double l;
			 if(x!=Double.MAX_VALUE){   // adjust the x0, y0
				 double dist=Tools.pointEuDist(x, y, x1, y1);
				
				 if(dist<stepDist)   //Figure 4 condition
					 continue;
				 if(Math.abs(x1-x0)<0.000001){  // handling a = infinity case
					 //x'==x0; so (x'-x)^2+(y'-y)^2 = stepDist*stepDist  <=> (x0-x)^2+(y'-y)^2 = stepDist*stepDist;
					 double A = 1;
					 double B = -2*y;
					 double C = (x0-x)*(x0-x)+y*y-stepDist*stepDist;
					 double b4ac=Math.sqrt(B*B-4*A*C);
					 double root1 = (-B+b4ac)/(2*A);   
					 double root2 = (-B-b4ac)/(2*A);
					 if(root1>=Math.min(y0, y1)&&root1<=Math.max(y0,y1)){
						 y0 = root1;
					 }
					 else if(root2>=Math.min(y0, y1)&&root2<=Math.max(y0,y1)){
						 y0 = root2;
					 }
					 else{
						 double test_dist = Tools.pointEuDist(x, y, x0, y0);
						 
						 System.out.println("stepDist = "+stepDist);
						 System.out.println("x_y1_dist = "+dist);
						 System.out.println("x_y0_test_dist = "+test_dist);
						 System.out.println("x = "+x);
						 System.out.println("y = "+y);
						 
						 System.out.println("A:B:C = "+A+" "+B+" "+C);
						 System.out.println("B^2-4AC = "+(B*B-4*A*C));
//						 System.out.println("sqrt_bb4ac = "+sqrt_bb4ac);
						System.out.println(" y computation error::::::::  y0="+y0+"    y1="+y1+"   root1="+root1+"   root2="+ root2);
					
					 }
					 
				 }
				 else{
				 double a = (y1-y0)/(x1-x0);
				 double b = y0-a*x0;
				 // y' = ax_1'+b;
				 //Math.sqrt((x-x')^2+(y-y')^2) = stepDist;
				 // solve the above equations about x' and y'
				 double A = 1+a*a;
				 double B = 2*a*b-2*x-2*a*y;
				 double C = x*x+y*y+b*b-2*y*b-stepDist*stepDist;
		
				 double sqrt_bb4ac = Math.sqrt(B*B-4*A*C);
		
				 double root1 = (sqrt_bb4ac-B)/(2*A);
				 double root2 = (-sqrt_bb4ac-B)/(2*A);
				 if(root1-Math.min(x0, x1)>-0.000001&&root1-Math.max(x0,x1)<0.000001){
					 x0 = root1;
				 }
				 else if(root2-Math.min(x0, x1)>-0.000001&&root2-Math.max(x0,x1)<0.000001){
					 x0 = root2;
				 }
				 else{
					 
					 double test_dist = Tools.pointEuDist(x, y, x0, y0);
					 
					 System.out.println("stepDist = "+stepDist);
					 System.out.println("dist = "+dist);
					 System.out.println("test_dist = "+test_dist);
					 System.out.println("x = "+x);
					 System.out.println("y = "+y);
					 
					 System.out.println("A:B:C = "+A+" "+B+" "+C);
					 System.out.println("B^2-4AC = "+(B*B-4*A*C));
//					 System.out.println("sqrt_bb4ac = "+sqrt_bb4ac);
					System.out.println("computation error::::::::  x0="+x0+"    x1="+x1+"   root1="+root1+"   root2="+ root2);
				 }
				 y0 = a*x0+b;
				 }
			 }
			 x = x0;
			 y = y0;
			 rescaleX.add(x0);
			 rescaleY.add(y0);
			 Integer[] map = {trajId,posInTraj};
			 posInTraj++;
			 whole2separateTrajMap.add(map);
				 l = Tools.pointEuDist(x0, y0, x1, y1);
			 for(int i=1; i*stepDist<=l+0.00001; i++){
				 x = i*stepDist*(x1-x0)/l+x0;
				 y = i*stepDist*(y1-y0)/l+y0;
				 rescaleX.add(x);
				 rescaleY.add(y);
				 Integer[] map1 = {trajId,posInTraj};
				 posInTraj++;
				 whole2separateTrajMap.add(map1);
			 }
			  
			  
		  }
		  
		  System.out.println("latOri = "+latOri.subList(0, 100));
		  System.out.println("lonOri = "+lonOri.subList(0, 100));
		  System.out.println("stepDist = "+stepDist);
		  System.out.println("original length : rescaled length = "+latOri.size()+" : "+rescaleX.size());
		  System.out.println("whole2separateTrajMap.size() = "+whole2separateTrajMap.size());
		  System.out.println("rescaleX = "+rescaleX.subList(0, 100));
		  System.out.println("rescaleY = "+rescaleY.subList(0, 100));
		  /*
		  for(int i = 0; i<whole2separateTrajMap.size(); i++)
		  System.out.println(i+"("+whole2separateTrajMap.get(i)[0]+","+whole2separateTrajMap.get(i)[1]+")  ");
		  */
		  
		  
		  
		
	}


	

	// Draw raw trajectory on left map
	private void leftPanelRaw(){
		System.out.println("rawtrajX.size = "+rawtrajX.size());
		/*
		for (int i =0; i<rawtrajX.size();i++){
			if(rawtrajX.get(i).size()==0||rawtrajX.get(i).size()==0)
			System.out.println(rawtrajX.get(i));
			Route singleRoute = new Route(rawtrajX.get(i),rawtrajY.get(i));
		
			rawRoutes.add(singleRoute);
			
		}*/
		int start = 0;
		int trajId = -1001;
		int posInTraj = 0;
		for(int i = 0; i<rescaleX.size();i++){
		//	System.out.println("i = "+i+"  x = "+rescaleX.get(i));
			if(rescaleX.get(i)<=-1000&&i-start>1){
				/*
				if(rescaleY.get(i-1)<-500)
				System.out.println("start = "+ start+"i = "+i+ "rescaleY.get(i-1) = "+rescaleY.get(i-1));
				*/
				
				List<Double> subX = rescaleX.subList(start, i);
			//	System.out.println("start = "+ start+"i = "+i+ "sublistX = "+subX);
				
				
				Route singleRoute = new Route(rescaleX.subList(start, i),rescaleY.subList(start, i));
				rescaleRoutes.add(singleRoute);
				int j = i+1;
				while(j+1<rescaleX.size()&&(rescaleX.get(j)<-1000||rescaleX.get(j+1)<-1000))
					j++;
				
				start = j;
				i = j;
				//System.out.println("nextstart = "+ start);
				if(i<rescaleX.size()){
				trajId = rescaleX.get(i).intValue();
				posInTraj = 0;
				}
			}
			/*
			else{
				
				
				Integer[] map = {trajId,posInTraj};
				whole2separateTrajMap.add(map);
			}
			*/
		}
	}
	
	/*
	private void drawRawTrajectories() {
		System.out.println("drawRawTraj: "+rawAllIntervals.size());
		
		
	  			
	  		for (int k=0;k<rawAllIntervals.size();k++)
	  				  {
	  					  Route singleRoute = new Route();
	  					  int startPos = rawAllIntervals.get(k).getStartPos();
	  						int endPos = rawAllIntervals.get(k).getEndPos();
	  						
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
	*/

	private void runSequitur() {
			  chartData = new MotifChartData(this.dataFileName, lat, lon, 1, alphabetSize); //PAA is always 1.
		
			  filter = new ArrayList<Integer>();
		
			
			  rules = new GrammarRules();
				try{
				  SAXRecords saxFrequencyData = null;
				  saxFrequencyData = SequiturFactory.entries2SAXRecords(trimedTrack);
				  System.out.println("Input String Length: " + Tools.countSpaces(saxFrequencyData.getSAXString(SPACE)));
				  consoleLogger.trace("String: " + saxFrequencyData.getSAXString(SPACE));
				  System.out.println("String: "+ saxFrequencyData.getSAXString(SPACE));//.substring(0,1000));
				  consoleLogger.debug("running sequitur...");
				  
				  SAXRule sequiturGrammar = SequiturFactory.runSequitur(saxFrequencyData.getSAXString(SPACE));
				//  System.out.println("sequiturGrammar: "+sequiturGrammar.toGrammarRulesData().getRuleRecord(1));
				  consoleLogger.debug("collecting grammar rules data ...");
				 // GrammarRules rules1 = sequiturGrammar.toGRD();
				 // System.out.println("rules size: "+ rules1.size());			 
		          rules = sequiturGrammar.toGrammarRulesData();
		        //  rules.setParsedString();
		          realRuleSize = rules.size();
		     //     allRules.add(rules);
		          System.out.println("real rules size: "+ realRuleSize);
		          //debug
		          
		          
		          consoleLogger.debug("mapping rule intervals on timeseries ...");
		          
		         
		          
		          consoleLogger.debug("done ...");
		          
		       // SequiturFactory.updateRuleIntervals(filteredRules, saxFrequencyData, lat.size());
		    	  //        filteredRules = sequiturGrammar.toFilteredGrammarRulesData(filteredRuleMap);
		    	  //        for(int i=0;i<filteredRules.size();i++)
		    	  //      	  System.out.println(filteredRules.get(i));
		    	  SequiturFactory.updateRuleIntervals(rules,saxFrequencyData,rescaleX.size());
		          
		          
		          
		          chartData.setGrammarRules(rules);
		          
		          System.out.println("chartData size: "+ chartData.getRulesNumber());
		          System.out.println("rules size: "+ rules.size());
				  motifList = new ArrayList<Motif>();
				  motifListRuleMap = new ArrayList<Integer>();
		
			  }
			  catch (TSException e){
				  this.log("error while processing data "+StackTrace.toString(e));
				  e.printStackTrace();
			  }
			  
			  		
	}

	
	 
	  
		  /*    	20170228
		   * Generate All Motifs and record them on files respectively.
		   */
	/*
		private void plotOnMap() {
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
			    			mapOriginRules.add(set);
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
			    	
			    	for(int stepDist : clusters.get(i)){
			    	
			    		int rule = stepDist; //filter.get(stepDist);
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
			  	
				  {//(Tools.countSpaces(chartData.getRule(i).getExpandedRuleString())>minBlocks){
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
	
	/*
		private void finalCluster() {
			//currentClusters = new HashMap<String, Cluster>();
			for(int i = 0; i<r0.length; i++	){
				String s = r0[i];
				if(!Tools.isNumeric(s) && s!=null){
					if(r0[i].charAt(0)=='R'){
						s = "I"+iteration+"stepDist"+r0[i].substring(1);
						
						cluster = new Cluster(s);
						//cluster.addRule(rules.get(Integer.valueOf(r0[i].substring(1))));
						cluster.addRule(rules.get(Integer.valueOf(r0[i].substring(1))));
						//r0[i]=s;
					    currentClusters.put(s, cluster);	
					}

				}
			}
			RuleDistanceMatrix rdm = new RuleDistanceMatrix(blocks,currentClusters,minBlocks,maxPointErrorDistance);
			
		
	}
	*/

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
	        		  	int length2 = Tools.countSpaces(rules.get(currentRule).getExpandedRuleString());
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
			
	    	if(Tools.isNumeric(s)&&Integer.valueOf(s)>=0){
	    	
	    		int numStartPos;
	    		int numEndPos;
	    		if(c<1||Tools.isNumeric(r0[c-1])||r0[c-1]==null){
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
			
			if(!Tools.isNumeric(r0[i])){
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

	private void drawOnMap(){   //iterative sequitur
/*
		Comparator<Double> doubleComparator = new Comparator<Double>() {
	        @Override public int compare(Double s1, Double s2) {
	            return s1.compareTo(s2);
	        }           
	    };
	    		for(int i = 0; i<rules.size(); i++){
	    		  int startPos = rules.get(i).
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
		    		  if(Tools.isNumeric(r0[nextNonTerminal])) // negative
		    			  numEndPos = mapToOriginalTS.get(nextNonTerminal-nullCounter)-1;
		    		  else
		    			  numEndPos = mapToOriginalTS.get(nextNonTerminal-nullCounter);
		    		  
		    		  
		    		  currentAmountDistance = currentAmountDistance+allDiscordDistances.get(numStartPos+","+numEndPos);
		    		
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
		    	  totalRuleLength = totalRuleLength + Tools.countSpaces(RuleDistanceMatrix.parseRule(pair.getKey()));

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
		  					  	
		  			
		  				  else
		  					  i1++;
		  				
		  			  }
		  			int ruleCoverCount = 0;
		  			
		  			
		  			
		  			for (int a = 0;a<ruleCovered.length;a++){
		  				if(ruleCovered[a])
		  					ruleCoverCount++;
		  			}

		  			  drawAnomaly();
		  			  	  
		  			  
				*/
		  }

	private int getNextNonTerminal(int i) {
		int j = i+1;
	//	System.out.println("j: "+j);
		while(j<r0.length&&Tools.isNumeric(r0[j])&&Integer.valueOf(r0[j])>=0)
		{
			j++;
			
		}
		
		return j;
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
		  File evalFile = new File("./evaluation/"+"evaluate_"+(int)(maxPointErrorDistance*1000)+"_"+alphabetSize+"_"+minBlocks+"_"+resampleRate+"_"+lat.size()+"_"+fileNameOnly);
		  FileWriter fr = new FileWriter(evalFile);
		  String sb1; // = new StringBuffer();
		  //sb1.append(fileNameOnly+",");
		  sb1 = (fileNameOnly+","+maxPointErrorDistance+","+alphabetSize+","+minBlocks+","+resampleRate+","+runTime+","+avgIntraDistance+","+avgIntraDistanceStdDev+","+ minInterDistance+","+avgSilhouetteCoefficient+","+routes.size()+","+lat.size()+','+totalSubTrajectory+","+coverCount+","+immergableRuleCount+"\n");
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
		    	
		    	for(int stepDist : clusters.get(i)){
		    	
		    		int rule = stepDist; //filter.get(stepDist);
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
		  	
			  {//(Tools.countSpaces(chartData.getRule(i).getExpandedRuleString())>minBlocks){
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
					//if(distance[x][i]>(maxPointErrorDistance*2))
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
					//if(distance[i][y]>(maxPointErrorDistance*2))
					
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
			//	if(distance[i][j]>(maxPointErrorDistance*2))
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
			double avgDistance = Tools.avg(pairwiseDistances);
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
				 pairwiseInterDistances.add(Tools.avg(getSimilaritiesInterRules(allRules.get(i),allRules.get(j))));
			
			}
			if(pairwiseInterDistances.size()>0)
				allMinimalInterDistances.add(min(pairwiseInterDistances));

		}
		
//		System.out.println("average distances among all motifs: "+Tools.avg(allDistances));
	//	System.out.println("average standard deviation among all motifs: "+Tools.avg(allStdDev));
		//sb.append(Tools.avg(allDistances)+","+Tools.avg(allStdDev)+"\n");
		
		//return sb.toString();
		result[0] = Tools.avg(allDistances);  // avg intra distances
		result[1] = Tools.avg(allStdDev);  // avg std. dev. intra distances
		result[2] = Tools.avg(allMinimalInterDistances);  // avg minimal inter distances
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
		result[3] = Tools.avg(silhouetteCoefficients);
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
				pairwiseInterDistances.add(Tools.avgDTWDistance(eBlocks, rule1.get(i),rule2.get(j)));
			}
		
		return pairwiseInterDistances;
	}

	private double dev(ArrayList<Double> pairwiseDistances) {
		double avg = Tools.avg(pairwiseDistances);
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
				double similarity = Tools.avgDTWDistance(eBlocks, allTracks.get(i),allTracks.get(j));
				pairwiseDistance.add(similarity);
			}
		}
		return pairwiseDistance;
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
    //2017 new
    public static ArrayList<Route> getRescaleTrajectory(){
    	System.out.println("getrawTrajectory: "+rescaleRoutes.size());
    	return rescaleRoutes; 
    }
    public static ArrayList<Route> getAnomaly(){
    	return anomalyRoutes; 
    }
    public static ArrayList<Route> getDiscord(){
    	return discordRoutes; 
    }
    public static ArrayList<Route> getTrueAnomaly(){
    	return trueAnomalyRoutes; 
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
	 

	public static int getTrajectory(int endPos) {
		int i = endPos;
		Double traj;
		while(ncLat.get(i)>-999)
			i++;
		traj = -1000-ncLat.get(i);
		return traj.intValue();
	}
	

	public static Instance trajToInstance(Integer traj) {
		// TODO Auto-generated method stub
		return null;
	}

}
