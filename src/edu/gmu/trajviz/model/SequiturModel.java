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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Observable;
import java.util.concurrent.TimeUnit;

import org.slf4j.LoggerFactory;

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
	public final static int EVAL_RESOLUTION = 500;
	final static Charset DEFAULT_CHARSET = StandardCharsets.UTF_8;
	public final static String EVALUATION_HEAD = "DataName,PaaSize,AlphabetSize,MinBlocks,NCThreshold,RunningTime,AvgDistance,AvgeStdDev,MinInterDistance,SilhouetteCoefficient,TotalRules,TotalDataPoints, TotalSubTrajectories,CoveredPoints, ImmergableRuleCount\n";
//	public static final int ALPHABETSIZE = 50;
//	public static final int CONTINUALBLOCKTHRESHOLD = 10;
	//public static final int paaSize = 10;
	private static final String SPACE = " ";
	private static final String CR = "\n";
	private static final int STEP = 2;
	
//	private static final int NOISYELIMINATIONTHRESHOLD = 5;
	private ArrayList<HashSet<Integer>> clusters;
	private ArrayList<Integer> filter;
	private String dataFileName, fileNameOnly;
	private static double lat_center;
	private double latMax;
	private double latMin;
	private double lonMin;
	private double lonMax;
	public int trajCounter;
	private static double lon_center;
	//The outer arrayList includes all rules, the inner arrayList includes all route under the same rule
	private static ArrayList<ArrayList<Route>> routes;  
	public ArrayList<Double> lat;
	public ArrayList<Double> paaLat;
	public ArrayList<Double> paaLon;
	public ArrayList<Double> lon;
	private MotifChartData chartData;
	private double runTime = -1;
	//private static GrammarRules filteredRules;
	public ArrayList<NumerosityReductionMapEntry> trimedTrack;
	// index is the rule# after filtering, Integer value is the actual rule number
//	private ArrayList<Integer> filteredRuleMap = new ArrayList<Integer>(); 
	public ArrayList<Integer> words;
	public Blocks blocks, eBlocks; 
	private static Logger consoleLogger;
	  private static Level LOGGING_LEVEL = Level.DEBUG;
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
		  // read the input
		  // init the data array
		  ArrayList<Double> data = new ArrayList<Double>();
		  ArrayList<Double> data1 = new ArrayList<Double>();
		  ArrayList<Integer> status = new ArrayList<Integer>();    // taxi loading status
		  ArrayList<Integer> timeAsUnixEpoc = new ArrayList<Integer>();   //number of seconds since Jan. 1 1970 midnight GMT, if the time is in milliseconds, it will need to be converted to seconds,otherwise it may over Integer's limite. 
		  
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
					  int value3 = Integer.parseInt(lineSplit[3]);
					  
					  if((lineCounter<=1)||(Math.abs(value3-timeAsUnixEpoc.get(timeAsUnixEpoc.size()-1))<=300))
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
			  this.log("error while trying to read data from " + this.dataFileName + ":\n" + stackTrace);
		  }
		  
		
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
			  this.lat.add(temp_latitude);
			  this.lon.add(temp_longitude);
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
		  System.out.println("lonMax:  "+lonMax+"       lonMin: "+lonMin);
		  System.out.println("latMax:  "+latMax+"       latMin: "+latMin);
		  System.out.println("Number of trajectories: "+trajCounter);
		  consoleLogger.debug("loaded " + this.lat.size() + " points and "+trajCounter+" Trajecoties... ");
		  this.log("loaded " + this.lat.size() + " points from " + this.dataFileName);
		  
		  
		  
		  setChanged();
		  notifyObservers(new SequiturMessage(SequiturMessage.TIME_SERIES_MESSAGE, this.lat,this.lon));
		  
	}
	  public static double getLatitudeCenter(){
		  return lat_center;
	  }
	  public static double getLongitudeCenter(){
		  return lon_center;
	  }
	  @SuppressWarnings("rawtypes")
	
	  
	  
	  
	  public synchronized void processData(double minLink, int alphabetSize, int minBlocks, int noiseThreshold)throws IOException{
		  StringBuffer sb = new StringBuffer();
		  if (null == this.lat ||null == this.lon|| this.lat.size()==0 || this.lon.size()==0 ){
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
		  long beginTime = System.currentTimeMillis();
			 // beginTime  =  System.nanoTime();
			  System.out.println("begin time: "+beginTime);
		  routes = new ArrayList< ArrayList<Route>>();
		  paaLat = new ArrayList<Double>();
		  paaLon = new ArrayList<Double>();
		  ArrayList<Double> latBuffer=new ArrayList<Double>();
		  ArrayList<Double> lonBuffer=new ArrayList<Double>();
		  double avgLat;
		  double avgLon;
		  /*
		   * use the centroid(paaLat,paaLon) to represent the data
		   */
		//  System.out.println("paaSize: "+paaSize);
		 /*
		  if(paaSize==1)
		  {
			  paaLat = lat;
			  paaLon = lon;
			  
		  }
		  else
		  */
		   // System.out.println("Should not see this msg.");
		  
		  for(int i=0;i<lat.size();i++){
			/*  if((i%paaSize)==0){
			  
			  }
			  */
			  latBuffer.add(lat.get(i));
			  lonBuffer.add(lon.get(i));
			  if((i+1==lat.size())||((i+1)%1)==0){
				  //compute the avg of the buffered arrayList into paaLat and paaLon
				  avgLat = avg(latBuffer);
				  avgLon = avg(lonBuffer);
				  for(int j=0; j<latBuffer.size();j++){
					  paaLat.add(avgLat);
					  paaLon.add(avgLon);
				  }
				  // refresh
				  latBuffer=new ArrayList<Double>();
				  lonBuffer=new ArrayList<Double>();
			  }
			  
		  }
		  
		/*  
		  System.out.println("oriLat" + lat);
		  System.out.println("paaLat:  " +paaLat);
		  System.out.println("oriLon" + lon);
		  System.out.println("paaLon:  " +paaLon);
		  */
		  blocks = new Blocks(alphabetSize,latMin,latMax,lonMin,lonMax);
		  
		  words = new ArrayList<Integer>();
		  // add all points into blocks.
		  Integer previousId=(Integer)(-1);
		
		//  HashMap<Integer,Integer> trackMap = new HashMap<Integer, Integer>();
		  trimedTrack = new ArrayList<NumerosityReductionMapEntry>();
		  
		  
		  for (int i = 0; i<paaLat.size();i++){
			  Location loc = new Location(paaLat.get(i),paaLon.get(i));
			//  blocks.addPoint2Block(loc); this should not work here because the point will change if it is a noisy point.
			  Integer id = new Integer(blocks.findBlockIdForPoint(loc));
			  if(isNoise(id,i,noiseThreshold)){
				 // lat.set(i, lat.get(i-1));
				  paaLat.set(i, paaLat.get(i-1));
				 // lon.set(i, lon.get(i-1));
				  paaLon.set(i, paaLon.get(i-1));
				  id = previousId;
			  }
			  
			  words.add(id);
		//	  System.out.println("previousId, id:  "+previousId+",   "+id+"        i:   "+i+"   Lat,Lon: "+paaLat.get(i)+","+paaLon.get(i));
			  
			  if (!id.equals(previousId))
			  {
				  NumerosityReductionMapEntry<Integer, Integer> entry = new NumerosityReductionMapEntry<Integer, Integer>(new Integer(i),id);
		//		  System.out.println("entry: "+i+","+id);
				  trimedTrack.add(entry);
				  //put the new <index,id> pair into a map 
				  //NumerosityReductionMapEntry entry = new NumerosityReductionMapEntry(i,id);
			//	  trackMap.put(i, id);
				  previousId = id;
			  }
			  
			  				  
		  }
		  
		  /*
		   * Following is put the cleaned location data into block again
		   */
		  /*
		  for(int i=0; i<20;i++)
		  System.out.println("Orignal String: " + words.get(i));
		  
		  System.out.println("StringTrimedTrack:  "+trimedTrack.size());
		  for(int i=0; i<trimedTrack.size();i++){
			  System.out.println(i+" : "+trimedTrack.get(i).getValue()+" ");
		  }
		  System.out.println();
		  
		  */
		  this.chartData = new MotifChartData(this.dataFileName, paaLat, paaLon, 1, alphabetSize); //PAA is always 1.
	//	  GrammarRules filteredRules = new GrammarRules();
		  
	//	   filteredRuleMap = new ArrayList<Integer>(); // index is the rule# after filtering, Integer value is the actual rule number. 
		  clusters = new ArrayList<HashSet<Integer>>();
		  filter = new ArrayList<Integer>();
		  HashMap<Integer, Integer> clusterMap = new HashMap<Integer,Integer>();
		   long buildMatrixTime = 0;
		   long clusterTime = 0;
		  try{
			  SAXRecords saxFrequencyData = null;
			  saxFrequencyData = SequiturFactory.entries2SAXRecords(trimedTrack);
			//  System.out.println("String: " + saxFrequencyData.getSAXString(SPACE));
			  consoleLogger.trace("String: " + saxFrequencyData.getSAXString(SPACE));
			//  System.out.println("String: "+ saxFrequencyData.getSAXString(SPACE));
			  consoleLogger.debug("running sequitur...");
			  
			  SAXRule sequiturGrammar = SequiturFactory.runSequitur(saxFrequencyData.getSAXString(SPACE));
			//  System.out.println("sequiturGrammar: "+sequiturGrammar.toGrammarRulesData().getRuleRecord(1));
			  consoleLogger.debug("collecting grammar rules data ...");
			 // GrammarRules rules1 = sequiturGrammar.toGRD();
			 // System.out.println("rules size: "+ rules1.size());			 
	          GrammarRules rules = sequiturGrammar.toGrammarRulesData();
	          System.out.println("rules size: "+ rules.size());
	          //debug
	          
	          consoleLogger.debug("mapping rule intervals on timeseries ...");
	          SequiturFactory.updateRuleIntervals(rules, saxFrequencyData, lat.size());
	          HashMap<String, Integer> hm = new HashMap<String, Integer>();
	          String rule0 = rules.get(0).getRuleString();
	          String[] r0 = rule0.split(" ");
	          for(int i = 0; i<rules.size();i++){
	        	  String key = rules.get(i).getRuleName();
	        	//  System.out.println(rules.get(i));
	        	  hm.put(key, 0);
	          }
	          /*
	          for(int i = 0; i<r0.length;i++){
	        	          
	          System.out.print(r0[i]+" ");
	          }
	          */
	          System.out.println();
	          for(int i=0;i<r0.length;i++){
	        	  if(r0[i].contains("R"))
	        		  hm.put(r0[i], hm.get(r0[i])+1);
	          }
	          for(int i = 1; i<rules.size();i++){
	        	  
	        	  String key = rules.get(i).getRuleName();
	        	  rules.get(i).setFrequencyInR0((hm.get(key)).intValue());
	          }
	      /*    
	          for(int i=0;i<rules.size();i++){
	        	  System.out.println("Rule number: "+rules.getRuleRecord(i).getRuleNumber()+" Fre in R0: "+rules.get(i).frequencyInR0()+" "+rules.get(i)+" StringOccurance: "+rules.getRuleRecord(i).occurrencesToString()+" Rule String: "+rules.getRuleRecord(i).getExpandedRuleString()+" Rule Positions: "+rules.getRuleRecord(i).getRuleIntervals());
	          }
	         
	          */
	          
	          /*
	           * Postprocessing merge, connect
	           */
	          for (int i = 0; i<rules.size();i++){
					if ((rules.get(i).frequencyInR0()>2&&rules.get(i).getRuleYield()>=minBlocks)||
							(rules.get(i).frequencyInR0()>1&&rules.get(i).getRuleIntervals().size()>2&&rules.get(i).getRuleYield()>=minBlocks))
						{
						//HashSet<Integer> set = new HashSet<Integer>();
		//				System.out.println("Yield: "+rules.get(i).getRuleYield()+" string: "+rules.get(i).getExpandedRuleString());
						filter.add(i);
						if(rules.get(i).getRuleIntervals().size()<=2)
							System.out.println("Bug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"+i);
						}
					
				}
	          System.out.println("filter Size = "+filter.size());
	          
	          //HashMap<Integer,ArrayList<Integer>> mergeRecord = new HashMap<Integer, ArrayList<Integer>>();
	          long t1s = System.currentTimeMillis();
	          RuleDistanceMatrix rdm = new RuleDistanceMatrix(blocks,rules, filter,minBlocks, minLink); 
	          long t1e = System.currentTimeMillis();
	          buildMatrixTime = t1e-t1s;
	          
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
	        		  mergableCount++;
	        	//	  merge(rules,rdm.filter.get(pair.getLine()),rdm.filter.get(pair.getCol()));
	        		  if(clusterMap.containsKey(pair.getLine())||clusterMap.containsKey(pair.getCol()))
	        		  {
	        			  if(!clusterMap.containsKey(pair.getLine())){
	        				  clusters.get(clusterMap.get(pair.getCol())).add(pair.getLine());
	        				  clusterMap.put(pair.getLine(), clusterMap.get(pair.getCol()));
	        			//	  System.out.println("Adding Line  to a cluster, Line:"+pair.getLine()+" Colu:"+pair.getCol()+clusters.get(clusterMap.get(pair.getCol())));
	        				  //System.out.println("Map:"+clusterMap);
	        				  
	        			  }
	        			  else if(!clusterMap.containsKey(pair.getCol())){
	        				  clusters.get(clusterMap.get(pair.getLine())).add(pair.getCol());
	        				  clusterMap.put(pair.getCol(), clusterMap.get(pair.getLine()));
	        			//	  System.out.println("Adding Colum to a cluster,Colum:"+pair.getCol()+" Colu:"+pair.getCol()+clusters.get(clusterMap.get(pair.getLine())));
	        				  //System.out.println("Map:"+clusterMap);
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
	        					  clusterMap.put(v, clusterMap.get(pair.getLine()));
	        				//	  clusters.get(clusterMap.get(pair.getLine())).add(v);
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
	        			  set.add(pair.getLine());
	        			  set.add(pair.getCol());
	        			  clusters.add(set);
	        			  clusterMap.put(pair.getLine(), clusters.size()-1);
	        			  clusterMap.put(pair.getCol(), clusters.size()-1);
	        			 // System.out.println("Created a cluster: "+clusters.get(clusters.size()-1));
	        		//	  System.out.println("Map:"+clusterMap);
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
	          
	          long t2e = System.currentTimeMillis();
	          clusterTime = t2e -t2s;

	          /*
	          for(int i = 0; i<clusters.size();i++){
	        	  System.out.println("i = "+i+" : "+clusters.get(i));
	          }
	          
	          */
	          
	          
	       //   double[][] ruleDistances
	      //   GrammarRules preprocessedRules = new GrammarRules();
	         /*
	         for (int i = 0; i<rules.size();i++)
	        	 {
	        	 	if(rules.get(i).frequencyInR0()>2)
	        	 		preprocessedRules.addRule(rules.get(i));
	        	 }
	        	 */
	         
	     //    while
	     /*
	         while(rdm.getMinDistance()<0.1){
	        	 rules.merge(rdm.getMinPair()[0],rdm.getMinPair()[1]);
	        	 
	        	 rdm = new RuleDistanceMatrix(blocks, rules);
	         }
	       */  
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          
	          long endTime = System.currentTimeMillis();
	          System.out.println("end time: "+endTime);
	          runTime = (endTime-beginTime)/1000.0;
	          
	          consoleLogger.debug("done ...");
	          
	          
	          
	          
	          /*
	          
	          //filter the rules
	          for(int i=0; i<rules.size();i++){
	        	  if((countSpaces(rules.getRuleRecord(i).getExpandedRuleString())>minBlocks)&&(rules.getRuleRecord(i).getRuleIntervals().size()>2)){
	        	//	  filteredRules.addRule(rules.getRuleRecord(i));
	        		  filteredRuleMap.add(i);
	        	  }
	          }
	          System.out.println("filterrulemap:   "+filteredRuleMap);
	          */
	           // SequiturFactory.updateRuleIntervals(filteredRules, saxFrequencyData, lat.size());
	  //        filteredRules = sequiturGrammar.toFilteredGrammarRulesData(filteredRuleMap);
	  //        for(int i=0;i<filteredRules.size();i++)
	  //      	  System.out.println(filteredRules.get(i));
	  //        SequiturFactory.updateRuleIntervals(filteredRules,saxFrequencyData,lat.size());
	            this.chartData.setGrammarRules(rules);
	   //       this.chartData.setGrammarRules(filteredRules);
	          System.out.println("chartData size: "+ chartData.getRulesNumber());
			  
			  
		  }
		  catch (TSException e){
			  this.log("error while processing data "+StackTrace.toString(e));
			  e.printStackTrace();
		  }
		  this.log("processed data, painting on map");
		  consoleLogger.info("process finished");
		  setChanged();
		  /*
		   * Generate All Motifs and record them on files respectively.
		   */
		 // String header = "type,latitude,longitude";
//		  System.out.println("Total rules:"+chartData.getRulesNumber());
		  
		//  ArrayList<SAXMotif> allMotifs = chartData.getAllMotifs();
		//  for (int i=1; i<chartData.getRulesNumber();i++){
		    // create merged rule interval data structure corresponding to "clusters" 
		    ArrayList<ArrayList<RuleInterval>> ruleIntervals = new ArrayList<ArrayList<RuleInterval>>();
	          ArrayList<HashSet<Integer>> mapToOriginRules = new ArrayList<HashSet<Integer>>();
	        System.out.println("cluster map size = "+ clusterMap.size());
		    System.out.println("clusterMap:   "+clusterMap);
		    int totalRuleCount = 0;
		    int immergableRuleCount = 0;
		    for(int i=0;i<filter.size();i++){
		    	if(!clusterMap.containsKey(i)&&chartData.getRulePositionsByRuleNum(filter.get(i)).size()>=minBlocks)
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
		    	
		    		int rule = filter.get(r);
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
		    int coverCount = 0;
		    for (int i = 0; i<isCovered.length;i++)
		    	isCovered[i] = false;
		    int totalSubTrajectory = 0;
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
					/*
					 * Generating evaluation file
					 */
					  
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
						/*
						if(k==0){
							route0Id = new Integer(blocks.findBlockIdForPoint(new Location(paaLat.get(startPos),paaLon.get(startPos))));
						}
						if(k==1){
							route1Id = new Integer(blocks.findBlockIdForPoint(new Location(paaLat.get(startPos),paaLon.get(startPos))));
						}
						*/
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
				  /* bug fixed!!!
				   if(!route0Id.equals(route1Id))
				   {
					   //route.remove(0);
					   System.out.println("i: "+i);
					  	System.out.println(chartData.getRule(filteredRuleMap.get(i)).getExpandedRuleString());
					  	System.out.println(chartData.getRulePositionsByRuleNum(filteredRuleMap.get(i)));
					//  	System.out.println(singleRoute.);
				   }
				    */	
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
		  System.out.println("cover rate: " +coverCount/isCovered.length);
		  /*
		   * Generate All Motifs and record them on files respectively.
		   */
		  
		  
		  
		  
		 
		  
		  //test
		/*  blocks.printBlockMap();
		  for (int i = 0; i<words.size(); i++)
		  {
			  System.out.print("  "+words.get(i));
		  }
		  */
		  /*r
		  System.out.println("trackMap:");
		  System.out.println(trackMap.toString());
		 // System.out.println(map2String(trackMap));
		  System.out.println("Postions:\t"+getTrimedPositions(trimedTrack).toString()+"\t");
		  System.out.println("TrimedStrs:\t"+getTrimedIds(trimedTrack));
		  */
			System.out.println("build matrix: "+(double)(buildMatrixTime/1000.0));
		  System.out.println("Clustering Time: "+(clusterTime/1000.0));

		  System.out.println("running time: "+runTime);
		  ArrayList<Integer> frequency = new ArrayList<Integer>();
		  /*
		  for (int i=0;i<ruleIntervals.size();i++)
			  frequency.add(ruleIntervals.get(i).size());
			  */
		  notifyObservers(new SequiturMessage(SequiturMessage.CHART_MESSAGE, this.chartData, ruleIntervals, mapToOriginRules));//, frequency ));
	  
		  
		  
		  /*
		   * evaluation
		   */
	
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
	  
	

	private boolean isMergable(double[][] distance, ArrayList<HashSet<Integer>> families, int x, int y, HashMap<Integer, Integer> map, double minLink) {
		//boolean mergable = true;
		if(map.containsKey(x)||map.containsKey(y)){
			if(!map.containsKey(x)){
				for(int i: families.get(map.get(y)))
				
					if(distance[x][i]>(minLink*2))
						return false;
			}
			else if(!map.containsKey(y)){
				for(int i: families.get( map.get(x)))
					if(distance[i][y]>(minLink*2))
						return false;
			}
			else
			{	
			for (int i : families.get(map.get(x)))
				for(int j : families.get(map.get(y)))

				{
		//		int xSibling = families.get(x).get(i);
		//		int ySibling = families.get(y).get(j);
				if(distance[i][j]>(minLink*2))
					return false;
				}
			}
		}
		else if(distance[x][y]>(minLink*2))
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
			
			if(allDistances.get(i)<allMinimalInterDistances.get(i)){
				sc = 1 - allDistances.get(i)/allMinimalInterDistances.get(i);
			}
			else{
				if(allDistances.get(i)==0)
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
		if(list == null)
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
		  {	Location loc = new Location(paaLat.get(i+j),paaLon.get(i+j));
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
	  private static int countSpaces(String str) {
	    int counter = 0;
	    for (int i = 0; i < str.length(); i++) {
	      if (str.charAt(i) == ' ') {
	        counter++;
	      }
	    }
	    return counter;
	  }

}
