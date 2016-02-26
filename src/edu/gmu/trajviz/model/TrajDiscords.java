package edu.gmu.trajviz.model;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Scanner;


public class TrajDiscords{
  public static String fileName = "1.txt";
  public static final int MIN = 2;
  public static final int MAX = 16;
  public static final int TOP_K = 200;
  public static ArrayList<Double> lats,lons;
  public static final long BASE_TIME = 1256171117;
  public HashMap<Integer,Double> minDistances;
 // public static PriorityQueue<DiscordRecord> pq;
  
  public static HashMap<Integer,SortedDiscodRecordArray> allDiscords;
 // public static ArrayList<PriorityQueue<Integer>> allDiscords;
  public static HashMap<Integer, SortedDiscodRecordArray> getAllDiscords(){
	
	 lats = SequiturModel.ncLat;
	 lons = SequiturModel.ncLon;
	 allDiscords = new HashMap<Integer,SortedDiscodRecordArray>();
	 for (int i = MIN; i<=MAX; i++)  //loop for all length discords
		 {
	//	 	pq = new PriorityQueue<DiscordRecord>();
		    
		// 	topDiscords = new SortedArray(TOP_K);
		 	allDiscords.put(i, findDiscords(i,10));
		 	
		 }
	
	 for (int i = MIN; i<=MAX; i++)  //loop for all length discords
	 {

	 	System.out.println(allDiscords.get(i));
	 	
	 }
	return allDiscords;
  }
  /*
   * rank is the rank th nearest distance for inner loop
   */
  private static SortedDiscodRecordArray findDiscords(int size, int rank) {
	  if(rank<1)
		  throw new IllegalArgumentException("rank = "+rank);
	SortedDiscodRecordArray topDiscords = new SortedDiscodRecordArray(TOP_K);
	
  
	for(int i = 0; i<lats.size()-size; i++ )    //outer loop
	{   
		double[] minDistances = new double[rank+1];
		for (int index = 0;index<=rank; index++){
		minDistances[index] = Double.MAX_VALUE;
	}   
		int start = i;
		int end = i+size-1;
		
		if(isCross(start,end))
			continue;
		double currentMinDist = Double.MAX_VALUE;
		double secondMinDist = Double.MAX_VALUE;
		for (int j=0; j<lats.size()-size; j++) //inner loop
		{
			
			int s = j;
			int e = j+size-1;
			if((s>start-size+1&&s<end)||(isCross(s,e)))
				continue;
			else
			{
				
				double dist = pairEuDist(lats,lons,start,end,s,e);
				minDistances[rank] = dist;
				Arrays.sort(minDistances);
			//	System.out.println("dist = "+dist);
				
			/*
				if(dist<currentMinDist)
					{
					secondMinDist = currentMinDist;
					
					
					currentMinDist = dist;
					//System.out.println("currentMinDist = "+currentMinDist);
					}
					
					*/
			}
			
		}
		DiscordRecord discord;
		
			discord = new DiscordRecord(start,end,minDistances[rank-1]);
			SequiturModel.allDiscordDistances.put(start+","+end, minDistances[0]);
		//else
			//discord = new DiscordRecord(start,end,secondMinDist);
	//	System.out.println("current Discord: "+discord);
		topDiscords.add(discord);
	}
	return topDiscords;
}
private static boolean isCross(int start, int end) {
for (int x = start; x<=end; x++)
	if (x%17==16)
		return true;
return false;
}

private static double pairEuDist(ArrayList<Double> array1, ArrayList<Double> array2, int s, int e, int s1,int e1) {
	if(array1.size()!=array2.size()||e-s!=e1-s1)
		throw new ArrayIndexOutOfBoundsException("ts1.length = "+array1.size()+"   ts2.length = "+array2.size());
	double eDist = 0;
	
	for(int i = 0; i<=e-s; i++){
		//if(e1>4415)
		//System.out.println("e:s:i = "+(e+i)+":"+s+":"+i+":"+(e1+i)+":"+s1+":"+array1.size()+":"+array2.size());
	
		eDist = eDist+euDist(array1.get(s+i),
			  array2.get(s+i),
			  array1.get(s1+i),
			  array2.get(s1+i));
			 
	/*
		eDist = eDist+distFrom(array1.get(s+i),
				  array2.get(s+i),
				  array1.get(s1+i),
				  array2.get(s1+i));
				  */
	}
	return eDist;
}

/*
 * Euclidean Distance
 */
private static double euDist(double x1, double y1 , double x2, double y2) {
	if(x1<-999||x2<-999)
		return Double.MAX_VALUE;
	return Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

public static double distFrom(double lat1, double lon1, double lat2, double lon2) {
	    double earthRadius = 6371.0; //3958.75 miles or 6371.0 kilometers
	    double dLat = Math.toRadians(lat2-lat1);
	    double dLng = Math.toRadians(lon2-lon1);
	    double sindLat = Math.sin(dLat / 2);
	    double sindLng = Math.sin(dLng / 2);
	    double a = Math.pow(sindLat, 2) + Math.pow(sindLng, 2)
	            * Math.cos(Math.toRadians(lat1)) * Math.cos(Math.toRadians(lat2));
	    double c = 2 * Math.atan2(Math.sqrt(a), Math.sqrt(1-a));
	    double dist = earthRadius * c;

	    return dist;
	    }


public static int getTrajectory(int endPos) {
	int i = endPos;
	Double traj;
	while(lats.get(i)>-999)
		i++;
	traj = -1000-lats.get(i);
	return traj.intValue();
}
}