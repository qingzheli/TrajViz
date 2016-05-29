package edu.gmu.trajviz.model;

import edu.gmu.trajviz.logic.Cluster;
import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.util.Tools;

public class RoutePair implements Comparable<RoutePair> {
	public String r1;
	public String r2;
	public int length;
	public double dist;
	private int[] t1;
	private int[] t2;
	public Route route1;
	public Route route2;
	public boolean isTrivial;
	public double threshold;
	public RoutePair(String s1, String s2,double threshold){
		r1 = s1;
		r2 = s2;
		t1 = Tools.parseTrajId(r1);
		t2 = Tools.parseTrajId(r2);
		this.threshold = threshold;
		
		  
		  if(t1[0]!=t2[0]){
		  if (t1[2]!=t2[2] )
			throw new IllegalArgumentException(r1+" is not comparable with "+ r2+"          t1[2] = "+t1[2]+"   t2[2] = "+ t2[2]);
		  length = t1[2];
		  route1 = new Route();
		  route2 = new Route();
		  for(int i = 0; i<length; i++){
			  double lat1 = SequiturModel.oldtrajX.get(t1[0]).get(t1[1]+i);
			  double lon1 = SequiturModel.oldtrajY.get(t1[0]).get(t1[1]+i);
			  route1.addLocation(lat1, lon1);
			  double lat2 = SequiturModel.oldtrajX.get(t2[0]).get(t2[1]+i);
			  double lon2 = SequiturModel.oldtrajY.get(t2[0]).get(t2[1]+i);
			  route2.addLocation(lat2, lon2);
		  }
		 // dist = Tools.closeRouteEuDist(route1, route2,threshold);
		  dist = Tools.routeEuDist(route1, route2);
		  }
		  else
			  {
			  	dist = Double.MAX_VALUE;
			  	isTrivial = true;
			  }
	}
	/*
	public RoutePair(String s1, String s2){
		r1 = s1;
		r2 = s2;
		t1 = Tools.parseTrajId(r1);
		t2 = Tools.parseTrajId(r2);
		
		
		  
		  if(t1[0]!=t2[0]){
		  if (t1[2]!=t2[2] )
			throw new IllegalArgumentException(r1+" is not comparable with "+ r2+"          t1[2] = "+t1[2]+"   t2[2] = "+ t2[2]);
		  length = t1[2];
		  route1 = new Route();
		  route2 = new Route();
		  for(int i = 0; i<length; i++){
			  double lat1 = SequiturModel.oldtrajX.get(t1[0]).get(t1[1]+i);
			  double lon1 = SequiturModel.oldtrajY.get(t1[0]).get(t1[1]+i);
			  route1.addLocation(lat1, lon1);
			  double lat2 = SequiturModel.oldtrajX.get(t2[0]).get(t2[1]+i);
			  double lon2 = SequiturModel.oldtrajY.get(t2[0]).get(t2[1]+i);
			  route2.addLocation(lat2, lon2);
		  }
		  dist = Tools.routeEuDist(route1, route2);
		  }
		  else
			  {
			  	dist = Double.MAX_VALUE;
			  	isTrivial = true;
			  }
	}
	*/
	public String toString(){
		return "distance ( "+r1+", "+r2+" ) = "+dist+"\n";
	}
	public boolean isTrivial(){
		if(t1[0]==t2[0])
			return true;
		else
			return false;
	}
	
	@Override
	public int compareTo(RoutePair pair2) {
		if (this.length!=pair2.length )
			throw new IllegalArgumentException(this+" is not comparable with "+ pair2);
		 
		 
		  if(this.dist<pair2.dist)
			  return -1;
		  if(this.dist>pair2.dist)
			  return 1;
		  	return 0;
	}
	public boolean isClose(double threshold) {
		
		for(int index = 0; index<length; index++){
			if (Tools.euDist(route1.getLats().get(index),route1.getLons().get(index),route2.getLats().get(index),route2.getLons().get(index))>threshold)
				return false;
			
		}	
		return true;
	}

}
