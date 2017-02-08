package edu.gmu.trajviz.model;

import edu.gmu.trajviz.logic.Cluster;
import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.util.Tools;

public class RoutePair implements Comparable<RoutePair> {
	public String r1;
	public String r2;
	public int length;
	public double dist;
	public int[] t1;
	public int[] t2;
	private Route route1;
	private Route route2;
	public boolean isTrivial;
	//public double threshold;
	public RoutePair(String s1, String s2){
		r1 = s1;
		r2 = s2;
		t1 = Tools.parseTrajId(r1);
		t2 = Tools.parseTrajId(r2);
		//this.threshold = threshold;
		
		  
		  if(!Tools.isTrivialMatch(s1, s2)){
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
		  dist = Tools.routeSqrEuDist(route1, route2,SequiturModel.R);
		  }
		  else
			  {
			  	dist = Double.MAX_VALUE;
			  	isTrivial = true;
			  }
	}
	public RoutePair(int[] n1, int[] n2, double dist2) {
		t1 = n1;
		t2 = n2;
		dist = dist2;
		r1 = "T"+n1[0]+"S"+n1[1]+"L"+n1[2];
		r2 = "T"+n2[0]+"S"+n2[1]+"L"+n2[2];
		length = n1[2];
	}
	public RoutePair(int[] n1, int[] n2) {
		t1 = n1;
		t2 = n2;
		
		r1 = "T"+n1[0]+"S"+n1[1]+"L"+n1[2];
		r2 = "T"+n2[0]+"S"+n2[1]+"L"+n2[2];
		length = n1[2];
		if(!Tools.isTrivialMatch(r1, r2)){
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
			  dist = Tools.routeSqrEuDist(route1, route2,SequiturModel.R);
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
	/*
	public boolean isClose(double threshold) {
		
		for(int index = 0; index<length; index++){
			if (Tools.euDist(route1.getLats().get(index),route1.getLons().get(index),route2.getLats().get(index),route2.getLons().get(index))>threshold)
				return false;
			
		}	
		return true;
	}
	*/
	public boolean hasNextPair() {
		if((t1[1]+t1[2]+1)<SequiturModel.oldtrajX.get(t1[0]).size()&&(t2[1]+t2[2]+1)<SequiturModel.oldtrajX.get(t2[0]).size())
			return true;
		else
			return false;
	}
	public RoutePair nextPair() {
		if(this.hasNextPair()){
		int[] n1 = {t1[0],t1[1]+1,t1[2]};
		int[] n2 = {t2[0],t2[1]+1,t2[2]};
		/*
		double xhead = SequiturModel.oldtrajX.get(t1[0]).get(t1[1])-SequiturModel.oldtrajX.get(t2[0]).get(t2[1]);
		double yhead = SequiturModel.oldtrajY.get(t1[0]).get(t1[1])-SequiturModel.oldtrajY.get(t2[0]).get(t2[1]);
		double xtail = SequiturModel.oldtrajX.get(n1[0]).get(n1[1]+n1[2])-SequiturModel.oldtrajX.get(n2[0]).get(n2[1]+n2[2]);
		double ytail = SequiturModel.oldtrajY.get(n1[0]).get(n1[1]+n1[2])-SequiturModel.oldtrajY.get(n2[0]).get(n2[1]+n2[2]);
		double nextdist = dist-(xhead*xhead+yhead*yhead)+(xtail*xtail+ytail*ytail);
		*/
		RoutePair next = new RoutePair(n1,n2,0);
		return next;
		}
		else
			return null;
	}

}
