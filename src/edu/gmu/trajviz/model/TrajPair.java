package edu.gmu.trajviz.model;

import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.sax.datastructures.Motif;
import edu.gmu.trajviz.util.Tools;

public class TrajPair implements Comparable<TrajPair> {
	public Integer r1;
	public Integer r2;
	public int length;
	public double dist;
	
	public Route route1;
	public Route route2;
	public boolean isTrivial;
	public TrajPair(int traj1, int traj2){
		r1 = traj1;
		r2 = traj2;
		
		  
		  if(r1!=r2){
		  
		  route1 = new Route();
		  route2 = new Route();
		  for(int i = 0; i<SequiturModel.oldtrajX.get(r1).size(); i++){
			  double lat1 = SequiturModel.oldtrajX.get(r1).get(i);
			  double lon1 = SequiturModel.oldtrajY.get(r1).get(i);
			  route1.addLocation(lat1, lon1);
		  }
		  for(int i = 0; i<SequiturModel.oldtrajX.get(r2).size(); i++){
			  double lat2 = SequiturModel.oldtrajX.get(r2).get(i);
			  double lon2 = SequiturModel.oldtrajY.get(r2).get(i);
			  route2.addLocation(lat2, lon2);
		  }
		  
			
		  dist = Tools.RouteSITED(route1, route2,SequiturModel.R);
		  }
		  else
			  {
			  	dist = Double.MAX_VALUE;
			  	isTrivial = true;
			  }
	}
	public String toString(){
		return "distance ( "+r1+", "+r2+" ) = "+dist+"\n";
	}
	public boolean isTrivial(){
		if(r2==r1)
			return true;
		else
			return false;
	}
	
	@Override
	public int compareTo(TrajPair pair2) {
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

}
