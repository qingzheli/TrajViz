package edu.gmu.trajviz.sax.datastructures;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.model.SequiturModel;
import edu.gmu.trajviz.util.Tools;
//-qz
public class Motif {
	public int id; //motif id should be same to rule # before merge. 
	public int length;  // # of points
	public Route repRoute;
	public ArrayList<Double> repLineX,repLineY;
	private ArrayList<Route> routes; //<StartPositionInRescaledX, Route>
	public SortedSet<Integer> startPositions; 
	public SortedSet<Integer> trajIds;
	//private int firstTraj;
	//private int firstStartPos;
	//private boolean isSecond;
	//public HashSet<String> trajIds;
	
	public Motif(int ruleNumber,Route firstRoute, int startPosition){
		id = ruleNumber;
		length = firstRoute.size();
		routes = new ArrayList<Route>();
		routes.add(firstRoute);
		startPositions = new TreeSet<Integer>();
		startPositions.add(startPosition);
		trajIds = new TreeSet<Integer>();
		trajIds.add(SequiturModel.whole2separateTrajMap.get(startPosition)[0]);
		
		
	}
	/*
	public Cluster(int length,String s1,String s2){
		this.length = length;
		trajIds = new HashSet<String>();
		trajIds.add(s1);
		trajIds.add(s2);
		routes = new ArrayList<Route>();
		firstTraj = -1;
		firstStartPos = -1;
		isSecond = true;
	}
*/
	public int size(){
		return routes.size();
	}
	
	/*
	 * if repLine==null, repLine  = comming subtrajectory;
	 */
	public boolean add(Integer startPosition, Route route){
		
		
		if((!trajIds.contains(SequiturModel.whole2separateTrajMap.get(startPosition)[0]))||isNonTrivial(startPosition)){
			trajIds.add(SequiturModel.whole2separateTrajMap.get(startPosition)[0]);
			startPositions.add(startPosition);
			routes.add(route);
		}
		return true;
	
	}
	private boolean isNonTrivial(Integer startPosition) {
		
		SortedSet<Integer>candidatePositions = startPositions.subSet(startPosition-length, startPosition+length);
		if(candidatePositions.size()>0)
			return false;
		/*
		Iterator it = startPositions.iterator();
		while(it.hasNext()){
			int next = (int) it.next();
			if(Math.abs(startPosition-next)<length)
				return false;
		}
		*/
		return true;
	}
/*
	public boolean addBest(int traj, int s){
		String name = "T"+traj+"S"+s+"L"+length;
		if(repLineX==null || repLineY ==null){
			repLineX = new ArrayList<Double>();
			repLineY = new ArrayList<Double>();
			firstTraj = traj;
			firstStartPos = s;
			Route route = new Route();
			for(int i = s; i<s+length; i++){
				Double x = SequiturModel.oldtrajX.get(traj).get(i);
				Double y = SequiturModel.oldtrajY.get(traj).get(i);
			repLineX.add(x); 
			repLineY.add(y);
			
			route.addLocation(x,y);
			}
			trajIds.add(name);
			routes.put(name, route);
			repRoute = new Route(repLineX,repLineY);
			return true;
		}
		else if(nonTrivial(traj)){
			if(isSecond){
				for(int k = firstStartPos; k<firstStartPos+length; k++){
					
				
				SequiturModel.isAnomaly.get(firstTraj).set(k, false);
				}
				isSecond = false;
			}
			Route route = new Route();
			for(int i = s; i<s+length; i++){
				Double x = SequiturModel.oldtrajX.get(traj).get(i);
				Double y = SequiturModel.oldtrajY.get(traj).get(i);
				
				SequiturModel.isAnomaly.get(traj).set(i, false);
				repLineX.set(i-s, (repLineX.get(i-s)*trajIds.size()+x)/(trajIds.size()+1));
				repLineY.set(i-s, (repLineY.get(i-s)*trajIds.size()+y)/(trajIds.size()+1));
				route.addLocation(x,y);
			}
			trajIds.add(name);
			routes.put(name, route);
			repRoute = new Route(repLineX,repLineY);
		//	System.out.println("addToCluster: "+trajIds);
			return true;
		}
		else
			return false;
	}
	*/


public ArrayList<Route> getRoutes(){
	return routes;
}
private boolean isClose(int traj, int s) {
		
		for(int index = 0; index<length; index++){
			
		//	System.out.println("Eudist: "+SequiturModel.euDist(repLineX.get(index),repLineY.get(index),SequiturModel.oldtrajX.get(traj).get(s+index),SequiturModel.oldtrajY.get(traj).get(s+index)));
		//	if (Tools.euDist(repLineX.get(index),repLineY.get(index),SequiturModel.oldtrajX.get(traj).get(s+index),SequiturModel.oldtrajY.get(traj).get(s+index))>SequiturModel.maxPointErrorDistance)
			if (Tools.pointEuDist(repLineX.get(index),repLineY.get(index),SequiturModel.oldtrajX.get(traj).get(s+index),SequiturModel.oldtrajY.get(traj).get(s+index))>SequiturModel.distCut*length*SequiturModel.maxPointErrorDistance)

			return false;
			
		}	
		return true;
	}

public String toString(){
	return id+": size = "+routes.size();
}

/*
 * distance based
 * if a trajectory in A in cluster A is close to cluster B's repTraj, get a into b.
 */
/*
public void merge(Cluster c2,double threshold,HashMap<String, Cluster> findCluster) {
	if(this.getSize()>=c2.getSize()){
	Iterator it = c2.trajIds.iterator();
	ArrayList<String> removeCandidate = new ArrayList<String>();
	while (it.hasNext()){
		String name = (String) it.next();
		if(nonTrivial(Tools.parseTrajId(name)[0])){
	    Route route = c2.routes.get(name);
	  //  if(Tools.closeRouteEuDist(route, this.repRoute, R)/(this.getSize()+1)<=Tools.closeRouteEuDist(route, c2.respRoute, R)/(c2.getSize()+1)){
	 //   if(Tools.closeRouteEuDist(route, this.repRoute, R)<=Tools.closeRouteEuDist(route, c2.repRoute, R)){	
	    //if(Tools.closeRouteEuDist(route, this.repRoute, threshold)<=threshold)
	    if(Tools.routeSqrEuDist(route, this.repRoute,SequiturModel.R)<=threshold)
	    {	
  
	    	int[] subTraj = Tools.parseTrajId(name);
	    	this.add(subTraj[0], subTraj[1]);
	    	removeCandidate.add(name.toString());
	    //	System.out.println("cluster merged!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	    }
	    else
	    {
	    //	System.out.println("failed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
	    }
		}
	}
	if(c2.getSize()==1){
		String name = (String) c2.trajIds.toArray()[0];
		int[] subTraj = Tools.parseTrajId(name);
		this.add(subTraj[0], subTraj[1]);
		removeCandidate.add(name.toString());
	}
	for(int i = 0; i<removeCandidate.size(); i++){
		findCluster.put(removeCandidate.get(i), this);
		c2.remove(removeCandidate.get(i));
	}
	}
	else{
		Iterator it = this.trajIds.iterator();
		ArrayList<String> removeCandidate = new ArrayList<String>();

		while (it.hasNext()){
			String name = (String) it.next();
		    Route route = this.routes.get(name);
		  //  if(Tools.closeRouteEuDist(route, c2.repRoute, R)/(c2.getSize()+1)<=Tools.closeRouteEuDist(route, this.repRoute, R)/(this.getSize()+1)){
		 //   if(Tools.closeRouteEuDist(route, c2.repRoute, R)<=Tools.closeRouteEuDist(route, this.repRoute, R)){
		//    if(Tools.closeRouteEuDist(route, c2.repRoute, threshold)<=threshold){
		    if(Tools.routeSqrEuDist(route, this.repRoute,SequiturModel.R)<=threshold){
		    	int[] subTraj = Tools.parseTrajId(name);
		    	c2.add(subTraj[0], subTraj[1]);
		    	removeCandidate.add(name.toString());
		    //	System.out.println("cluster merged!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		    }
		    else
		    {
		    //	System.out.println("failed!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
		    }
		}
		if(this.getSize()==1){
			String name = (String) this.trajIds.toArray()[0];
			int[] subTraj = Tools.parseTrajId(name);
			c2.add(subTraj[0], subTraj[1]);
			removeCandidate.add(name.toString());
		}
		for(int i = 0; i<removeCandidate.size(); i++){
			findCluster.put(removeCandidate.get(i), c2);
			this.remove(removeCandidate.get(i));
		}
	
	}
	
	
	
}

/*
private void remove(String name) {
	
	Route route = this.routes.get(name);
	if(trajIds.size()>1)
		for(int i = 0; i<length; i++){
			repLineX.set(i, (repLineX.get(i)*trajIds.size()-route.getLats().get(i))/trajIds.size()-1);
			repLineY.set(i, (repLineY.get(i)*trajIds.size()-route.getLons().get(i))/trajIds.size()-1);
		}
	this.trajIds.remove(name);
	this.routes.remove(name);
}
@Override
public int compareTo(Object o) {
Cluster c = (Cluster) o;
if(this.getSize()>c.getSize())
	return 1;
else if(this.getSize()<c.getSize())
	return -1;
else	
	return 0;
}
public boolean findTraj(int[] i) {
	Iterator it = trajIds.iterator();
	while(it.hasNext()){
		String name = (String) it.next();
		int[] result =Tools.parseTrajId(name);
		if((i[0] == result[0])&&(Math.abs(i[1]-result[1])<length))
			{
				System.out.println("my best result: = "+name);
				return true;
			}
	}
	return false;
}
	
}*/
/*
public class Cluster {
public String name;
public HashSet<GrammarRuleRecord> ruleSet;
//public ArrayList<RuleInterval> intervals;
public Cluster(String name){
	this.name = name;
	ruleSet = new HashSet<GrammarRuleRecord>();
	//intervals = new ArrayList<RuleInterval>();	
}
public void addRule(GrammarRuleRecord rule){
	ruleSet.add(rule);
//	intervals.addAll(rule.getR0Intervals());
}
public HashSet<GrammarRuleRecord> getRules(){
	return ruleSet;
}
/*
public ArrayList<RuleInterval> getIntervals(){
	return intervals;
}
*/
}

