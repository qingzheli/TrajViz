package edu.gmu.trajviz.model;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.HashMap;

import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.util.Tools;

public class Center {
public Integer blockid;
public int traj;
public int s;
public int e;
public double x;
public double y;
public HashMap<Integer, Double> neighborMap;

public Center(int trajectory, int start, int end, double x, double y){
	traj = trajectory;
	s = start;
	e = end;
	this.x = x;
	this.y = y;
	blockid = -1;
	neighborMap = new HashMap<Integer,Double>();
}

public Center(int trajectory, int start,double x, double y){
	this.x = x;
	this.y = y;
	s =start;
	traj = trajectory;
	blockid=-1;
	
}
public double squareDistance(Center c2){
	return (x-c2.x)*(x-c2.x)+(y-c2.y)*(y-c2.y);
}
/*
public boolean isClose(Center c2, double threshold, int slidePoint){
	if(this.traj==c2.traj)
		return false;
	boolean answer = false;
//	int point = (int)(length-l)/2;
	double minDist = Double.MAX_VALUE;
	if(slidePoint<=1)
		slidePoint = 0;
	for(int x = s; x<=s+slidePoint; x++){
		Route r1 = Tools.getSubroute(traj, x, l);
		Route r2 = Tools.getSubroute(c2.traj, c2.s, l+slidePoint);
		double dist = Tools.RouteSITED(r1, r2);
		if(dist<=threshold&&minDist>dist){
			minDist = dist;
			answer = true;
		}
	}
	
	if(answer){
		addNeighbor(c2,minDist);
		c2.addNeighbor(this, minDist );
	}
	return answer;
	
}
private void addNeighbor(Center c2, double minDist) {
	if(neighborMap.containsKey(c2.traj)){
		if(minDist<neighborMap.get(c2.traj)){
			neighborMap.put(c2.traj, minDist);
			
		}
	}
	else
		{
			neighborMap.put(c2.traj, minDist);	
		}
	
}
*/
public String toString(){
	NumberFormat formatter = new DecimalFormat("#0.00");
	return "T"+traj+"S"+s+"E"+e+"("+formatter.format(x)+","+formatter.format(y)+")";
	//return blockid.toString();
	/*
	if(e>0)
	return blockid+"= ("+x+","+y+")  s"+s+",e"+e;
	else
	//	return blockid+"= ("+x+","+y+")  s"+s+",e"+e;
		return blockid.toString();
		*/
}

}
