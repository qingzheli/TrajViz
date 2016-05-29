package edu.gmu.trajviz.model;

import java.util.ArrayList;

public class SingleClass implements Comparable{
public String label;
public double left;
public double right;
//public int leftBoundary;
//public int rightBoundary;
//public ArrayList<Integer> family;
public ArrayList<Integer> member;
public SingleClass(String name, ArrayList<Integer> trajs){
	member = trajs;
	label = name;
//	family = new ArrayList<Integer>();
	double max = 0;
	double min = Double.MAX_VALUE;
	
	for(int i = 0; i<member.size();i++){
		int traj = member.get(i);
		double dist = SequiturModel.travelDistance.get(traj);
		if(max<dist)
			max = dist;
		if(min>dist){
			min = dist;
		}
	}
	left = min;
	right = max;
}
@Override
public int compareTo(Object o) {
	SingleClass c = (SingleClass)o;
	if(left<c.left)
		return -1;
	if(left==c.left)
		return 0;
	return 1;
}
@Override
public String toString() {
	return label+"["+left+","+right+"]";
}
}
