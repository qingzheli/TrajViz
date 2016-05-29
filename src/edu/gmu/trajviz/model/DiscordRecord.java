package edu.gmu.trajviz.model;

import edu.gmu.trajviz.util.Tools;

public class DiscordRecord implements Comparable<DiscordRecord> {
private String id;
private int length;
private int startPos;
private int endPos;
private Double minDist;
private int[] sub; 
public DiscordRecord(int start, int end, double currentMinDist) {
	length = end-start+1;
	id = start+","+length;
	this.startPos =start;
	this.endPos = end;
	this.minDist = currentMinDist;
	}
public DiscordRecord(String trajId, double currentMinDist) {
	sub = Tools.parseTrajId(trajId);
	length = sub[2];
	id = trajId;
	this.startPos =sub[1];
	this.endPos = sub[1]+length-1;
	this.minDist = currentMinDist;
	}
public int getLength(){
return this.length;
}
public String getId(){
	return this.id;
	
}
public int getStartPosition(){
	return this.startPos;
}
public int getEndPosition(){
	return this.endPos;
}
public int getMinDistance(){
	return this.getMinDistance();
}
@Override
public int compareTo(DiscordRecord o) {
	if(o==null)
		return -1;
	return this.minDist.compareTo(o.minDist);
}
public boolean hasOverlap(DiscordRecord record) {
	
	
	if(record == null)
		return false;
	if(length!=record.length){
		throw new IllegalArgumentException("length1 != length2");
	}
	if(sub[0]==record.sub[0]&&((startPos>=record.startPos&&startPos<record.endPos)||(endPos>record.startPos&&endPos<=record.endPos)))
	return true;
	else
		return false;
}
@Override
public String toString(){
//	return id+": ["+startPos+","+endPos+"], minDist = "+ minDist+ " Trajectory: "+TrajDiscords.getTrajectory(endPos)+", length = "+length+"\n";
	return id+","+minDist+"\n";
}
}
