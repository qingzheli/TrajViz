package edu.gmu.trajviz.sax.datastructures;

public class Interval {
	private int startIdx;
	private int endIdx;
	public Interval(int s, int e){
		startIdx = s;
		endIdx = e;
	}
	public int getStartIdx(){
		return startIdx;
	}
	public int getEndIdx(){
		return endIdx;
	}
	@Override
	public String toString(){
		return startIdx+"-"+endIdx;
	}

}
