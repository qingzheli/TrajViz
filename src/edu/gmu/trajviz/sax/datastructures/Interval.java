package edu.gmu.trajviz.sax.datastructures;

public class Interval implements Comparable<Interval> {
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

	@Override
	public int compareTo(Interval o) {
		if(startIdx==o.startIdx&&endIdx==o.endIdx)
			return 0;
		else{
			if(startIdx==o.startIdx)
				return endIdx>o.endIdx?1:(-1);
			else
				return startIdx>o.startIdx? 1:(-1);
		}
	}
	public static Interval mergeInterval(Interval first, Interval second){
		//System.out.println("first   : "+first);
		//System.out.println("second  : "+second);
		if((first.startIdx>=second.startIdx&&first.startIdx<=(second.endIdx+1))||(second.startIdx>=first.startIdx&&second.startIdx<=(first.endIdx+1))){
			Interval result = new Interval(Math.min(first.startIdx, second.startIdx),Math.max(first.endIdx, second.endIdx));
		//System.out.println("result  : "+result);
			return result;
		}
		else
			return null;
		
		
		
		/*
		if((first.endIdx)==(second.startIdx-1))
			return new Interval(first.startIdx, second.endIdx);
		else if(second.endIdx == (first.startIdx-1))
			return new Interval(second.startIdx,first.endIdx);
		else
		return null;
		*/
	}

}
