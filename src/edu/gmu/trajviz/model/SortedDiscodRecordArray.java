package edu.gmu.trajviz.model;
public class SortedDiscodRecordArray{
	private DiscordRecord[] sortedArray;
	//private Double[] sortedArray;
	private int length;
	public SortedDiscodRecordArray(int k){
		sortedArray = new DiscordRecord[k];
		length = k;
		//sortedArray = new Double[k];
	}

	public void add(DiscordRecord record){
		
		for(int i = 0; i<sortedArray.length; i++){
			if(record.hasOverlap(sortedArray[i]))
			{
				if(record.compareTo(sortedArray[i])>0)
					{
						sortedArray[i] = record;
						return;
					}
				else return;
				
			}
			if(sortedArray[i]==null||record.compareTo(sortedArray[i])>0)
			{
				for(int j = sortedArray.length-1; j>i; j-- )
					sortedArray[j]=sortedArray[j-1];
				sortedArray[i] = record;
				return;
			}
		}
		
	}
	/*
public int findDiscordByTrajId(String id){
		
		for(int i = 0; i<sortedArray.length; i++){
			if(TrajDiscords.getTrajectory(sortedArray[i].getEndPosition())==id){
				return i+1;
			}
			
		}
		return -1;
		
	}
	*/
	@Override
	public String toString(){
		StringBuffer sb = new StringBuffer();
		for(int i = 0; i<sortedArray.length;i++)
		 sb.append(sortedArray[i]);
		return sb.toString();
	}
	
	/*
	public void add(Double record){
		
		for(int i = 0; i<sortedArray.length; i++){
			
			if(sortedArray[i]==null||record.compareTo(sortedArray[i])<0)
			{
				for(int j = sortedArray.length-1; j>i; j-- )
					sortedArray[j]=sortedArray[j-1];
				sortedArray[i] = record;
				for(int k = 0; k< sortedArray.length; k++)
					System.out.println(i+":"+sortedArray[k]);
				return;
			}
		}
		
	}
	
	public static void main(String[] args){
		Double[] a = {1.0,3.0,9.0,5.0,7.0,8.0,2.0,4.0,1.0,6.0};
		SortedArray sa = new SortedArray(7);
		for (int i = 0; i< a.length;i++){
			sa.add(a[i]);
		}
		for(int i = 0; i< a.length; i++)
		System.out.println(a[i]);
		for(int i = 0; i< sa.sortedArray.length; i++)
			System.out.println(i+":"+sa.sortedArray[i]);
	}
*/
}
