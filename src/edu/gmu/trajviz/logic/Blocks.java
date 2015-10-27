package edu.gmu.trajviz.logic;

import java.util.ArrayList;

public class Blocks {
	public int n;   //means a area with n*n blocks
	public int nLat;
	public int nLon;
	public double latMin;
	public double lonMin;
	public double latMax;
	public double lonMax;
	public double latCut;
	public double lonCut;
	private int size;
	public ArrayList<Block> blocks;
	public Blocks(){
		n = 0;
		size = 0;
		blocks = new ArrayList<Block>();
	}
	public Blocks(int n, double laMin, double laMax, double loMin, double loMax){
		
		this.n = n;
		size = n*n;
		latMin = laMin;
		//deal with the max boundary problem, e.g. 8/1=8, however, the latId are from 0 to 7.
		latMax = laMax+0.0000001;
		lonMin = loMin;
		lonMax = loMax+0.0000001;
		latCut = (latMax - latMin)/n;
		lonCut = (lonMax - lonMin)/n;
		blocks = new ArrayList<Block>();
		for (int i=0;i<size;i++)
		{
			blocks.add(new Block(i,n,latCut,lonCut,latMin,lonMin));
		}
	}
	public void addPoint2Block(Location point){
		blocks.get(findBlockIdForPoint(point)).addPoint(point);
		
	}
	public int findBlockIdForPoint(Location point){
		 
		//System.out.println("latID: "+(Math.floor((point.latitude-latMin)/latCut))+" lonID: "+(Math.floor((point.longitude-lonMin)/lonCut)));
		if(point.latitude<-90||point.longitude<-180)
			return (int)(point.latitude);
		return (int)(Math.floor((point.latitude-latMin)/latCut))*n+(int)(Math.floor((point.longitude-lonMin)/lonCut));
	}
	public Block findBlockById(int id){
		return blocks.get(id);
	}
	public Block findBlockByLocation(Location point){
		return blocks.get(findBlockIdForPoint(point));
	}
	public void printBlockMap(){
		int i = 0;
		for (int k = 0; k<n; k++)
		{
			for (int j = 0; j<n; j++)
				{
					System.out.print(blocks.get(i).id+"\t");
					i++;
				}
		System.out.println();
		}
	}
	
	/*
	 *  compute distance between cells
	 */
	public double distance(Integer block1, Integer block2) {
		int latBlock1 = block1/n;
		int lonBlock1 = block1%n;
		int latBlock2 = block2/n;
		int lonBlock2 = block2%n;
		double distance = Math.sqrt((latBlock1-latBlock2)*(latBlock1-latBlock2)*latCut+(lonBlock1-lonBlock2)*(lonBlock1-lonBlock2)*lonCut);
		return distance;
	}
}
