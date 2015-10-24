package edu.gmu.trajviz.logic;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.PriorityQueue;

import edu.gmu.trajviz.gi.GrammarRules;
import edu.gmu.trajviz.model.SequiturModel;

public class RuleDistanceMatrix {

public double[][] matrix;
private double minDistance;
private int[] minPair = new int[2];
private GrammarRules rules;
public RuleDistanceMatrix(){}
public ArrayList<Integer> filter;
public PriorityQueue<PairDistance> pq; 
private PairDistanceComparator comparator;
@SuppressWarnings("unchecked")
public RuleDistanceMatrix(Blocks blocks, GrammarRules rules, ArrayList<Integer> filter){
	
	this.filter = filter;
	/*
	for (int i=0; i<filter.size();i++)
	System.out.println(i+" : "+filter.get(i)+" ");
	*///System.out.println();
	matrix = new double[filter.size()][filter.size()];
	//matrix[0][0] =100;// Double.MAX_VALUE;
	comparator = new PairDistanceComparator();
	pq= new PriorityQueue<PairDistance>(filter.size()*filter.size()/2, comparator);
	minDistance = 100;//Double.MAX_VALUE;
	minPair[0] = 0;
	minPair[1] = 0;
	//int line = 0;
	//int col = 0;
	for(int i = 0; i<filter.size();i++ )
		for(int j = i+1; j<filter.size();j++){
			if(rules.get(filter.get(i)).frequencyInR0()>2&&rules.get(filter.get(j)).frequencyInR0()>2)
			{
				String rule1 = rules.getRuleRecord(filter.get(i)).getExpandedRuleString();
				String rule2 = rules.getRuleRecord(filter.get(j)).getExpandedRuleString();
		//	matrix[0][j] = 100;//Double.MAX_VALUE;
		//	matrix[i][0] = 100;//Double.MAX_VALUE;
			matrix[i][j] = avgDTWDistance(blocks, toArrayList(rule1), toArrayList(rule2));
			matrix[j][i] = matrix[i][j];
			if(matrix[i][j]>0&&matrix[i][j]<SequiturModel.MINLINK){//&&matrix[i][j]<minDistance){
				pq.add(new PairDistance(i,j,matrix[i][j]));
			/*
				minDistance = matrix[i][j];
				minPair[0] = i;
				minPair[1] = j;
				*/
			}
			}
			//if(rules.get(i).)
		}
	this.rules = rules;
	printMatrix(matrix);
}
private void printMatrix(double[][] matrix) {
	NumberFormat formatter = new DecimalFormat("#0.00");     
/*	
	for(int i=0; i<matrix.length; i++)
	{
		for(int j = 0; j<matrix[0].length;j++)
		{
			System.out.print(formatter.format(matrix[i][j])+"   ");
		}
	System.out.println();	

	}
	*/
	System.out.println("minDistance"+pq.peek().getDistance());
	System.out.println("minPair: "+pq.peek().getLine()+", "+filter.get(pq.peek().getLine())+";"+pq.peek().getCol()+","+filter.get(pq.peek().getCol()));
	System.out.println("rule1 :"+rules.getRuleRecord(filter.get(pq.peek().getLine()))+" Expand String: "+rules.get(filter.get(pq.peek().getLine())).getExpandedRuleString());
	System.out.println("rule2 :"+rules.getRuleRecord(filter.get(pq.peek().getCol()))+" Expand String: "+rules.get(filter.get(pq.peek().getCol())).getExpandedRuleString());
	System.out.println("Matrix size: "+filter.size());
}
public double[][] getMatrix(){
	return matrix;
}
public double getMinDistance(){
	return minDistance;
}
public int[] getMinPair(){
	return minPair;
}
public static ArrayList<Integer> toArrayList(String rule) {
	String[] strArray = rule.split(" ");
	ArrayList<Integer> al = new ArrayList<>(); 
	for (int i = 0; i<strArray.length;i++){
		al.add(Integer.valueOf(strArray[i]));
	}
	return al;
}

private double avgDTWDistance(Blocks blocks,ArrayList<Integer> s,
		ArrayList<Integer> t) {
	
//	System.out.print("s::::::::::size:"+s.size());
	/*
	for(int i=0; i<s.size();i++)
		System.out.print(" "+s.get(i));
		*/
//	System.out.println();
//	System.out.print("t::::::::::size:"+t.size()+"   ");
	/*
	for(int i=0; i<t.size();i++)
		System.out.print(" "+t.get(i));
		*/
//	System.out.println();
	int n = s.size();
	int m = t.size();
	double[][] DTW = new double[n+1][m+1];
	double cost = 0;
	for(int i=0;i<n;i++)
		DTW[i+1][0]=Double.MAX_VALUE;
	for(int i=0;i<m;i++)
		DTW[0][i+1] = Double.MAX_VALUE;
	DTW[0][0] = 0;
	for (int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			cost = blocks.distance(s.get(i),t.get(j));
		//	System.out.println("cost_"+i+","+j+": "+cost);
			DTW[i+1][j+1]=cost+minimum(DTW[i][j+1],		// insertion
									   DTW[i+1][j], 	// deletion
									   DTW[i][j]);	// match
		}
	}
//	System.out.println("DTW:::::"+DTW[n][m]);
	int step = 1;
	int x = n;
	int y = m;
	while(!((x==1)&&(y==1))){
		step = step + 1;
		switch(min(DTW[x-1][y-1],DTW[x-1][y],DTW[x][y-1])){
		case 1: x--; y--; break;
		case 2: x--; break;
		case 3: y--; break;
		default: System.out.println("Error!!!!");
		}
		
	}
//	System.out.println("step: "+step);
	double avg = DTW[n][m]/step;
//	System.out.println("avgDTW:::::"+avg);
	return avg;
}

private int min(double d, double e, double f) {
	if(d<=e&&d<=f)
		return 1;
	if(e<=d&&e<=f)
		return 2;
	if(f<=e&&f<=d)
		return 3;
	return 0;
}

private double minimum(double a, double b, double c) {
	return Math.min(Math.min(a, b), Math.min(b, c));
}
}