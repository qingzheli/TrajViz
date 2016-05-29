package edu.gmu.trajviz.model;

import java.util.ArrayList;
import java.util.HashMap;

import libsvm.LibSVM;
import net.sf.javaml.classification.Classifier;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.Instance;

public class Group {
	//public Integer label;
	public double left;
	public double right;
	public ArrayList<Integer> trajs;
	public double leftBoundary;
	public double rightBoundary;
	public ArrayList<SingleClass> family;
	public LibSVM svm;
	HashMap <String, ArrayList<Integer>> labelMap = new HashMap<String, ArrayList<Integer>>();
	//public ArrayList<Integer> member;
	public Group(){
		family = new ArrayList<SingleClass>();
		leftBoundary = -1;
		rightBoundary = -1;
		//label = name;
//		family = new ArrayList<Integer>();
		/*
		double max = 0;
		double min = Double.MAX_VALUE;
		
		for(int i = 0; i<family.size();i++){
			SingleClass singleClass = family.get(i);
			
			if(max<singleClass.right)
				max = singleClass.right;
			if(min>singleClass.left){
				min = singleClass.left;
			}
		}
		left = min;
		right = max;
		*/
		left = Double.MAX_VALUE;
		right = -2;
		svm = new LibSVM();
		trajs = new ArrayList<Integer>();
	}
	/*
	public void add(SingleClass sc){
		family.add(sc);
		if(left>sc.left)
			left = sc.left;
		if(right<sc.right)
			right = sc.right;
		for(int i = 0; i<sc.member.size();i++){
		String label = sc.label;
		  if(labelMap.containsKey(label)){
			  labelMap.get(label).add(i);
		  }
		  else
		  {
			  ArrayList<Integer> member = new ArrayList<Integer>();
			  member.add(i);
			  labelMap.put(label,member);	
			  
		  }
		}
	}
	*/
	public void add(SingleClass sc){
		family.add(sc);
		if(left>sc.left)
			left = sc.left;
		if(right<sc.right)
			right = sc.right;
		for(int i = 0; i<sc.member.size();i++){
		String label = sc.label;
		  if(labelMap.containsKey(label)){
			  labelMap.get(label).add(i);
		  }
		  else
		  {
			  ArrayList<Integer> member = new ArrayList<Integer>();
			  member.add(i);
			  labelMap.put(label,member);	
			  
		  }
		}
	}
	
	public String toString(){
		//return "["+left+","+right+"] boundary = ["+leftBoundary+","+rightBoundary+"],"+"Group family = "+family+"\n";
		return "["+left+","+right+"]" +SequiturModel.travelDistance.get(trajs.get(0))+"Group trajs = "+trajs+"\n";
	}
	public void add(int i) {
		trajs.add(i);
	}
	public void train() {
		Dataset datasets = new DefaultDataset();
		for(int i = 0; i<trajs.size();i++){
			
				Integer traj = trajs.get(i);
				Instance instance = SequiturModel.trajToInstance(traj);
				datasets.add(instance);
			}
		svm.buildClassifier(datasets);
	}
}
