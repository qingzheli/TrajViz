package edu.gmu.trajviz.logic;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Stack;

import edu.gmu.trajviz.model.Center;
import edu.gmu.trajviz.model.RoutePair;
import edu.gmu.trajviz.model.SequiturModel;
import edu.gmu.trajviz.util.Tools;

public class Block {
	public Blocks blocks;
	public int id;
	public int latId;
	public int lonId;
	public double latBlockMin;
	public double lonBlockMin;
	public double latBlockMax;
	public double lonBlockMax;
	public Block n; //north
	public Block e; //east
	public Block w; //west
	public Block s;
	public Block ne;
	public Block nw;
	public Block se;
	public Block sw;
	public ArrayList<Block> nearbyBlocks;
	public ArrayList<Location> points;
	public ArrayList<Center> centers;
	public ArrayList<String> residualSet;
	public int size;  // a huge area divided into size*size blocks
	public Block(){
		// just an empty block;
		id = 0;
		latId = 0;
		lonId = 0;
		size = 0;
		points = new ArrayList<Location>();
		centers = new ArrayList<Center>();
		residualSet = new ArrayList<String>();
			
	}
	/* doesn't use for now
	public Block(int latId, int lonId, int size){
		this.latId = latId;
		this.lonId = lonId;
		this.size = size;
		id = latId*size+lonId;
		points = new ArrayList<Location>();
	}
	*/
	public Block(int id, int size,double latCut, double lonCut,double latMin, double lonMin,Blocks blks){
		this.blocks = blks;
		this.id = id;
		this.size = size;
		this.latId = id/size;
		this.lonId = id%size;
		latBlockMin = latMin + this.latId*latCut;
		latBlockMax = latBlockMin + latCut;
		lonBlockMin = lonMin + this.lonId*lonCut;
		lonBlockMax = lonBlockMin + lonCut;
		points = new ArrayList<Location>();
		centers = new ArrayList<Center>();
		residualSet = new ArrayList<String>();
		SequiturModel.allSubseq.put(id, new ArrayList<String>());
	//	setNearbyBlocks();
		
	}
	public void addPoint(Location point){
		points.add(point);
	}
	public void addCenter(Center center){
		centers.add(center);
		
	}
	/*
	public void addSubtrajectory(ArrayList<Double> subx, ArrayList<Double> suby){
		
	}
	*/
	public void findAnomaly(){
		
		if(centers.size()>0){
			
		
		double r = SequiturModel.R;
	//	System.out.println("r = "+r);
		
	//	SequiturModel.allSubseq.put(id, new ArrayList<String>());
		 
	//	int step = (int)(len*SequiturModel.minLink);
	/*	int len = SequiturModel.minBlocks;
		ArrayList<String> residualSet = new ArrayList<String>();
		for(int i = 0; i<centers.size(); i++){
			Center center = centers.get(i);
			System.out.println("anomaly_center = "+center);
			for(int s = center.s; s<=center.e; s= s+1){
				String subseqId = "T"+center.traj+"S"+s+"L"+len;
			//	System.out.println("anomaly_center_subseqId = "+subseqId);
				SequiturModel.allSubseq.get(id).add(subseqId);
				residualSet.add(subseqId);
			}
		}
		*/
		 for(int i = 0; i<SequiturModel.allSubseq.get(id).size()&&residualSet.size()>0; i++){
			 String subseq1 = SequiturModel.allSubseq.get(id).get(i);
			 if(residualSet.contains(subseq1)){
				// int[] sub1 = Tools.parseTrajId(subseq1); 
				 for(int j = 0; j<SequiturModel.allSubseq.get(id).size(); j++){
					 String subseq2 = SequiturModel.allSubseq.get(id).get(j);
				//	 int[] sub2 = Tools.parseTrajId(subseq2);
				//	 if(sub1[0]!=sub2[0]){
					 if(!Tools.isTrivialMatch(subseq1, subseq2)){					 
						 RoutePair pair = new RoutePair(subseq1,subseq2,r);
						
						 if(pair.dist<=r){
							
							 residualSet.remove(subseq1);
							 residualSet.remove(subseq2);
							 setisAnomalyFalse(subseq1);
							 setisAnomalyFalse(subseq2);
							 
						 }
				 }
				 
				 }
			 }
		 }
		
		
		
		 
		 /*
		  * check nearby blocks
		  */
		 /*
		 Stack<Integer> removeCandidates = new Stack<Integer>();
		 if(residualSet.size()>0){
			 for(int i=0; i<residualSet.size(); i++){
				 
				 	
				 	right
				 
			 }
		 }
		 */
		}
	}
	public void checkNearby(){
		double r = SequiturModel.R;
		 for(int index = 0; index<nearbyBlocks.size()&&residualSet.size()>0;index++){
			 Block nearbyBlock = nearbyBlocks.get(index);
			 for(int i = 0; i<SequiturModel.allSubseq.get(nearbyBlock.id).size()&&residualSet.size()>0; i++){
				 String subseq1 = SequiturModel.allSubseq.get(nearbyBlock.id).get(i);
				 
				//	 int[] sub1 = Tools.parseTrajId(subseq1);
					 Stack<Integer> removeCandidates = new Stack<Integer>();
					 for(int j = 0; j<residualSet.size(); j++){
						 String subseq2 = residualSet.get(j);
				//		 int[] sub2 = Tools.parseTrajId(subseq2);
			//			 if(sub1[0]!=sub2[0]){
						 if(!Tools.isTrivialMatch(subseq1, subseq2)){					 
							 RoutePair pair = new RoutePair(subseq1,subseq2,r);
							
							 if(pair.dist<=r){
								
								 nearbyBlock.residualSet.remove(subseq1);
								 //residualSet.remove(subseq2);
								 removeCandidates.push(j);
								 setisAnomalyFalse(subseq1);
								 setisAnomalyFalse(subseq2);
								 SequiturModel.count_works = SequiturModel.count_works+1;
							 }
					 }
						 
					 
					 }
					 while(!removeCandidates.isEmpty()){
						 residualSet.remove(removeCandidates.pop());
					 }
				 
			 }
		 }
	}
	public ArrayList<Block> setNearbyBlocks() {
		nearbyBlocks = new ArrayList<Block>();
		int nid = id-size,sid =id+size,wid = id-1,eid = id+1,neid = id-size+1,nwid = id-size-1,seid = id+size+1,swid = id+size-1;
		if(latId==0)
		{
			nid = -1;
			neid = -1;
			nwid = -1;
		}
		if(latId == blocks.nLat-1){
			sid = -1;
			seid = -1;
			swid = -1;
		}
		if(lonId==0){
			wid = -1;
			swid = -1;
			nwid = -1;
		}
		if(lonId == blocks.nLon-1){
			eid = -1;
			seid = -1;
			neid = -1;
		}
	//	System.out.println("id = "+id);
		if(eid !=-1){
			// eid = id+1;
			
			 e = blocks.findBlockById(eid);
			 nearbyBlocks.add(e);
		//	 centers.addAll(e.centers);
			 if(nid!=-1){
			//	neid = id - size+1;
				ne = blocks.findBlockById(neid);
				nearbyBlocks.add(ne);
		//		centers.addAll(ne.centers);
			 }
			 if(sid!=-1){
			 //	seid = id+size+1;
				se = blocks.findBlockById(seid);
				nearbyBlocks.add(se);
		//		centers.addAll(se.centers);
			 }
		}
		if(wid != -1){
		//	wid = id - 1;
			w = blocks.findBlockById(wid);
			nearbyBlocks.add(w);
		//	centers.addAll(w.centers);
			if(nid!=-1){
			//	nwid = id-size-1;
				nw = blocks.findBlockById(nwid);
				nearbyBlocks.add(nw);
	//			centers.addAll(nw.centers);
			}
			if(sid!=-1){
			//	swid = id+size-1;
				sw = blocks.findBlockById(swid);
				nearbyBlocks.add(sw);
	//			centers.addAll(sw.centers);
			}
		}
		if(nid!=-1){
		//	nid = id-size;
			n = blocks.findBlockById(nid);
			nearbyBlocks.add(n);
	//		centers.addAll(n.centers);
		}
		if(sid!=-1){
		//	sid = id + size;
			s = blocks.findBlockById(sid);
			nearbyBlocks.add(s);
	//		centers.addAll(s.centers);
		}
	//	blocks.printBlockMap();
	//	System.out.println("nearbyBlocks.size = "+nearbyBlocks.size());
	//	System.out.println(nearbyBlocks);
		
		return nearbyBlocks;
	}
	public static void setisAnomalyFalse(String subseq) {
		int[] sub = Tools.parseTrajId(subseq);
		int length = 0;
		if((sub[1]+sub[2])==SequiturModel.oldtrajX.get(sub[0]).size())
			length = sub[2];
		else
			length = 2;
		for(int s = 0; s<sub[2]; s++){
	//	for(int s = 0; s<length; s++){
			int i = sub[1]+s;
			SequiturModel.isAnomaly.get(sub[0]).set(i, false);
		}
		
	}
	public void findHierarchicalMotifs() {
		if(centers.size()>0){
		SequiturModel.allSubseq.put(id, new ArrayList<String>());
		
		PriorityQueue<RoutePair> mergablePair = new PriorityQueue<RoutePair>();
		HashMap<String, Cluster> clusterMap = new HashMap<String, Cluster>();
		/*
	     SequiturModel.allSubseq = new HashMap<Integer,ArrayList<String>>();
		  for(int len = minLength; len<=maxLength; len = minLength +len){
			  allTrajClusters.put(len, new  ArrayList<Cluster>());
			  allSubseq.put(len, new ArrayList<String>());
			  mergablePair.put(len, new PriorityQueue<RoutePair>());
			  clusterMaps.put(len, new HashMap<String, Cluster>());
		  }
		  */
		//iterator all centers subtract subseqId
		double r = SequiturModel.R;
		int len = SequiturModel.minBlocks; 
		int step = (int)(len*SequiturModel.minLink);
		/*
		if(centers.size()>0){
			len = centers.get(0).e-centers.get(0).s+1;
		}
		else
			len = -1;
			*/
		  for(int i = 0; i<centers.size(); i++){
			 Center center = centers.get(i);
		//	 System.out.println("center = "+center);
			//	  for(int s = center.s; s<=center.e; s++ ){
			 for(int s = center.s; s<=center.e; s=s+step ){
					  String subseqId = "T"+center.traj+"S"+s+"L"+len;
					  SequiturModel.allSubseq.get(id).add(subseqId);
				  }
			 
		  }
		  
			  for(int i = 0; i<SequiturModel.allSubseq.get(id).size(); i++){
				 String subseq1 = SequiturModel.allSubseq.get(id).get(i);
				//  System.out.println("Threshold = " + R);
				  for(int j = i+1; j<SequiturModel.allSubseq.get(id).size();j++){
					//  System.out.println("R" + R);
					  
					  String subseq2 = SequiturModel.allSubseq.get(id).get(j);
					  if(Tools.parseTrajId(subseq1)[0]!=Tools.parseTrajId(subseq2)[0]){
					  RoutePair pair = new RoutePair(subseq1,subseq2,r);
					  if(pair.dist<=r){
					 // if(!pair.isTrivial&&pair.dist<=R){
					//	  System.out.println(R +" Pair: "+pair);
						  mergablePair.add(pair);
					  }
					  }
				  }
				  
			  }
	//		 System.out.println("PQ = "+mergablePair);
		//	  System.out.println("len = "+ len+ "size = "+ mergablePair.size()	 + " peak = "+mergablePair.peek() );
		  
		  
		  
		  
		  //  start hierarchical clustering
		  
		  
			  PriorityQueue<RoutePair> pq = mergablePair;
			  HashMap<String, Cluster> findCluster = clusterMap;
			  while(!pq.isEmpty()){
				  RoutePair pair = pq.poll();
				  if(!findCluster.containsKey(pair.r1)&&!findCluster.containsKey(pair.r2)) // if neither belongs to a cluster, create new cluster
				  {
					Cluster cluster = new Cluster(len);
					int[] t1 = Tools.parseTrajId(pair.r1);
					int[] t2 = Tools.parseTrajId(pair.r2);
					cluster.add(t1[0], t1[1]);
					cluster.add(t2[0], t2[1]);
				  findCluster.put(pair.r1, cluster);
				  findCluster.put(pair.r2, cluster);
				  SequiturModel.allTrajClusters.get(id).add(cluster);
		//		  System.out.println("Clusters "+(SequiturModel.allTrajClusters.get(id).size()-1)+"\n"+cluster);
				  }
				  else if(!findCluster.containsKey(pair.r1)&&findCluster.containsKey(pair.r2)){    // if r1 is not in cluster r2 is in a cluster
					  int[] t1 = Tools.parseTrajId(pair.r1);
					  Cluster cluster = findCluster.get(pair.r2);
					  cluster.add(t1[0], t1[1]);
					  findCluster.put(pair.r1, cluster);
					  
				  }
				  else if(findCluster.containsKey(pair.r1)&&!findCluster.containsKey(pair.r2)){   // if r2 is not in cluster r1 is in a cluster
					  int[] t2 = Tools.parseTrajId(pair.r2);
					  Cluster cluster = findCluster.get(pair.r1);
					  cluster.add(t2[0], t2[1]);
					  findCluster.put(pair.r2, cluster);
					  
				  }
				  else  
				  {
					 Cluster c1 = findCluster.get(pair.r1);
					 Cluster c2 = findCluster.get(pair.r2);
					 if(c1!=c2){
						// System.out.println("Before Merege C1: "+c1.trajIds);
						// System.out.println("Before Merege C2: "+c2.trajIds); 
						 c1.merge(c2,r,findCluster);
						 if(c1.getRoutes().size()<1){
							 
							 SequiturModel.allTrajClusters.get(id).remove(c1);
							 
						 }
						 if(c2.getRoutes().size()<1){
							 SequiturModel.allTrajClusters.get(id).remove(c2);
						 }
						 
				/*		
						 if(c1!=null)
							 System.out.println("After  Merege C1: "+c1.trajIds);
						 if(c2!=null)
							 System.out.println("After  Merege C2: "+c2.trajIds);
				*/		
					 }
				  }
				  
			  }
			  
			  
		  
		  
		}  
	}
	public void setResidualSet() {
		int len = SequiturModel.minBlocks;
		for(int i = 0; i<centers.size(); i++){
			Center center = centers.get(i);
//			System.out.println("anomaly_center = "+center);
			for(int s = center.s; s<=center.e; s= s+1){
				String subseqId = "T"+center.traj+"S"+s+"L"+len;
		//		System.out.println("set_center_subseqId = "+subseqId);
				SequiturModel.allSubseq.get(id).add(subseqId);
				residualSet.add(subseqId);
			}
		}
	
	
				
	}
	public String toString(){
		return id+"("+latId+","+lonId+") ";
	}
}
