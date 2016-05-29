package edu.gmu.trajviz.model;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import edu.gmu.trajviz.logic.Blocks;
import edu.gmu.trajviz.logic.Cluster;
import edu.gmu.trajviz.logic.Location;
import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.util.Tools;

public class Evaluation{
public static HashMap<Integer, Double> silCoefMap(HashMap<Integer, ArrayList<Cluster>> allMotifs){
	HashMap<Integer, Double> coefMap = new HashMap<Integer,Double>();
	Iterator it = allMotifs.entrySet().iterator();
	while(it.hasNext())
	{
		Entry<Integer, ArrayList<Cluster>> entry = (Entry<Integer, ArrayList<Cluster>>) it.next();
		Integer len = entry.getKey();
		ArrayList<Cluster> clusters = entry.getValue();
		Double avgSilcoef = getAvgSilcoef(clusters);
		coefMap.put(len, avgSilcoef);
	}
	
	//for
	return coefMap;
}
private static double getAvgSilcoef(ArrayList<Cluster> clusters) {
	double[] result = new double[4];  
	StringBuffer sb = new StringBuffer();
	ArrayList<Double> allSilcoef = new ArrayList<Double>();
	for (int i=0;i<clusters.size();i++){   // iterate all motifs
		Cluster cluster = clusters.get(i);
		Iterator it = cluster.trajIds.iterator();
		while(it.hasNext()){
			String name = (String) it.next();
			double a = getAvgIntraDistance(name,cluster);	
			double b = getMinInterDistance(i,name,clusters);
		//	System.out.println("a = "+a);
		//	System.out.println("b = "+b);
			double silcoef = (b-a)/Math.max(a, b);
			allSilcoef.add(silcoef);
	//	double avgIntraDistance = getAvgPairDistances(clusters.get(i)); //intra-distances
		}
	}
	System.out.println("=============================================================");
	System.out.println(allSilcoef);
	double avgSilcoef = avg(allSilcoef);
	return avgSilcoef;
}
	/*
	// evaluate distances inter-rules
	
		ArrayList<Double> pairInterDistances = new ArrayList<Double>();
		for(int i = 0; i<clusters.size();i++){
			for(int j=i+1; j<clusters.size();j++){
			
			 pairInterDistances.add(getMinDist(clusters.get(i),clusters.get(j)));
		
		}
	
	}
	
//	System.out.println("average distances among all motifs: "+avg(allDistances));
//	System.out.println("average standard deviation among all motifs: "+avg(allStdDev));
	//sb.append(avg(allDistances)+","+avg(allStdDev)+"\n");
	
	//return sb.toString();
	result[0] = avg(allDistances);  // avg intra distances
	result[1] = avg(allStdDev);  // avg std. dev. intra distances
	result[2] = avg(allMinimalInterDistances);  // avg minimal inter distances
	ArrayList<Double> silhouetteCoefficients = new ArrayList<Double>();
	for (int i = 0; i<routes.size(); i++)
	{
		
		double sc;
		
		if(allMinimalInterDistances.size()>0&&allDistances.get(i)<allMinimalInterDistances.get(i)){
	//		if(allDistances.get(i)<allMinimalInterDistances.get(i)){

			sc = 1 - allDistances.get(i)/allMinimalInterDistances.get(i);
		}
		else{
			if(allDistances.get(i)==0||allMinimalInterDistances.size()==0)
				sc = 1;
			else
			sc = allMinimalInterDistances.get(i)/allDistances.get(i) - 1;
			
		}
//		System.out.println("compare: "+allDistances.get(i)+"/"+allMinimalInterDistances.get(i)+" = "+sc);
		silhouetteCoefficients.add(sc);
	}
	result[3] = avg(silhouetteCoefficients);
	return result;
  }
  */
private static double getMinInterDistance(int i, String name, ArrayList<Cluster> clusters) {
	Route route = clusters.get(i).getRoutes().get(name);
	double min = Double.MAX_VALUE;
	for(int k = 0; k<clusters.size(); k++){
		double dist = 0;
		if(i!=k){
			Iterator it = clusters.get(k).trajIds.iterator();
			
			while(it.hasNext()){
				String a = (String) it.next();
				 dist = dist+Tools.routeEuDist(route,clusters.get(k).getRoutes().get(a));
				
			}
			double avgDist = dist/clusters.get(k).getRoutes().size();
			if(avgDist<min)
					min = avgDist;
		}
	}
	return min;
}
private static double getAvgIntraDistance(String name, Cluster cluster) {
	double sum = 0;
	int count = 0;
	Iterator it = cluster.trajIds.iterator();
	while(it.hasNext()){
		String next = (String) it.next();
		if (!name.equals(next)){
			sum = sum+Tools.routeEuDist(cluster.getRoutes().get(name),cluster.getRoutes().get(next));
		//	System.out.println("Sum = "+sum);
			count++;
		}
	}
//	System.out.println("Count = "+count);
//	System.out.println("Size  = "+cluster.getSize());
	double avg = sum/count;
	return avg;
}
private static Double getMinDist(Cluster c1, Cluster c2) {
	Double min = Double.MAX_VALUE;
	
	for(int i = 0; i<c1.getRoutes().size();i++)
		for(int j = 0; j<c2.getRoutes().size();j++){
			double dist = Tools.routeEuDist(c1.getRoutes().get(i),c2.getRoutes().get(j));
			if(dist<min){
				min = dist;
			}
		}
	return min;
}
/*
private static Double getAvgPairDistances(Cluster c) {
	ArrayList<Double> pairDistances = new ArrayList<Double>();
	HashMap<String, Route> c.getRoutes();
	ArrayList<Route> routes = 
	Iterator it1 = routes.entrySet().iterator();
	Iterator it2 = routes.entrySet().iterator();
	for(int i = 0; i<routes.size(); i++){
		for(int j = i+1; j<routes.size();j++){
			pairDistances.add(Tools.routeEuDist(routes.get(i),routes.get(j)));
		}
	}
	double sum = 0;
	for(int i = 0; i<pairDistances.size();i++){
		sum = sum+pairDistances.get(i);
	}
	
	return sum/pairDistances.size();
}
*/
/*
 * requires: two routes has the same length
 * return avg Euclidean distance between two routes
 */

private static Double avg(ArrayList<Double> list){
	Double sum= new Double(0);
	for (int i = 0; i< list.size(); i++){
		if(Double.isNaN(list.get(i)))
			throw new NullPointerException();
		sum = sum + list.get(i);
//		System.out.print(list.get(i)+" ");
	}
//	System.out.println();
	
//	System.out.println("sum = "+sum+" avg = "+sum/list.size());
	if(list.size()>0)
	return sum/list.size();
	else
		return -88888.0;
}
private static double dev(ArrayList<Double> pairwiseDistances) {
	double avg = avg(pairwiseDistances);
	double sum = 0;
	for (int i=0; i<pairwiseDistances.size(); i++)
		sum = sum + (pairwiseDistances.get(i)-avg)*(pairwiseDistances.get(i)-avg);
	return Math.sqrt(sum);
}

}