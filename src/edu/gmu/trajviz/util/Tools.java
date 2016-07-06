package edu.gmu.trajviz.util;

import edu.gmu.trajviz.logic.Route;
import edu.gmu.trajviz.model.Center;
import edu.gmu.trajviz.model.SequiturModel;

public class Tools {
	
	public static double pointEuDist(double x1, double y1 , double x2, double y2) {
		
	//	System.out.println("x1,y1,x2,y2 : "+x1+","+y1+","+x2+","+y2);
		return Math.sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	}
	
	/*
	 *  if r1 and r2 are close to each other (all euDist between pairs of points are <= R
	 *  	 return avg euDist;
	 *  else
	 *  	return Double.MAXVALUE;
	 */
	/*
	public static Double closeRouteEuDist(Route r1, Route r2, double threshold) {
		if(r1.getLats().size()!=r2.getLats().size())
		{
			throw new IllegalArgumentException("r1.size != r2.size    ->"+r1.getLats().size()+ " : " +r2.getLats().size());
		}
		double dist = 0;
		for(int i = 0; i<r1.getLats().size(); i++){
			double pointDist = Tools.euDist(r1.getLats().get(i), r1.getLons().get(i), r2.getLats().get(i), r2.getLons().get(i));
			if(pointDist<=threshold)
				dist = dist+ pointDist;
			else
				{
				//	System.out.println(R+"   PointDist: "+pointDist);
					return Double.MAX_VALUE;
				}
		}
		return dist/r1.getLats().size();
	}
	*/
	public static Double routeSqrEuDist(Route r1, Route r2,double r) {
		if(r1.getLats().size()!=r2.getLats().size())
		{
			throw new IllegalArgumentException("r1.size != r2.size    ->"+r1.getLats().size()+ " : " +r2.getLats().size());
		}
		double dist = 0;
		for(int i = 0; i<r1.getLats().size(); i++){
			if(dist>r)
				{
					dist = Double.MAX_VALUE;
		//			System.out.println("pruning power============== "+(0.0+i)/r1.getLats().size());
					return dist;
				}
			dist+=(r1.getLats().get(i)-r2.getLats().get(i))*(r1.getLats().get(i)-r2.getLats().get(i))+(r1.getLons().get(i)-r2.getLons().get(i))*(r1.getLons().get(i)-r2.getLons().get(i));

			
		}
	//	System.out.println("!!!!!!!!!!!!!!!!!!NO pruning power!!!!!!!!!!!!!!!!!!!! ");
		return dist;
	}
	public static Double routeHuDist(Route r1, Route r2) {
		if(r1.getLats().size()!=r2.getLats().size())
		{
			throw new IllegalArgumentException("r1.size != r2.size    ->"+r1.getLats().size()+ " : " +r2.getLats().size());
		}
		double dist = 0;
		for(int i = 0; i<r1.getLats().size(); i++){
			double pointDist = Tools.pointEuDist(r1.getLats().get(i), r1.getLons().get(i), r2.getLats().get(i), r2.getLons().get(i));
	//		System.out.println("PointeuDist = "+pointDist);
			dist = dist+ pointDist;
		}
		
		return dist/r1.getLats().size();
	}
	public static int[] parseTrajId(String s) {
		int[] subTraj = new int[3];
		subTraj[0] = Integer.parseInt(s.substring(1, s.indexOf("S")));
		subTraj[1] = Integer.parseInt(s.substring(s.indexOf("S")+1, s.indexOf("L")));
		//subTraj[2] = subTraj[1]+Integer.parseInt(s.substring(s.indexOf("L")+1));
		subTraj[2] = Integer.parseInt(s.substring(s.indexOf("L")+1));
		return subTraj;
		
	}
	public static boolean isTrivialMatch(String sub1,String sub2){
		int[] s1 = parseTrajId(sub1);
		int[] s2 = parseTrajId(sub2);
		if(s1[2]!=s2[2])
			throw new IllegalArgumentException("sub1.length !=sub2.length    "+s1[2]+":"+s2[2]);
		//System.out.println("s1 = "+s1[0]+",   sub2 = "+s2[0] );

		if(s1[0]!=s2[0])
			return false;
		else{	
	//		System.out.println("sub1 = "+sub1+",   sub2 = "+sub2 + "Distance = "+ Math.abs(s1[1]-s2[1]));

			if(Math.abs(s1[1]-s2[1])<s1[2])
				return true;
			else
				return false;
		}
	}
	public static Route getSubroute(int traj, int s, int length){
		Route route = new Route();
		for(int i = s; i<s+length; i++){
			
			Double x = SequiturModel.oldtrajX.get(traj).get(i);
			Double y = SequiturModel.oldtrajY.get(traj).get(i);
		
		
		
		route.addLocation(x,y);
		}
		return route;
	}
	public static double RouteSITED(Route route1, Route route2, double r) {
		Route r1, r2;
		
		if(route1.getLats().size()<=route2.getLats().size()){
			r1 = route1;
			r2 = route2;
		}
		else{
			r1 = route2;
			r2 = route1;
		}
		int p = r1.getLats().size();
		int q = r2.getLats().size();
		
		double minDist = Double.MAX_VALUE;
	//	System.out.println(p+" p:q "+q);
		for(int i = 0; i<=q-p; i++){
			Route subR2 = new Route();
			for(int j = i; j<i+p; j++){
				subR2.addLocation(r2.getLats().get(j), r2.getLons().get(j));
			}
	//		SequiturModel.allEuDistanceCalled++;
			double currentDist = routeSqrEuDist(r1,subR2,r);
		//	System.out.println("currentdist = "+currentDist);
			
			if(currentDist<minDist){
				
				minDist = currentDist;
			}
		}
	//	return minDist*q/p;  //penalty
		return minDist;
	}
	

	

	public static Route getRoute(int i) {
		Route route = new Route();
		for(int j = 0; j<SequiturModel.oldtrajX.get(i).size(); j++){
			double x = SequiturModel.oldtrajX.get(i).get(j);
			double y = SequiturModel.oldtrajY.get(i).get(j);
			route.addLocation(x, y);
		}
		return route;
	}
	
	
}
