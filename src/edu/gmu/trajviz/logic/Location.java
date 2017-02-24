package edu.gmu.trajviz.logic;

public class Location {
	private double x;
	private double y;
	public Location(double lat, double lon)
	{
		x = lat;
		y = lon;
	}

	public Location(float lat, float lon)
	{
		x = (double)lat;
		y = (double)lon;
	}
	public double getX(){
		return x;
	}
	public double getY(){
		return y;
	}
	
}
