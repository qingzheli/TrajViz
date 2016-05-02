package edu.gmu.trajviz.logic;

public class Location {
	double x;
	double y;
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
	
}
