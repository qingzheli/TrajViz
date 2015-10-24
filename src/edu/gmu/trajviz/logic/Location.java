package edu.gmu.trajviz.logic;

public class Location {
	double latitude;
	double longitude;
	public Location(double lat, double lon)
	{
		latitude = lat;
		longitude = lon;
	}

	public Location(float lat, float lon)
	{
		latitude = (double)lat;
		longitude = (double)lon;
	}
}
