package edu.gmu.trajviz.logic;

import java.util.ArrayList;
import java.util.List;

import edu.gmu.trajviz.model.SequiturModel;

public class Route{
	private ArrayList<Double> lat;
	private ArrayList<Double> lon;
	private Location startLocation,endLocation;
	public Route(){
		lat = new ArrayList<Double>();
		lon = new ArrayList<Double>();
		}
	public ArrayList<Double> getLats(){
		return lat;
	}
	public ArrayList<Double> getLons(){
		return lon;
	}
	public Route(List<Double> latitudes, List<Double> longitudes)
	{
		this();
		lat.addAll(latitudes);
		lon.addAll(longitudes);
		startLocation = new Location(lat.get(0),lon.get(0));
		endLocation = new Location(lat.get(lat.size()-1),lon.get(lon.size()-1));
	}
	public Block getStartBlock(){
		return SequiturModel.blocks.findBlockByLocation(startLocation);
	}
	public Block getEndBlock(){
		return SequiturModel.blocks.findBlockByLocation(endLocation);
	}
	public void addLocation(double latitude, double longitude){
		lat.add(latitude);
		lon.add(longitude);
		endLocation = new Location(latitude,longitude);
	}
	public void print(){
		for (int i = 0; i<lat.size();i++){
			System.out.println(lat.get(i)+", "+lon.get(i));
		}
	}
	@Override
	public String toString(){
		return lat.toString()+"\n"+lon.toString();
	}
}