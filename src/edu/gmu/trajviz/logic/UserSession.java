package edu.gmu.trajviz.logic;

public class UserSession {
	  public static final int DEFAULT_PAA_SIZE = 1;
	  public static final int DEFAULT_ALPHABET_SIZE = 100;
	  public static final int DEFAULT_MINIMUM_BLOCKS = 5;
	  public static final int DEFALULT_NOISE_POINT_THRESHOLD =1;
	  private int paaSize;
	  private int alphabetSize;
	  private int minBlocks;
	  private int noisePointsThreshold;
  public UserSession(){
	  super();
	  this.paaSize = DEFAULT_PAA_SIZE;
	  this.alphabetSize = DEFAULT_ALPHABET_SIZE;
	  this.minBlocks = DEFAULT_MINIMUM_BLOCKS;
	  this.noisePointsThreshold = DEFALULT_NOISE_POINT_THRESHOLD;
	  
  }
  public int getPAA(){
	  return this.paaSize;
  }
  public int getAlphabet(){
	  return this.alphabetSize;
  }
  public int getMinBlocks(){
	  return this.minBlocks;
  }
  public int getNoisePointThreshold(){
	  return this.noisePointsThreshold;
  }
  
  public void setNoisePointThreshold(int count){
	  this.noisePointsThreshold = count;
  }
  
  public void setPAA(int paaSize){
	  this.paaSize = paaSize;
  }
  public void setAlphabet(int alphabet){
	  this.alphabetSize = alphabet;
  }
  public void setMinBlocks(int minBlocks){
	  this.minBlocks = minBlocks;
  }
}
