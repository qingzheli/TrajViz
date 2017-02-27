package edu.gmu.trajviz.logic;

public class UserSession {
	  public static final double DEFAULT_MAX_ERROR_STEP = 5;
	  public static final int DEFAULT_ALPHABET_SIZE = 100;
	  public static final int DEFAULT_MINIMUM_BLOCKS = 3;
	  public static final int DEFALULT_RESAMPLING_RATE = 1000;
	  private double maxErrorSteps;
	  private int alphabetSize;
	  private int minBlocks;
	  private int resamplingRate;
  public UserSession(){
	  super();
	  this.maxErrorSteps = DEFAULT_MAX_ERROR_STEP;
	  this.alphabetSize = DEFAULT_ALPHABET_SIZE;
	  this.minBlocks = DEFAULT_MINIMUM_BLOCKS;
	  this.resamplingRate = DEFALULT_RESAMPLING_RATE;
	 // log("maxErrorSteps: " + maxErrorSteps + ", blockSize: "
	 //           + blockSize + ", Minimal Continous Blocks: " + minBlocks+", resamplingRate:"+resamplingRate);
	  
  }
  public double maxErrorSteps(){
	  return this.maxErrorSteps;
  }
  public int getBlockSize(){
	  return this.alphabetSize;
  }
  public int getMinBlocks(){
	  return this.minBlocks;
  }
  public int getResamplingRate(){
	  return this.resamplingRate;
  }
  
  public void setResamplingRate(int count){
	  this.resamplingRate = count;
  }
  
  public void setMaxErrorSteps(double minLink){
	  this.maxErrorSteps = minLink;
  }
  public void setAlphabet(int alphabet){
	  this.alphabetSize = alphabet;
  }
  public void setMinBlocks(int minBlocks){
	  this.minBlocks = minBlocks;
  }
}
