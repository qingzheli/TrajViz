package edu.gmu.trajviz.logic;

import java.util.ArrayList;
import java.util.HashSet;

import edu.gmu.trajviz.gi.GrammarRuleRecord;
//-qz
public class Cluster {
public String name;
public HashSet<GrammarRuleRecord> ruleSet;
//public ArrayList<RuleInterval> intervals;
public Cluster(String name){
	this.name = name;
	ruleSet = new HashSet<GrammarRuleRecord>();
	//intervals = new ArrayList<RuleInterval>();	
}
public void addRule(GrammarRuleRecord rule){
	ruleSet.add(rule);
//	intervals.addAll(rule.getR0Intervals());
}
public HashSet<GrammarRuleRecord> getRules(){
	return ruleSet;
}
/*
public ArrayList<RuleInterval> getIntervals(){
	return intervals;
}
*/
}
