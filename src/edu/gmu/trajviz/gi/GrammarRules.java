package edu.gmu.trajviz.gi;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.Map;
import java.util.SortedMap;
import java.util.SortedSet;
import java.util.TreeMap;

public class GrammarRules implements Iterable<GrammarRuleRecord> {

  private SortedMap<Integer, GrammarRuleRecord> rules;
  private GrammarRuleRecord[] sortedRulesByLength;

  public GrammarRules() {
    super();
    this.rules = new TreeMap<Integer, GrammarRuleRecord>();
  }

  public void addRule(GrammarRuleRecord arrRule) {
    int key = arrRule.getRuleNumber();
    this.rules.put(key, arrRule);
  }
/*
 * -qz override
 */
  public void addRule(GrammarRuleRecord arrRule, int ruleNumber) {
	    int key = ruleNumber;
	    this.rules.put(key, arrRule);
	  }
  
  public GrammarRuleRecord getRuleRecord(Integer ruleIdx) {
    return this.rules.get(ruleIdx);
  }

  @Override
  public Iterator<GrammarRuleRecord> iterator() {
    return rules.values().iterator();
  }

  public GrammarRuleRecord get(Integer ruleIndex) {
    return rules.get(ruleIndex);
  }

  public int size() {
    return this.rules.size();
  }
  public void sortByLength() {
	  sortedRulesByLength = new GrammarRuleRecord[rules.size()];
	  
	  rules.values().toArray(sortedRulesByLength);
	  Arrays.sort(sortedRulesByLength);
  }
  public GrammarRuleRecord[] getSortedRulesByLength(){
	  return sortedRulesByLength;
  }



  
}
