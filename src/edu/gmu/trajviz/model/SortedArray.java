package edu.gmu.trajviz.model;

import java.util.ArrayList;
import java.util.Collections;

public class SortedArray<T> extends ArrayList<T> {

    @SuppressWarnings("unchecked")
    public void insert(T value) {
        add(value);
        Comparable<T> cmp = (Comparable<T>) value;
        for (int i = size()-1; i > 0 && cmp.compareTo(get(i-1)) < 0; i--)
            Collections.swap(this, i, i-1);
    }
}