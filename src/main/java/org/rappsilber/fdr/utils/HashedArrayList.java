/*
 * Copyright 2015 Lutz Fischer <lfischer at staffmail.ed.ac.uk>.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.rappsilber.fdr.utils;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;

/**
 *
 * @author lfischer
 */
public class HashedArrayList<T> extends ArrayList<T> {
    private HashMap<T,T> hash;

    public HashedArrayList() {
        hash = new HashMap<T,T>();
    }

    public HashedArrayList(int initialCapacity) {
        super(initialCapacity);
        hash = new HashMap<T,T>(initialCapacity);
    }

    public HashedArrayList(Collection<? extends T> c) {
        super(c);
        hash = new HashMap<T,T>(c.size());
        for (T e: c) {
            hash.put(e,e);
        }
    }
    
    
    
    public boolean add(T e) {
        hash.put(e,e);
        return super.add(e);
    }
    
    
    public void add(int index, T e) {
        hash.put(e,e);
        super.add(index, e);
    }
    
    public boolean addAll(Collection<? extends T> c) {
        boolean ret =super.addAll(c);
        for (T e: c) {
            hash.put(e,e);
        }
        return ret;
    } 

    public T set(int index, T e) {
        hash.put(e,e);
        T ret = super.set(index, e);
        hash.remove(ret);
        return ret;
    }
    
    public T remove(int index) {
        T ret = super.remove(index);
        hash.remove(ret);
        return ret;
    }
    
    public boolean remove(Object e) {
        boolean ret = super.remove(e);
        hash.remove(e);
        return ret;
    }
    
    public boolean removeAll(Collection c) {
        boolean ret = super.removeAll(c);
        return ret && hash.keySet().removeAll(c);
    }

    public boolean retainAll(Collection c) {
        boolean ret = super.retainAll(c);
        return ret && hash.keySet().retainAll(c);
    }
 

    public boolean contains(Object e) {
        return hash.containsKey(e);
    }

    public T get(T e) {
        return hash.get(e);
    }
    
    
    public void clear() {
        super.clear();
        hash.clear();
    }
}
