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
package org.rappsilber.fdr.result;

import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import org.rappsilber.fdr.utils.FDRSelfAdd;
import org.rappsilber.fdr.utils.HashedArrayList;

/**
 *
 * @author lfischer
 */
public class FDRResultLevel<T extends FDRSelfAdd> extends HashMap<Integer,SubGroupFdrInfo<T>> implements Iterable<T> {
    public boolean isDirectional = false;
    private int within;
    private int between;
    private int linear;

    public int getInputCount() {
        int c=0;
        for (SubGroupFdrInfo<T> g: this.values())
            c+=g.inputCount;
        return c;
    }

    public int getResultCount() {
        int c=0;
        for (SubGroupFdrInfo<T> g: this.values())
            c+=g.filteredResult.size();
        return c;
    }

    public double getLowerFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: this.values()) {
            fdr+=g.results.size()*g.lowerFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public double getHigherFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: this.values()) {
            fdr+=g.results.size()*g.higherFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public double geFirstParsingFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: this.values()) {
            fdr+=g.results.size()*g.firstPassingFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public double getTargetFDR() {
        int c=0;
        double fdr = 0d;
        for (SubGroupFdrInfo<T> g: this.values()) {
            fdr+=g.results.size()*g.targteFDR;
            c+=g.results.size();
        }
        return fdr/c;
    }

    public Iterator<T> iterator() {
        final Iterator<SubGroupFdrInfo<T>> gi = this.values().iterator();

        return new Iterator<T>() {
            SubGroupFdrInfo<T> ng = null;
            Iterator<T> ngi = null;
            public boolean hasNext() {
                if (ngi!= null  && ngi.hasNext())
                    return true;

                while (gi.hasNext()) {
                    ng=gi.next();
                    ngi = ng.iterator();
                    if (ngi.hasNext())
                        return true;
                }
                return false;
            }

            public T next() {
                return ngi.next();
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
        };

    }
    public Iterable<T> filteredRsults() {
        return new Iterable<T>() {

            @Override
            public Iterator<T> iterator() {
                return filteredIterator();
            }
            
        };
    }

    public Iterator<T> filteredIterator() {
        final Iterator<SubGroupFdrInfo<T>> gi = this.values().iterator();

        return new Iterator<T>() {
            SubGroupFdrInfo<T> ng = null;
            Iterator<T> ngi = null;
            public boolean hasNext() {
                if (ngi!= null  && ngi.hasNext())
                    return true;

                while (gi.hasNext()) {
                    ng=gi.next();
                    ngi = ng.filteredResult.iterator();
                    if (ngi.hasNext())
                        return true;
                }
                return false;
            }

            public T next() {
                return ngi.next();
            }

            public void remove() {
                throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
            }
        };

    }
    
    public boolean filteredContains(T e)  {
        for (SubGroupFdrInfo<T> g : this.values()) {
            if (g.filteredContains(e))
                return true;
        }
        return false;
    }

    public T filteredGet(T e)  {
        for (SubGroupFdrInfo<T> g : this.values()) {
            T ret = g.filteredGet(e);
            if (ret!=null) {
                return ret;
            }
        }
        return null;
    }    
    
    public void retainAll(Collection<T> k) {

        for (SubGroupFdrInfo<T> g: values()) {
            g.filteredResult = new HashedArrayList<T>(g.filteredResult);
            g.filteredResult.retainAll(k);
        }
    }

    /**
     * @return the within
     */
    public int getWithin() {
        return within;
    }

    /**
     * @param within the within to set
     */
    public void setWithin(int within) {
        this.within = within;
    }

    /**
     * @return the between
     */
    public int getBetween() {
        return between;
    }

    /**
     * @param between the between to set
     */
    public void setBetween(int between) {
        this.between = between;
    }

    /**
     * @return the linear
     */
    public int getLinear() {
        return linear;
    }

    /**
     * @param linear the linear to set
     */
    public void setLinear(int linear) {
        this.linear = linear;
    }



}
