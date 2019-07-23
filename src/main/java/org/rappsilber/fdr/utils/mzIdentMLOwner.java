/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package org.rappsilber.fdr.utils;

/**
 *
 * @author lfischer
 */
public class mzIdentMLOwner {
    String first;
    String last;
    String email;
    String org;
    String adress;

    public mzIdentMLOwner(String first, String last, String email, String org, String adress) {
        this.first = first;
        this.last = last;
        this.email = email;
        this.org = org;
        this.adress = adress;
    }
    
}
